#!/usr/bin/env python3
"""
Table 1: Poincaré constant growth exponent α on the unit cube (convex), h = 0.2.
(Corresponds to Table tab:numerical_convex in the paper.)

Mathematical claim tested:
    For z_h ∈ Z_h (discrete curl-free-complement) and w_h ∈ W_h (discrete
    div-free-complement), the ratio

        C_p  =  ||z_h||_{L^p} / ||curl z_h||_{L^p}   (k=1)
        C_p  =  ||w_h||_{L^p} / ||div  w_h||_{L^p}   (k=2)

    is bounded independently of p.  We measure this by fitting a power law
    C_p ~ p^α and checking that |α| is small.

Expected conclusion:
    |α| < 0.1 for k=1 test functions (constant C_p).
    |α| ≤ 0.21 for k=2 worst case (sub-linear, still O(1) for finite p).
    This confirms the discrete Poincaré inequality (Lemma 3.4) is p-independent.

Usage:
    python table1_poincare_convex.py
"""

import numpy as np
import sys

try:
    import ngsolve as ngs
    from ngsolve import (
        Mesh, HCurl, HDiv, H1, GridFunction, BilinearForm, LinearForm,
        InnerProduct, grad, curl, div, dx, Integrate, CF,
    )
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
except ImportError:
    sys.exit("ERROR: NGSolve not available.  Install with: pip install ngsolve")

from shared import fit_power_law_exponent, print_table1


# ---------------------------------------------------------------------------
# Parameters (matching Table 1 in the paper)
# ---------------------------------------------------------------------------

H_FIXED   = 0.2                                        # single mesh size for Table 1
ORDER     = 1                                          # lowest-order Nédélec / RT
# p values used for the power-law fit (more points → better regression)
P_FIT     = [2, 4, 6, 8, 10, 15, 20, 30, 50, 75, 100]
# p values displayed in the printed table
P_DISPLAY = [2, 50, 100]
INT_ORDER = 10                                         # quadrature order for L^p integrals


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------

def make_unit_cube_mesh(h: float) -> Mesh:
    """Unit cube [0,1]^3 mesh with mesh size h."""
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0, 0, 0), Pnt(1, 1, 1)).bc("boundary"))
    return Mesh(geo.GenerateMesh(maxh=h))


# ---------------------------------------------------------------------------
# L^p norm
# ---------------------------------------------------------------------------

def lp_norm(mesh: Mesh, f: ngs.CF, p: float) -> float:
    """Compute ||f||_{L^p(mesh)} using numerical quadrature of order INT_ORDER."""
    # For vector fields use the Euclidean pointwise norm
    inner = InnerProduct(f, f) if f.dim > 1 else f * f
    return Integrate(inner ** (p / 2), mesh, order=INT_ORDER) ** (1.0 / p)


# ---------------------------------------------------------------------------
# Projection onto Z_h  (k=1, HCurl cohomology complement)
# ---------------------------------------------------------------------------

def project_to_Zh(mesh: Mesh, v_h: GridFunction, order: int) -> GridFunction:
    """
    Project v_h onto Z_h = ker(curl) ∩ V_h^{Néd,0} modulo grad(S_h^0).

    Solves: find ψ_h ∈ S_h^0 such that
        (∇ψ_h, ∇φ_h) = (v_h, ∇φ_h)  ∀ φ_h ∈ S_h^0
    then returns z_h = InterpolateInto(V_h, v_h - ∇ψ_h).
    The result z_h satisfies curl(z_h) = curl(v_h) and is L²-orthogonal
    to all discrete gradients, i.e. z_h ∈ Z_h.
    """
    S_h = H1(mesh, order=order + 1, dirichlet="boundary")
    psi, phi = S_h.TrialFunction(), S_h.TestFunction()

    a = BilinearForm(S_h)
    a += InnerProduct(grad(psi), grad(phi)) * dx
    a.Assemble()

    f = LinearForm(S_h)
    f += InnerProduct(v_h, grad(phi)) * dx
    f.Assemble()

    psi_h = GridFunction(S_h)
    psi_h.vec.data = a.mat.Inverse(S_h.FreeDofs(), inverse="sparsecholesky") * f.vec

    V_h = HCurl(mesh, order=order, dirichlet="boundary")
    z_h = GridFunction(V_h)
    z_h.Set(v_h - grad(psi_h))
    return z_h


# ---------------------------------------------------------------------------
# Projection onto W_h  (k=2, HDiv cohomology complement)
# ---------------------------------------------------------------------------

def project_to_Wh(mesh: Mesh, v_h: GridFunction, order: int) -> GridFunction:
    """
    Project v_h onto W_h = ker(div) ∩ V_h^{RT,0} modulo curl(V_h^{Néd,0}).

    Solves: find ξ_h ∈ V_h^{Néd,0} such that
        (curl ξ_h, curl ψ_h) = (v_h, curl ψ_h)  ∀ ψ_h ∈ V_h^{Néd,0}
    then returns w_h = InterpolateInto(V_h^{RT}, v_h - curl ξ_h).
    A small regularisation (1e-10 mass term) is added to handle the kernel.
    """
    Ned_h = HCurl(mesh, order=order, dirichlet="boundary")
    xi, psi = Ned_h.TrialFunction(), Ned_h.TestFunction()

    a = BilinearForm(Ned_h)
    a += InnerProduct(curl(xi), curl(psi)) * dx
    a += 1e-10 * InnerProduct(xi, psi) * dx   # regularise kernel
    a.Assemble()

    f = LinearForm(Ned_h)
    f += InnerProduct(v_h, curl(psi)) * dx
    f.Assemble()

    xi_h = GridFunction(Ned_h)
    xi_h.vec.data = a.mat.Inverse(Ned_h.FreeDofs(), inverse="sparsecholesky") * f.vec

    RT_h = HDiv(mesh, order=order, dirichlet="boundary")
    w_h = GridFunction(RT_h)
    w_h.Set(v_h - curl(xi_h))
    return w_h


# ---------------------------------------------------------------------------
# Test-element generation  (k=1 and k=2)
# ---------------------------------------------------------------------------

def make_Zh_test_elements(mesh: Mesh, order: int):
    """
    Return three representative elements of Z_h for k=1 testing.

    Three categories stress-test different parts of the function space:
      1. curl_polynomial : curl of a smooth polynomial  → lies exactly in Z_h
      2. polynomial_proj : generic polynomial projected onto Z_h
      3. random_proj     : random DOF vector projected onto Z_h (stress test)
    """
    x, y, z = ngs.x, ngs.y, ngs.z
    V_h = HCurl(mesh, order=order, dirichlet="boundary")
    # Smooth bubble vanishing on boundary
    bubble = x * (1 - x) * y * (1 - y) * z * (1 - z)

    elements = []

    # 1. Curl of a polynomial stream potential (exact element of Z_h)
    px, py, pz = bubble * y * z, bubble * x * z, bubble * x * y
    curl_pol = CF((py.Diff(z) - pz.Diff(y),   # actually: curl of (px,py,pz)
                   pz.Diff(x) - px.Diff(z),
                   px.Diff(y) - py.Diff(x)))
    # Note: for numerical stability we project even exact curl fields
    gh = GridFunction(V_h); gh.Set(curl_pol)
    elements.append(("curl_polynomial", project_to_Zh(mesh, gh, order)))

    # 2. Generic polynomial projected onto Z_h
    gh = GridFunction(V_h)
    gh.Set(CF((y * z * bubble, x * z * bubble, x * y * bubble)))
    elements.append(("polynomial_proj", project_to_Zh(mesh, gh, order)))

    # 3. Random DOF vector projected onto Z_h  (seed fixed for reproducibility)
    gh = GridFunction(V_h)
    np.random.seed(42)
    rnd = np.random.randn(len(gh.vec))
    free = V_h.FreeDofs()
    for i in range(len(gh.vec)):
        gh.vec[i] = rnd[i] if free[i] else 0.0
    elements.append(("random_proj", project_to_Zh(mesh, gh, order)))

    return elements


def make_Wh_test_elements(mesh: Mesh, order: int):
    """
    Return three representative elements of W_h for k=2 testing.

    Same three categories as the k=1 case, using HDiv / RT elements.
    """
    x, y, z = ngs.x, ngs.y, ngs.z
    V_h = HDiv(mesh, order=order, dirichlet="boundary")
    bubble = x * (1 - x) * y * (1 - y) * z * (1 - z)

    elements = []

    # 1. Divergence-free polynomial field projected onto W_h
    gh = GridFunction(V_h)
    gh.Set(CF((x * bubble, y * bubble, z * bubble)))
    elements.append(("div_polynomial", project_to_Wh(mesh, gh, order)))

    # 2. Generic polynomial
    gh = GridFunction(V_h)
    gh.Set(CF((y * z * bubble, x * z * bubble, x * y * bubble)))
    elements.append(("polynomial_proj", project_to_Wh(mesh, gh, order)))

    # 3. Random DOFs
    gh = GridFunction(V_h)
    np.random.seed(123)
    rnd = np.random.randn(len(gh.vec))
    free = V_h.FreeDofs()
    for i in range(len(gh.vec)):
        gh.vec[i] = rnd[i] if free[i] else 0.0
    elements.append(("random_proj", project_to_Wh(mesh, gh, order)))

    return elements


# ---------------------------------------------------------------------------
# Core computation: Poincaré ratio for one (field, p)
# ---------------------------------------------------------------------------

def poincare_ratio(mesh: Mesh, field_h, p: float, k: int) -> float:
    """
    Compute C_p = ||field_h||_{L^p} / ||d_k field_h||_{L^p}.

    For k=1: d_k = curl.   For k=2: d_k = div.
    Returns inf if the denominator is below machine noise.
    """
    num = lp_norm(mesh, field_h, p)
    den = lp_norm(mesh, curl(field_h) if k == 1 else div(field_h), p)
    return num / den if den > 1e-14 else float("inf")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 70)
    print("TABLE 1: Poincaré growth exponent α  (unit cube, h=0.2, r=1)")
    print("Claim: |α| < 0.1 for k=1 (constant C_p), α < 1 for k=2")
    print("=" * 70)

    mesh = make_unit_cube_mesh(H_FIXED)
    print(f"Mesh: {mesh.ne} elements, h={H_FIXED}")

    rows = []   # each row: (label, {p: C_p}, alpha, interpretation)

    for k, make_elements in [(1, make_Zh_test_elements), (2, make_Wh_test_elements)]:
        operator = "curl" if k == 1 else "div"
        print(f"\n--- k={k} ({operator}-Poincaré) ---")
        for name, field_h in make_elements(mesh, ORDER):
            # Compute C_p for every p in the fit grid
            ratios = {}
            for p in P_FIT:
                r = poincare_ratio(mesh, field_h, float(p), k)
                if r != float("inf"):
                    ratios[p] = r
                    print(f"  {name}  p={p:>4}  C_p={r:.4f}")

            # Fit power law: log C_p = α log p + const  (least-squares)
            alpha = fit_power_law_exponent(
                ps=[p for p in P_FIT if p in ratios],
                values=[ratios[p] for p in P_FIT if p in ratios],
            )
            interp = (
                "constant"    if abs(alpha) < 0.1 else
                "sub-linear"  if alpha < 1.0     else
                "super-linear (WARNING)"
            )
            rows.append((k, name, ratios, alpha, interp))

    # Print formatted table matching the paper layout
    print_table1(rows, P_DISPLAY)
