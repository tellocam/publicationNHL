#!/usr/bin/env python3
"""
Table 2: Discrete Poincaré constant C_p across mesh sizes
         (k=1, unit cube, p=2 eigenvalue maximiser).
(Corresponds to Table tab:numerical_hindep in the paper.)

Mathematical claim tested:
    The discrete Poincaré constant

        C_p(h)  =  sup_{z_h ∈ Z_h \\ {0}}  ||z_h||_{L^p} / ||curl z_h||_{L^p}

    is bounded independently of the mesh size h.  We approximate C_p(h) from
    below by evaluating the ratio on the "p=2 eigenvalue maximiser": the
    element z_h* ∈ Z_h that maximises ||z_h||_{L^2} / ||curl z_h||_{L^2},
    found via inverse iteration on the generalised eigenvalue problem
        M z = λ K z   (M = L² mass matrix, K = curl-curl matrix on Z_h).
    The p=2 maximiser gives the tightest available lower bound on C_p.

Expected conclusion:
    Variation in C_p across h ∈ {0.30, 0.20, 0.15, 0.10, 0.075, 0.05}
    is below 6 % in the asymptotic regime (h ≤ 0.2), confirming that
    the discrete Poincaré constant is h-independent (Lemma 3.4).

Usage:
    python table2_poincare_h_independence.py
"""

import numpy as np
import sys

try:
    import ngsolve as ngs
    from ngsolve import (
        Mesh, HCurl, H1, GridFunction, BilinearForm, LinearForm,
        InnerProduct, grad, curl, dx, Integrate,
    )
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
except ImportError:
    sys.exit("ERROR: NGSolve not available.  Install with: pip install ngsolve")

from shared import fit_power_law_exponent, print_table2


# ---------------------------------------------------------------------------
# Parameters (matching Table 2 in the paper)
# ---------------------------------------------------------------------------

H_VALUES  = [0.30, 0.20, 0.15, 0.10, 0.075, 0.05]    # mesh refinement sequence
ORDER     = 1                                           # lowest-order Nédélec
# p values at which C_p is evaluated (display columns + fit grid)
P_VALUES  = [2, 10, 50, 100]
INT_ORDER = 10                                          # quadrature order
INV_ITER  = 20                                          # inverse-iteration steps


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
    """Compute ||f||_{L^p(mesh)} via numerical quadrature."""
    inner = InnerProduct(f, f) if f.dim > 1 else f * f
    return Integrate(inner ** (p / 2), mesh, order=INT_ORDER) ** (1.0 / p)


# ---------------------------------------------------------------------------
# Projection onto Z_h
# ---------------------------------------------------------------------------

def project_to_Zh(mesh: Mesh, v_h: GridFunction, order: int) -> GridFunction:
    """
    Project v_h onto Z_h = (grad S_h^0)^⊥ in V_h^{Néd,0}.

    Solves the Poisson problem for the potential ψ_h and returns
    z_h = v_h - grad ψ_h  interpolated back into V_h^{Néd,0}.
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
# p=2 Poincaré maximiser via inverse iteration
# ---------------------------------------------------------------------------

def compute_p2_maximiser(mesh: Mesh, order: int) -> GridFunction:
    """
    Find z_h* ∈ Z_h that maximises ||z_h||_{L^2} / ||curl z_h||_{L^2}.

    This is the eigenvector corresponding to the smallest eigenvalue of
        K z = λ M z
    restricted to Z_h, where M is the L² mass matrix and K the curl-curl
    stiffness matrix.  Equivalently, it is the maximiser of the Rayleigh
    quotient  ||z||² / ||curl z||²  over Z_h.

    Algorithm (inverse iteration on Z_h):
      1. Start with a random vector in Z_h.
      2. Multiply by M  (map z → M z).
      3. Solve K w = M z  (a small regularisation 1e-10 M is added to K
         to handle its kernel outside Z_h).
      4. Project w back onto Z_h.
      5. Normalise in L².
      6. Repeat INV_ITER times.

    Convergence: at each step the component along the dominant eigenvector
    (smallest eigenvalue of K w.r.t. M on Z_h) grows relative to all others,
    so the iterate converges to z_h*.
    """
    V_h = HCurl(mesh, order=order, dirichlet="boundary")
    u, v = V_h.TrialFunction(), V_h.TestFunction()

    # L² mass matrix on V_h
    mass = BilinearForm(V_h)
    mass += InnerProduct(u, v) * dx
    mass.Assemble()

    # Regularised curl-curl stiffness (invertible on all of V_h)
    K_reg = BilinearForm(V_h)
    K_reg += InnerProduct(curl(u), curl(v)) * dx
    K_reg += 1e-10 * InnerProduct(u, v) * dx
    K_reg.Assemble()
    inv_K = K_reg.mat.Inverse(V_h.FreeDofs(), inverse="sparsecholesky")

    # Initialise with a random vector projected onto Z_h
    v_h = GridFunction(V_h)
    np.random.seed(999)
    rnd = np.random.randn(len(v_h.vec))
    free = V_h.FreeDofs()
    for i in range(len(v_h.vec)):
        v_h.vec[i] = rnd[i] if free[i] else 0.0
    z_h = project_to_Zh(mesh, v_h, order)

    # Inverse iteration loop
    for step in range(INV_ITER):
        # Apply mass matrix
        rhs = z_h.vec.CreateVector()
        rhs.data = mass.mat * z_h.vec

        # Solve K w = M z_h  (moves weight toward small-eigenvalue component)
        w_h = GridFunction(V_h)
        w_h.vec.data = inv_K * rhs

        # Project back to Z_h and normalise
        z_h = project_to_Zh(mesh, w_h, order)
        norm_l2 = lp_norm(mesh, z_h, 2.0)
        if norm_l2 < 1e-14:
            raise RuntimeError(f"Inverse iteration collapsed at step {step}")
        z_h.vec.data = (1.0 / norm_l2) * z_h.vec

    return z_h


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 70)
    print("TABLE 2: h-independence of C_p  (k=1, unit cube, p=2 maximiser)")
    print("Claim: variation < 6% across h ∈ [0.05, 0.30] for h ≤ 0.20")
    print("=" * 70)

    rows = []   # each row: (h, ndof, {p: C_p}, alpha)

    for h in H_VALUES:
        print(f"\n--- h = {h} ---")
        mesh = make_unit_cube_mesh(h)
        V_h  = HCurl(mesh, order=ORDER, dirichlet="boundary")
        ndof = V_h.ndof
        nel  = mesh.ne
        print(f"  {nel} elements, {ndof} DOFs")

        # Compute the p=2 maximiser: the element of Z_h with largest C_2
        print("  Computing p=2 eigenvalue maximiser (inverse iteration)...")
        z_star = compute_p2_maximiser(mesh, ORDER)

        # Evaluate C_p = ||z*||_{L^p} / ||curl z*||_{L^p} for each display p
        ratios = {}
        for p in P_VALUES:
            num = lp_norm(mesh, z_star, float(p))
            den = lp_norm(mesh, curl(z_star), float(p))
            cp  = num / den if den > 1e-14 else float("inf")
            ratios[p] = cp
            print(f"  p={p:>4}: C_p = {cp:.4f}")

        # Fit power law across p values
        valid_ps  = [p for p in P_VALUES if ratios[p] != float("inf")]
        valid_cps = [ratios[p] for p in valid_ps]
        alpha = fit_power_law_exponent(valid_ps, valid_cps)

        rows.append((h, ndof, ratios, alpha))

    # Print formatted table matching the paper layout
    print_table2(rows, P_VALUES)

    # Summary statistics (asymptotic regime: h ≤ 0.20)
    asym_rows = [(h, nd, r, a) for h, nd, r, a in rows if h <= 0.20]
    c2_vals   = [r[2] for _, _, r, _ in asym_rows]
    variation = (max(c2_vals) - min(c2_vals)) / np.mean(c2_vals) * 100
    print(f"\nAsymptotic regime (h ≤ 0.20): C_2 variation = {variation:.1f}%")
    if variation < 6.0:
        print("PASS: C_p is h-independent (< 6% variation)")
    else:
        print(f"WARN: variation {variation:.1f}% exceeds 6% threshold")
