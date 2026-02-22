"""
shared.py — Shared utilities for NHL numerical experiments.

Provides mesh creation, L^p norm computation, and power-law fitting
helpers used across the table scripts (table1 through table4).

NGSolve is required for make_unit_cube_mesh and compute_Lp_norm.
fit_power_law / fit_power_law_exponent / fit_p_growth_exponent
are pure Python + NumPy and have no NGSolve dependency.
"""

from __future__ import annotations

import numpy as np
import sys

try:
    import ngsolve as ngs
    from ngsolve import Integrate, InnerProduct
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
except ImportError:
    sys.exit(
        "ERROR: NGSolve not available.  "
        "Install the numerical extras with:  pip install ngsolve"
    )


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------

def make_unit_cube_mesh(h: float, order: int = 1) -> ngs.Mesh:
    """
    Create a simplicial mesh of the unit cube [0,1]^3.

    Parameters
    ----------
    h : float
        Target mesh size (maxh parameter passed to Netgen).
    order : int
        Mesh curvature order (default 1 = flat faces).

    Returns
    -------
    ngs.Mesh
    """
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0, 0, 0), Pnt(1, 1, 1)).bc("boundary"))
    return ngs.Mesh(geo.GenerateMesh(maxh=h))


# ---------------------------------------------------------------------------
# L^p norm
# ---------------------------------------------------------------------------

def compute_Lp_norm(
    mesh: ngs.Mesh,
    f,
    p: float,
    order_integrate: int = 10,
) -> float:
    """
    Compute ||f||_{L^p(mesh)} = ( integral |f|^p dx )^{1/p}.

    Works for both scalar and vector NGSolve CoefficientFunctions.
    Vector fields use the Euclidean pointwise norm: |f|^2 = f · f.

    Parameters
    ----------
    mesh : ngs.Mesh
    f    : NGSolve CF or GridFunction (scalar or vector)
    p    : float  — exponent (>= 1)
    order_integrate : int — quadrature order passed to Integrate

    Returns
    -------
    float
    """
    f_sq = InnerProduct(f, f) if f.dim > 1 else f * f
    return float(Integrate(f_sq ** (p / 2), mesh, order=order_integrate)) ** (1.0 / p)


# ---------------------------------------------------------------------------
# Power-law fitting
# ---------------------------------------------------------------------------

def fit_power_law(p_values: list[float], c_values: list[float]) -> float:
    """
    Fit C_p ~ A * p^alpha via least-squares on a log-log scale.

    Filters out non-positive and non-finite values before fitting.

    Parameters
    ----------
    p_values : sequence of floats  — exponent values (>0)
    c_values : sequence of floats  — corresponding C_p values (>0)

    Returns
    -------
    float — growth exponent alpha (slope of log C_p vs log p).
             Returns nan if fewer than 2 valid pairs are available.
    """
    pairs = [
        (p, c) for p, c in zip(p_values, c_values)
        if p > 0 and c > 0 and np.isfinite(c)
    ]
    if len(pairs) < 2:
        return float("nan")
    log_p = np.array([np.log(p) for p, _ in pairs])
    log_c = np.array([np.log(c) for _, c in pairs])
    alpha = float(np.polyfit(log_p, log_c, 1)[0])
    return alpha


# Alias used by table1 / table2 scripts (same (ps, values) call signature)
def fit_power_law_exponent(ps: list[float], values: list[float]) -> float:
    """
    Alias of fit_power_law with an alternate name used by table1/table2.

    Returns 0.0 instead of nan when fewer than 2 points are available,
    matching the table1/table2 convention.
    """
    result = fit_power_law(ps, values)
    return result if np.isfinite(result) else 0.0


def fit_p_growth_exponent(
    log_p: list[float],
    log_r: list[float],
) -> float | None:
    """
    Log-log least-squares fit on pre-computed log arrays.

    Used by table3 (T4), which accumulates log(p) and log(ratio) values
    across all mesh sizes before fitting.

    Parameters
    ----------
    log_p : list of log(p) values
    log_r : list of log(ratio) values (same length as log_p)

    Returns
    -------
    float — slope (growth exponent alpha), or None if < 3 points.
    """
    if len(log_p) < 3:
        return None
    lp = np.array(log_p, dtype=float)
    lr = np.array(log_r, dtype=float)
    mask = np.isfinite(lp) & np.isfinite(lr)
    if mask.sum() < 3:
        return None
    A = np.vstack([lp[mask], np.ones(mask.sum())]).T
    alpha, _ = np.linalg.lstsq(A, lr[mask], rcond=None)[0]
    return float(alpha)


# ---------------------------------------------------------------------------
# Table formatting helpers (used by table1 and table2)
# ---------------------------------------------------------------------------

def print_table1(rows: list, p_display: list[int]) -> None:
    """
    Print Table 1 (Poincaré growth exponent alpha) in paper-matching format.

    Parameters
    ----------
    rows : list of (k, name, ratios_dict, alpha, interpretation)
           ratios_dict maps p -> C_p float
    p_display : p values to show as columns (e.g. [2, 50, 100])
    """
    header = f"\n{'Test function':<22}" + "".join(f"  C_{p:<5}" for p in p_display)
    header += f"  {'alpha':>8}  Interpretation"
    sep = "-" * len(header)

    print(f"\n{'TABLE 1 — Poincaré growth exponent alpha (h=0.2, r=1)':^70}")
    print(sep)
    print(header)
    print(sep)

    current_k = None
    for k, name, ratios, alpha, interp in rows:
        if k != current_k:
            label = "k=1 (curl-Poincaré)" if k == 1 else "k=2 (div-Poincaré)"
            print(f"  [{label}]")
            current_k = k
        cp_cols = "".join(
            f"  {ratios[p]:>6.3f}" if p in ratios else "     ---"
            for p in p_display
        )
        alpha_str = f"{alpha:+.3f}"
        print(f"  {name:<22}{cp_cols}  {alpha_str:>8}  {interp}")
    print(sep)


def print_table2(rows: list, p_display: list[int]) -> None:
    """
    Print Table 2 (h-independence of C_p) in paper-matching format.

    Parameters
    ----------
    rows : list of (h, ndof, ratios_dict, alpha)
           ratios_dict maps p -> C_p float
    p_display : p values to show as columns (e.g. [2, 10, 50, 100])
    """
    header = f"\n{'h':>6}  {'ndof':>6}" + "".join(f"  C_{p:<5}" for p in p_display)
    header += f"  {'alpha':>8}"
    sep = "-" * len(header)

    print(f"\n{'TABLE 2 — h-independence of C_p (k=1, p=2 maximiser)':^70}")
    print(sep)
    print(header)
    print(sep)

    for h, ndof, ratios, alpha in rows:
        cp_cols = "".join(
            f"  {ratios[p]:>6.3f}" if p in ratios else "     ---"
            for p in p_display
        )
        alpha_str = f"{alpha:+.3f}"
        print(f"  {h:>5.3f}  {ndof:>6}{cp_cols}  {alpha_str:>8}")
    print(sep)
