# Numerical Experiments — Nitsche Hodge-Laplace Optimal Error Estimates

This directory contains the Python scripts that reproduce the four numerical
tables in the paper:

> **Optimal error estimates for the Nitsche Hodge-Laplace method**
> (submitted, 2025)

The experiments confirm the discrete Poincaré inequalities, Galerkin L^p
stability properties, and convergence rate improvements stated in the paper.

---

## Requirements

| Software | Minimum version |
|----------|----------------|
| Python   | 3.10            |
| NGSolve  | 6.2             |
| NumPy    | 1.24            |
| SciPy    | 1.10            |

---

## Installation

Install all dependencies from the project root:

```bash
pip install .
```

Or directly into an existing environment:

```bash
pip install ngsolve numpy scipy
```

NGSolve bundles Netgen (the mesh generator); no separate Netgen installation
is needed.

---

## Running the Experiments

### All tables at once

```bash
cd numericalExperiments
python run_all.py
```

Results are printed to stdout in ASCII table format matching the paper.
Expected runtime: **30–90 minutes** on a modern workstation (dominated by
the inverse-iteration computation in Table 2 and the four-mesh sweep in
Table 3).

### Individual tables

```bash
python table1_poincare_convex.py
python table2_poincare_h_independence.py
python table3_galerkin_lp_stability.py
python table4_poincare_nonconvex.py
```

---

## File Overview

| File | Description |
|------|-------------|
| `shared.py` | Shared utilities: mesh creation, L^p norm computation, power-law fitting, table formatters |
| `run_all.py` | Runner: executes all four table scripts in sequence |
| `table1_poincare_convex.py` | Table 1 |
| `table2_poincare_h_independence.py` | Table 2 |
| `table3_galerkin_lp_stability.py` | Table 3 |
| `table4_poincare_nonconvex.py` | Table 4 |

---

## What Each Table Shows

**Table 1** — Poincaré growth exponent alpha on the unit cube (h = 0.2).
Fits C_p ~ p^alpha for both k=1 (curl-Poincaré) and k=2 (div-Poincaré)
test functions. Confirms |alpha| < 0.1 for k=1, validating that C_p is
essentially constant in p (Lemma 3.4).

**Table 2** — h-independence of the discrete Poincaré constant (k=1).
Evaluates C_p(h) across six mesh sizes h in {0.30, 0.20, 0.15, 0.10,
0.075, 0.05} using the p=2 eigenvalue maximiser. Confirms that variation
remains below 6% in the asymptotic regime, demonstrating that the Poincaré
constant is mesh-size-independent.

**Table 3** — Galerkin L^p stability ratio for the curl-curl problem (k=2).
Computes ||curl(psi_h)||_{L^p} / ||v||_{L^p} for four test functions across
mesh sizes h in {0.3, 0.2, 0.15, 0.1} and exponents p in {2,...,100}.
Confirms growth exponents alpha <= 0.001, establishing that the curl-curl
Galerkin projection is L^p-stable with constants independent of both h and p.

**Table 4** — Poincaré growth exponent alpha on non-convex domains.
Compares the unit cube (convex reference), a 3D L-shaped domain (edge
singularity), and the Fichera corner (vertex singularity) at h = 0.2.
All domains show |alpha| < 0.1, providing numerical evidence that the
p-independence result extends beyond convex Lipschitz domains.

---

## Reproducibility Notes

- Random seeds are fixed inside each script (numpy.random.seed) to ensure
  reproducible results for test functions constructed from random DOF vectors.
- All NGSolve solves use the `sparsecholesky` direct solver.
- L^p integrals use quadrature order 10; this is sufficient for polynomial
  degree r=1 fields with p up to 100.
