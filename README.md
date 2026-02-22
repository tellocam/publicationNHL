# Publication Companion: Improved A Priori Error Estimates for the Nitsche Hodge-Laplacian

This repository is the companion to the paper:

> **Improved A Priori Error Estimates for the Nitsche Hodge-Laplacian**
> Camilo Tello Fachin, Wouter Tonnon, Enrico Zampa (2026)

It contains two independent components: Python scripts that reproduce the four
numerical tables in the paper, and a Lean 4 formalization of the main proof
argument chain.

---

## Authors

- **Camilo Tello Fachin** (ZHAW Institute of Computational Physics, aramiko GmbH)
- **Wouter Tonnon** (ETH Zurich)
- **Enrico Zampa** (University of Vienna)

---

## Repository Structure

```
publicationNHL/
├── numericalExperiments/     # Python scripts reproducing Tables 1–4
│   ├── pyproject.toml        # Package metadata and dependencies
│   ├── shared.py             # Shared utilities (mesh, Lp norms, table formatters)
│   ├── run_all.py            # Runner: executes all four table scripts
│   ├── table1_poincare_convex.py
│   ├── table2_poincare_h_independence.py
│   ├── table3_galerkin_lp_stability.py
│   └── table4_poincare_nonconvex.py
└── formalization/            # Lean 4 formalization of Theorem 4.1 and Corollary 4.2
    ├── lakefile.toml         # Lake build configuration
    ├── lean-toolchain        # Lean version pin (leanprover/lean4:v4.27.0)
    ├── NitscheHodgeLaplacian.lean   # Root import
    ├── NitscheHodgeLaplacian/
    │   ├── Main.lean         # Proof chain: axioms → proved theorems → main results
    │   └── Foundations/      # Axiomatized PDE foundations
    │       ├── BVSpace.lean
    │       ├── DifferentialForms.lean
    │       ├── Discretization.lean
    │       ├── Domains.lean
    │       ├── Embedding.lean
    │       ├── LpSpaces.lean
    │       ├── Sobolev.lean
    │       ├── Trace.lean
    │       └── WeakSolutions.lean
    └── scripts/
        ├── generate_import_graph.sh   # Module import graph (requires graphviz)
        └── generate_proof_graph.py    # Proof-chain dependency graph (requires graphviz)
```

---

## Paper Summary

The paper proves improved a priori error estimates for the Nitsche method applied
to the Hodge-Laplacian. The key improvement over prior work (the thesis bound
h^{r-3/2} + h^{l-1}) is to h^{r-1/2} + h^{l-1/2} for k=1 (H(curl)) and
h^{r-1/2} + h^{l-2/3} for k=2 (H(div)). The proof turns on a novel bound on
the critical term from the BNB analysis, combined with Lp stability of the
hybrid projection and a p-independent discrete Poincare inequality.

---

## Numerical Experiments

The experiments confirm the discrete Poincare inequalities, Galerkin Lp
stability properties, and convergence rate improvements stated in the paper.

### Requirements

| Software | Minimum version |
|----------|----------------|
| Python   | 3.10            |
| NGSolve  | 6.2             |
| NumPy    | 1.24            |
| SciPy    | 1.10            |

NGSolve bundles Netgen (the mesh generator); no separate Netgen installation
is needed.

### Installation

From the repository root:

```bash
pip install ./numericalExperiments
```

Or install dependencies directly:

```bash
pip install ngsolve numpy scipy
```

### Running the experiments

Run all four tables in sequence:

```bash
cd numericalExperiments
python run_all.py
```

Expected runtime: **30–90 minutes** on a modern workstation (dominated by
the inverse-iteration computation in Table 2 and the four-mesh sweep in
Table 3). Results are printed to stdout in ASCII table format matching the paper.

Run individual tables:

```bash
python numericalExperiments/table1_poincare_convex.py
python numericalExperiments/table2_poincare_h_independence.py
python numericalExperiments/table3_galerkin_lp_stability.py
python numericalExperiments/table4_poincare_nonconvex.py
```

### Table Quick Reference

| Table | What it shows |
|-------|---------------|
| **Table 1** | Poincare growth exponent alpha on the unit cube (h=0.2). Confirms |alpha| < 0.1 for k=1, validating p-independence of the discrete Poincare constant (Lemma 3.4). |
| **Table 2** | h-independence of the discrete Poincare constant (k=1). Variation below 6% across six mesh sizes, demonstrating mesh-size independence. |
| **Table 3** | Galerkin Lp stability for the curl-curl problem (k=2). Growth exponents alpha <= 0.001 across four mesh sizes and exponents p in {2,...,100}. |
| **Table 4** | Poincare growth exponent on non-convex domains (unit cube, 3D L-shape, Fichera corner). All show |alpha| < 0.1, extending the p-independence result beyond convex domains. |

### Reproducibility Notes

- Random seeds are fixed inside each script (`numpy.random.seed`) to ensure
  reproducible results for test functions constructed from random DOF vectors.
- All NGSolve solves use the `sparsecholesky` direct solver.
- Lp integrals use quadrature order 10, sufficient for polynomial degree r=1
  fields with p up to 100.

---

## Formalization

The `formalization/` directory contains a Lean 4 formalization of the proof
argument chain for Theorem 4.1 and Corollary 4.2. The formalization captures
the logical structure of the proof: the BNB quasi-optimality framework, the
key bound on the critical term (Lemma 3.3), and the Lp approximation bound
(Corollary 3.7).

**Design note**: the formalization uses axiomatized foundations: opaque types
and explicit axioms for Sobolev spaces, differential forms, triangulations, and
FEEC spaces. This encodes the logical dependencies of the proof without fully
constructing all PDE objects within Lean's type theory. Every declaration is
tagged `[Established: ...]` or `[Novel: ...]` to indicate whether it relies on
prior results or constitutes a new contribution of the paper.

### Main results in `Main.lean`

| Declaration | Statement |
|-------------|-----------|
| `critical_term_lp_bound` | criticalTermSup T omega <= C * h^{l-1/2-1/p} * norm (chains Lemma 3.3 + Corollary 3.7) |
| `parametric_error_bound` | X-norm error <= C * (h^{r-1/2} + h^{l-1/2-1/p} * norm) |
| `nitsche_hodge_laplacian_k1` | Theorem 4.1(a): k=1, rate h^{r-1/2} + h^{l-1/2} |
| `nitsche_hodge_laplacian_k2` | Theorem 4.1(b): k=2, rate h^{r-1/2} + h^{l-2/3} (p=6) |
| `nitsche_hodge_laplacian_k2_conditional` | Corollary 4.2: k=2, rate h^{r-1/2} + h^{l-1/2} under discrete Hodge Lp stability |

### Building

Requires **Lean 4** and **elan** (the Lean version manager). See
[https://leanprover.github.io/lean4/doc/setup.html](https://leanprover.github.io/lean4/doc/setup.html).

```bash
# Install elan (if not already installed)
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh

# From the formalization directory:
cd formalization
lake exe cache get    # Download Mathlib build cache (recommended, avoids multi-hour build)
lake build            # Build the project
```

The build uses Lean 4 v4.27.0 (pinned in `lean-toolchain`) and Mathlib v4.27.0.

### Dependency graphs

Two scripts generate visual dependency graphs (require `graphviz`):

```bash
# Proof-chain graph (Main.lean axioms -> proved theorems -> main results)
python formalization/scripts/generate_proof_graph.py
# Output: formalization/graphs/proof_chain.pdf

# Module import graph
bash formalization/scripts/generate_import_graph.sh
```

Install graphviz if needed:

```bash
# Fedora/RHEL
sudo dnf install graphviz

# Ubuntu/Debian
sudo apt install graphviz

# macOS
brew install graphviz
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
