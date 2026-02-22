#!/usr/bin/env python3
"""
generate_proof_graph.py

Generates the proof-chain dependency graph for NitscheHodgeLaplacian/Main.lean
and renders it as a PDF using graphviz.

Usage:
    python generate_proof_graph.py [--output PATH]

Output:  ../graphs/proof_chain.pdf   (default)
Requires: graphviz installed (dot command in PATH)

The graph is hardcoded from the known proof structure in Main.lean (T8 source).
Node categories:
  - axiom/opaque inputs  → red  / dashed border
  - established results  → blue / solid
  - novel results        → green / solid
  - proved theorems      → gold  / solid
  - main theorems        → dark green / bold
  - open hypothesis      → orange / dashed
"""

import argparse
import os
import shutil
import subprocess
import sys
import textwrap


# ---------------------------------------------------------------------------
# Graph definition
# ---------------------------------------------------------------------------

def build_dot() -> str:
    """Return a DOT string encoding the Main.lean proof-chain graph."""

    dot = textwrap.dedent("""\
        digraph proof_chain {
            // Graph attributes
            graph [
                label = "NitscheHodgeLaplacian — Proof Chain Dependency Graph\\nMain.lean (Tello Fachin, Tonnon, Zampa 2026)"
                labelloc = "t"
                fontsize = 16
                fontname = "Helvetica"
                rankdir = "BT"
                bgcolor = "white"
                splines = "polyline"
                nodesep = 0.6
                ranksep = 1.0
            ]

            node [fontname = "Helvetica" fontsize = 10 style = "filled"]
            edge [color = "#444444" arrowsize = 0.75]

            // ----------------------------------------------------------------
            // LEGEND (invisible subgraph to fix legend nodes in top-left area)
            // ----------------------------------------------------------------
            subgraph cluster_legend {
                label = "Legend"
                fontsize = 11
                style = "rounded"
                color = "#aaaaaa"

                L_axiom   [label = "axiom / opaque input"  fillcolor = "#ffcccc" color = "#cc0000" style = "filled,dashed"]
                L_estab   [label = "established result"    fillcolor = "#cce0ff" color = "#0044cc" style = "filled"]
                L_novel   [label = "novel result"          fillcolor = "#ccffcc" color = "#007700" style = "filled"]
                L_proved  [label = "proved theorem"        fillcolor = "#fff3b0" color = "#996600" style = "filled"]
                L_main    [label = "main theorem / corollary" fillcolor = "#00aa55" fontcolor = "white" color = "#005522" style = "filled,bold"]
                L_open    [label = "open hypothesis"       fillcolor = "#ffe0b0" color = "#cc6600" style = "filled,dashed"]
            }

            // ----------------------------------------------------------------
            // TIER 0: Axiom / opaque inputs (the raw hypotheses)
            // ----------------------------------------------------------------
            subgraph cluster_inputs {
                label = "Axioms / Opaque Inputs"
                style = "rounded"
                color = "#cccccc"
                fontsize = 12

                XNorm [
                    label = "XNorm\\n[Established: Def. 3.1.7]\\nMesh-dependent product norm"
                    fillcolor = "#ffcccc"
                    color = "#cc0000"
                    style = "filled,dashed"
                    shape = "box"
                ]
                XNorm_nonneg [
                    label = "XNorm_nonneg\\n[Established]\\nNon-negativity of X-norm"
                    fillcolor = "#ffcccc"
                    color = "#cc0000"
                    style = "filled,dashed"
                    shape = "box"
                ]
                criticalTermSup [
                    label = "criticalTermSup\\n[Novel: proof of Thm 4.1]\\nsup_{y_h≠0} |⟨ω−Pω, d τ_h⟩| / ‖y_h‖_X"
                    fillcolor = "#ccffcc"
                    color = "#007700"
                    style = "filled,dashed"
                    shape = "box"
                ]
                lpApproxError [
                    label = "lpApproxError\\n[Novel: Cor. 3.7]\\n‖ω − Pω‖_{Lp}"
                    fillcolor = "#ccffcc"
                    color = "#007700"
                    style = "filled,dashed"
                    shape = "box"
                ]
                DiscreteHodgeLpStable [
                    label = "DiscreteHodgeLpStable\\n[Open: Issue #2]\\nDiscrete Hodge Lp stability\\n(open problem for k=2, p>6)"
                    fillcolor = "#ffe0b0"
                    color = "#cc6600"
                    style = "filled,dashed"
                    shape = "box"
                ]
            }

            // ----------------------------------------------------------------
            // TIER 1: Direct proof-chain axioms
            // ----------------------------------------------------------------
            bnb_with_noncritical [
                label = "bnb_with_noncritical\\n[Established: Lem. 3.5.2 + Thm. 3.3.7 + Thm. 3.5.3]\\nBNB quasi-optimality + non-critical terms"
                fillcolor = "#cce0ff"
                color = "#0044cc"
                style = "filled"
                shape = "box"
            ]
            key_bound [
                label = "key_bound\\n[Novel: Lemma 3.3]\\nAdapts [AFG12, Thm 3.6] to 3D\\ncriticalTermSup ≤ C·h^{−1/2−1/p}·lpApproxError"
                fillcolor = "#ccffcc"
                color = "#007700"
                style = "filled"
                shape = "box"
            ]
            lp_approximation [
                label = "lp_approximation\\n[Novel: Corollary 3.7]\\nLp stability of P + discrete Poincaré\\nlpApproxError ≤ C·h^l·‖ω‖_{W^{l,p}}"
                fillcolor = "#ccffcc"
                color = "#007700"
                style = "filled"
                shape = "box"
            ]
            k1_limit_bound [
                label = "k1_limit_bound\\n[Novel: Thm 4.1(a) limit step]\\np→∞ limit using [SW82, Thm. 5.1]\\nL∞ stability for k=1"
                fillcolor = "#ccffcc"
                color = "#007700"
                style = "filled"
                shape = "box"
            ]
            k2_conditional_lp_stability [
                label = "k2_conditional_lp_stability\\n[Novel: Corollary 4.2 axiom]\\nUnder DiscreteHodgeLpStable:\\nk=2 achieves same rate as k=1"
                fillcolor = "#ccffcc"
                color = "#007700"
                style = "filled"
                shape = "box"
            ]

            // ----------------------------------------------------------------
            // TIER 2: Private helper lemmas
            // ----------------------------------------------------------------
            absorb_constant [
                label = "absorb_constant\\n[helper]\\nAbsorb two constants via max"
                fillcolor = "#f5f5dc"
                color = "#888800"
                style = "filled"
                shape = "box"
                fontsize = 9
            ]
            exponent_k2 [
                label = "exponent_k2\\n[helper]\\nl − 1/2 − 1/6 = l − 2/3"
                fillcolor = "#f5f5dc"
                color = "#888800"
                style = "filled"
                shape = "box"
                fontsize = 9
            ]

            // ----------------------------------------------------------------
            // TIER 3: Proved intermediate theorems
            // ----------------------------------------------------------------
            critical_term_lp_bound [
                label = "critical_term_lp_bound\\n[Proved: chains key_bound + lp_approximation]\\ncriticalTermSup ≤ C·h^{l−1/2−1/p}·‖ω‖_{W^{l,p}}"
                fillcolor = "#fff3b0"
                color = "#996600"
                style = "filled"
                shape = "box"
                penwidth = 2
            ]
            parametric_error_bound [
                label = "parametric_error_bound\\n[Proved: chains bnb_with_noncritical + critical_term_lp_bound]\\n‖error‖_X ≤ C(h^{r−1/2} + h^{l−1/2−1/p}·‖ω‖_{W^{l,p}})"
                fillcolor = "#fff3b0"
                color = "#996600"
                style = "filled"
                shape = "box"
                penwidth = 2
            ]

            // ----------------------------------------------------------------
            // TIER 4: Main results
            // ----------------------------------------------------------------
            nitsche_hodge_laplacian_k1 [
                label = "nitsche_hodge_laplacian_k1\\n[Theorem 4.1(a)]\\nk=1 (H(curl)): rate h^{r−1/2} + h^{l−1/2}"
                fillcolor = "#00aa55"
                fontcolor = "white"
                color = "#005522"
                style = "filled,bold"
                shape = "box"
                penwidth = 3
                fontsize = 11
            ]
            nitsche_hodge_laplacian_k2 [
                label = "nitsche_hodge_laplacian_k2\\n[Theorem 4.1(b)]\\nk=2 (H(div)): rate h^{r−1/2} + h^{l−2/3}"
                fillcolor = "#00aa55"
                fontcolor = "white"
                color = "#005522"
                style = "filled,bold"
                shape = "box"
                penwidth = 3
                fontsize = 11
            ]
            nitsche_hodge_laplacian_k2_conditional [
                label = "nitsche_hodge_laplacian_k2_conditional\\n[Corollary 4.2]\\nk=2 conditional on DiscreteHodgeLpStable:\\nrate h^{r−1/2} + h^{l−1/2}"
                fillcolor = "#007744"
                fontcolor = "white"
                color = "#003322"
                style = "filled,bold"
                shape = "box"
                penwidth = 3
                fontsize = 11
            ]

            // ----------------------------------------------------------------
            // EDGES: proof dependencies
            // ----------------------------------------------------------------

            // Tier 1 axioms depend on opaque inputs
            bnb_with_noncritical -> XNorm
            bnb_with_noncritical -> criticalTermSup

            key_bound -> criticalTermSup
            key_bound -> lpApproxError

            lp_approximation -> lpApproxError

            k1_limit_bound -> XNorm
            k1_limit_bound -> XNorm_nonneg

            k2_conditional_lp_stability -> XNorm
            k2_conditional_lp_stability -> XNorm_nonneg
            k2_conditional_lp_stability -> DiscreteHodgeLpStable

            // Tier 2 helpers
            absorb_constant -> XNorm_nonneg

            // critical_term_lp_bound chains key_bound + lp_approximation
            critical_term_lp_bound -> key_bound
            critical_term_lp_bound -> lp_approximation

            // parametric_error_bound chains bnb + critical_term
            parametric_error_bound -> bnb_with_noncritical
            parametric_error_bound -> critical_term_lp_bound
            parametric_error_bound -> absorb_constant
            parametric_error_bound -> XNorm_nonneg

            // Main theorems
            nitsche_hodge_laplacian_k1 -> k1_limit_bound

            nitsche_hodge_laplacian_k2 -> parametric_error_bound
            nitsche_hodge_laplacian_k2 -> exponent_k2

            nitsche_hodge_laplacian_k2_conditional -> k2_conditional_lp_stability

            // ----------------------------------------------------------------
            // Rank hints for layered layout
            // ----------------------------------------------------------------
            { rank = same; XNorm; XNorm_nonneg; criticalTermSup; lpApproxError; DiscreteHodgeLpStable }
            { rank = same; bnb_with_noncritical; key_bound; lp_approximation; k1_limit_bound; k2_conditional_lp_stability }
            { rank = same; absorb_constant; exponent_k2; critical_term_lp_bound }
            { rank = same; parametric_error_bound }
            { rank = same; nitsche_hodge_laplacian_k1; nitsche_hodge_laplacian_k2; nitsche_hodge_laplacian_k2_conditional }
        }
    """)
    return dot


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate proof-chain dependency graph PDF for NitscheHodgeLaplacian/Main.lean"
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output PDF path (default: ../graphs/proof_chain.pdf relative to this script)"
    )
    parser.add_argument(
        "--dot-only",
        action="store_true",
        help="Write the DOT file and exit without rendering to PDF"
    )
    args = parser.parse_args()

    # Determine paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    graphs_dir = os.path.join(os.path.dirname(script_dir), "graphs")
    os.makedirs(graphs_dir, exist_ok=True)

    dot_path = os.path.join(graphs_dir, "proof_chain.dot")
    pdf_path = args.output if args.output else os.path.join(graphs_dir, "proof_chain.pdf")

    # Write DOT file
    dot_content = build_dot()
    with open(dot_path, "w", encoding="utf-8") as f:
        f.write(dot_content)
    print(f"DOT file written to: {dot_path}")

    if args.dot_only:
        print("--dot-only specified; skipping PDF rendering.")
        return

    # Check for graphviz
    if shutil.which("dot") is None:
        print("ERROR: graphviz is not installed or 'dot' is not in PATH.", file=sys.stderr)
        print("Install it with:", file=sys.stderr)
        print("  Fedora/RHEL:   sudo dnf install graphviz", file=sys.stderr)
        print("  Ubuntu/Debian: sudo apt install graphviz", file=sys.stderr)
        print("  macOS:         brew install graphviz", file=sys.stderr)
        sys.exit(1)

    # Render DOT -> PDF
    print(f"Rendering {dot_path} -> {pdf_path} ...")
    result = subprocess.run(
        ["dot", "-Tpdf", "-o", pdf_path, dot_path],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("ERROR: graphviz dot failed:", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"\nDone. Proof chain graph written to: {pdf_path}")


if __name__ == "__main__":
    main()
