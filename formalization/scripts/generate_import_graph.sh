#!/bin/bash
# generate_import_graph.sh
#
# Generates a module-level import dependency graph for the
# NitscheHodgeLaplacian Lean 4 formalization and renders it as a PDF.
#
# Usage: ./generate_import_graph.sh
# Output: ../graphs/import_graph.pdf
#
# Requires: graphviz (dot command)

set -euo pipefail

# ---------------------------------------------------------------------------
# Check for graphviz
# ---------------------------------------------------------------------------
if ! command -v dot &>/dev/null; then
    echo "ERROR: graphviz is not installed or 'dot' is not in PATH." >&2
    echo "Install it with:" >&2
    echo "  Fedora/RHEL:  sudo dnf install graphviz" >&2
    echo "  Ubuntu/Debian: sudo apt install graphviz" >&2
    echo "  macOS:         brew install graphviz" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Setup output directory
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRAPHS_DIR="$(dirname "$SCRIPT_DIR")/graphs"
mkdir -p "$GRAPHS_DIR"

DOT_FILE="$GRAPHS_DIR/import_graph.dot"
PDF_FILE="$GRAPHS_DIR/import_graph.pdf"

# ---------------------------------------------------------------------------
# Determine module source: try lake exe graph first, fall back to hardcoded
# ---------------------------------------------------------------------------

# Short module name helper (strips "NitscheHodgeLaplacian." prefix for labels)
# We use the full dotted name as the node ID, and a short label.

LEAN_ROOT="$(dirname "$SCRIPT_DIR")"  # parent of scripts/ == repo root

# Try lake exe graph if lake is available and we are inside the repo
USE_LAKE=false
if command -v lake &>/dev/null && [ -f "$LEAN_ROOT/lakefile.toml" ]; then
    echo "lake found; attempting 'lake exe graph' ..."
    # lake exe graph writes to stdout a DOT graph; redirect to file
    if lake --dir "$LEAN_ROOT" exe graph > "$DOT_FILE" 2>/dev/null; then
        USE_LAKE=true
        echo "lake graph succeeded; rendering PDF ..."
    else
        echo "lake exe graph failed or not available; falling back to hardcoded graph."
    fi
fi

# ---------------------------------------------------------------------------
# Hardcoded graph (used when lake is unavailable or fails)
# Derived from T8-lean-merge.md import analysis.
# ---------------------------------------------------------------------------
if [ "$USE_LAKE" = false ]; then
    echo "Generating hardcoded import graph from known module structure ..."

    cat > "$DOT_FILE" <<'EOF'
digraph NitscheHodgeLaplacian_imports {
    // Graph-level attributes
    graph [
        label = "NitscheHodgeLaplacian — Module Import Graph\n(Lean 4 / Mathlib v4.27.0)"
        labelloc = "t"
        fontsize = 16
        fontname = "Helvetica"
        rankdir = "BT"
        bgcolor = "white"
        splines = "ortho"
        nodesep = 0.5
        ranksep = 0.8
    ]

    node [
        shape = "box"
        style = "filled"
        fontname = "Helvetica"
        fontsize = 11
    ]

    edge [
        color = "#555555"
        arrowsize = 0.7
    ]

    // -----------------------------------------------------------------------
    // External dependency (Mathlib)
    // -----------------------------------------------------------------------
    Mathlib [
        label = "Mathlib\n(v4.27.0)"
        fillcolor = "#d0d0d0"
        style = "filled,dashed"
        shape = "ellipse"
        fontsize = 13
    ]

    // -----------------------------------------------------------------------
    // Root module
    // -----------------------------------------------------------------------
    Root [
        label = "NitscheHodgeLaplacian\n(root)"
        fillcolor = "#1a1a2e"
        fontcolor = "white"
        style = "filled,bold"
        fontsize = 13
    ]

    // -----------------------------------------------------------------------
    // Main theorem file
    // -----------------------------------------------------------------------
    Main [
        label = "Main\n(Theorems 4.1, Cor. 4.2)"
        fillcolor = "#16213e"
        fontcolor = "white"
        style = "filled,bold"
        fontsize = 12
    ]

    // -----------------------------------------------------------------------
    // Foundations modules — grouped by dependency tier
    // -----------------------------------------------------------------------

    // Tier 1: no inter-Foundations deps
    Domains [
        label = "Foundations.Domains"
        fillcolor = "#e8f4f8"
    ]
    WeakSolutions [
        label = "Foundations.WeakSolutions"
        fillcolor = "#e8f4f8"
    ]

    // Tier 2: depend on Domains only
    Sobolev [
        label = "Foundations.Sobolev"
        fillcolor = "#cce5ff"
    ]

    // Tier 3: depend on Domains + Sobolev
    LpSpaces [
        label = "Foundations.LpSpaces"
        fillcolor = "#b8daff"
    ]
    Trace [
        label = "Foundations.Trace"
        fillcolor = "#b8daff"
    ]

    // Tier 4: depend on Domains + Sobolev + LpSpaces
    DifferentialForms [
        label = "Foundations.DifferentialForms"
        fillcolor = "#99c9ff"
    ]
    BVSpace [
        label = "Foundations.BVSpace"
        fillcolor = "#99c9ff"
    ]

    // Tier 5: depend on Sobolev + LpSpaces + DifferentialForms + Domains
    Embedding [
        label = "Foundations.Embedding"
        fillcolor = "#70b0ff"
    ]

    // Tier 6: depend on Domains + DifferentialForms
    Discretization [
        label = "Foundations.Discretization"
        fillcolor = "#4499ee"
        fontcolor = "white"
    ]

    // -----------------------------------------------------------------------
    // Edges: Mathlib usage (every module imports Mathlib.Tactic implicitly)
    // Shown only for modules that explicitly list Mathlib imports in T8
    // -----------------------------------------------------------------------
    Domains        -> Mathlib
    Sobolev        -> Mathlib
    LpSpaces       -> Mathlib
    DifferentialForms -> Mathlib
    WeakSolutions  -> Mathlib
    BVSpace        -> Mathlib
    Discretization -> Mathlib
    Main           -> Mathlib  [label = "Analysis.SpecialFunctions.Pow.Real" fontsize = 8]

    // -----------------------------------------------------------------------
    // Edges: inter-Foundations imports
    // -----------------------------------------------------------------------
    Sobolev           -> Domains
    LpSpaces          -> Domains
    LpSpaces          -> Sobolev
    Trace             -> Sobolev
    DifferentialForms -> Domains
    DifferentialForms -> Sobolev
    DifferentialForms -> LpSpaces
    Embedding         -> Domains
    Embedding         -> Sobolev
    Embedding         -> LpSpaces
    Embedding         -> DifferentialForms
    BVSpace           -> Domains
    BVSpace           -> LpSpaces
    BVSpace           -> Sobolev
    Discretization    -> Domains
    Discretization    -> DifferentialForms

    // -----------------------------------------------------------------------
    // Edges: Main and Root
    // -----------------------------------------------------------------------
    Main -> Discretization

    Root -> Domains
    Root -> Sobolev
    Root -> LpSpaces
    Root -> Embedding
    Root -> Trace
    Root -> WeakSolutions
    Root -> DifferentialForms
    Root -> BVSpace
    Root -> Discretization
    Root -> Main

    // -----------------------------------------------------------------------
    // Rank hints for cleaner layout
    // -----------------------------------------------------------------------
    { rank = same; Domains; WeakSolutions }
    { rank = same; Sobolev }
    { rank = same; LpSpaces; Trace }
    { rank = same; DifferentialForms; BVSpace }
    { rank = same; Embedding; Discretization }
    { rank = same; Main }
    { rank = same; Root }
}
EOF

fi  # end USE_LAKE=false block

# ---------------------------------------------------------------------------
# Render DOT → PDF
# ---------------------------------------------------------------------------
echo "Rendering $DOT_FILE → $PDF_FILE ..."
dot -Tpdf -o "$PDF_FILE" "$DOT_FILE"

echo ""
echo "Done. Import graph written to: $PDF_FILE"
