"""
table5_comparison.py
====================
Prints Table 5 from the paper: "Comparison with prior results"
(remark:comparison in 05_error_analysis.tex).

Row-by-row mathematical meaning
--------------------------------
Row 1 – TelloFachin2025, Thm 3.5.3 (2D/3D, Lipschitz, k=1,2):
    The baseline Nitsche Hodge-Laplace estimate.  The factor h^{r-3/2}
    comes from using a suboptimal inverse estimate to bound the critical
    pairing term, losing one full power of h.  The factor h^{l-1}
    similarly loses half a power relative to the approximation order.

Row 2 – TelloFachin2025 Thm 3.5.5 + AFG2012 (2D, Lipschitz, k=1):
    The two-dimensional result that achieves the near-optimal rate
    h^{r-1/2} + h^l (no half-power loss in the l-term) by exploiting
    the stronger 2D elliptic regularity and the Bramble-Xu superapprox-
    imation argument, which is available in 2D for all polynomial degrees.

Row 3 – New result, Thm (a) (3D, convex, k=1, unconditional):
    Main contribution for k=1 in three dimensions.  Convexity gives a
    p-independent Poincare-Friedrichs constant (arXiv:2508.06741), which
    allows taking p -> inf in the L^p stability of the scalar Galerkin
    problem (Schatz-Wahlbin 1982).  Rate: h^{r-1/2} + h^{l-1/2}.

Row 4 – New result, Thm (b) (3D, convex, k=2, unconditional):
    Main contribution for k=2 in three dimensions.  The L^p stability of
    the discrete Hodge decomposition holds unconditionally only up to p=6
    (exponent 2/3 = 1/2 + 1/6 in l - 1/2 - 1/p with p=6).
    Rate: h^{r-1/2} + h^{l-2/3}.

Row 5 – New result, Cor. (conditional, 3D, convex, k=2):
    If the L^p stability of the discrete Hodge decomposition holds with
    h- and p-independent constant for all p in [2, inf) — satisfied e.g.
    by second-kind (Nedelec II / BDM) elements — then p -> inf recovers
    the optimal l-exponent, matching Row 3.
    Rate: h^{r-1/2} + h^{l-1/2}.
"""

def make_table():
    """Return the comparison table as a list of (source, setting, k, error) rows."""
    rows = [
        (
            "TelloFachin2025, Thm 3.5.3",
            "2D/3D, Lipschitz",
            "1, 2",
            "h^{r-3/2} + h^{l-1}",
        ),
        (
            "TelloFachin2025, Thm 3.5.5 + AFG2012",
            "2D, Lipschitz",
            "1",
            "h^{r-1/2} + h^l",
        ),
        (
            "NEW: Thm (a)  [unconditional]",
            "3D, convex",
            "1",
            "h^{r-1/2} + h^{l-1/2}",
        ),
        (
            "NEW: Thm (b)  [unconditional]",
            "3D, convex",
            "2",
            "h^{r-1/2} + h^{l-2/3}",
        ),
        (
            "NEW: Cor.     [cond. on discrete Hodge L^p stab.]",
            "3D, convex",
            "2",
            "h^{r-1/2} + h^{l-1/2}",
        ),
    ]
    return rows


def print_table(rows):
    headers = ("Source", "Setting", "k", "Error bound")
    col_widths = [
        max(len(headers[i]), max(len(r[i]) for r in rows))
        for i in range(4)
    ]

    sep = "+-" + "-+-".join("-" * w for w in col_widths) + "-+"
    fmt = "| " + " | ".join(f"{{:<{w}}}" for w in col_widths) + " |"

    print(sep)
    print(fmt.format(*headers))
    print(sep)
    for row in rows:
        print(fmt.format(*row))
    print(sep)


if __name__ == "__main__":
    print(__doc__)
    print("Table 5 — Comparison with prior results (remark:comparison)\n")
    print_table(make_table())
