#!/usr/bin/env python3
"""
run_all.py — Run all numerical experiment tables in sequence.

Each table script is imported and its main logic is invoked.
Results are printed to stdout.

Usage:
    python run_all.py

Expected total runtime: 30–90 minutes on a modern workstation.
The bottleneck is Table 2 (inverse iteration) and Table 3 (multi-mesh solves).

Tables
------
  Table 1 — Poincaré growth exponent alpha (unit cube, h=0.2)
  Table 2 — h-independence of the discrete Poincaré constant (k=1)
  Table 3 — Galerkin L^p stability ratio (k=2, curl-curl)
  Table 4 — Poincaré constant on non-convex domains (cube, L-shaped, Fichera)
"""

import importlib
import sys
import time


def run_module(module_name: str, description: str) -> None:
    """Import a table script and execute its __main__ block via runpy."""
    import runpy

    banner = f" {description} "
    print("\n" + "=" * 70)
    print(banner.center(70, "="))
    print("=" * 70)

    t0 = time.perf_counter()
    try:
        runpy.run_module(module_name, run_name="__main__", alter_sys=True)
    except SystemExit as exc:
        # A sys.exit(0) from the script is not an error
        if exc.code not in (None, 0):
            print(f"ERROR: {module_name} exited with code {exc.code}", file=sys.stderr)
            sys.exit(exc.code)
    elapsed = time.perf_counter() - t0
    print(f"\n[{module_name} completed in {elapsed:.1f}s]")


def main() -> None:
    tables = [
        ("table1_poincare_convex",         "Table 1 — Poincaré growth exponent alpha"),
        ("table2_poincare_h_independence", "Table 2 — h-independence of C_p"),
        ("table3_galerkin_lp_stability",   "Table 3 — Galerkin L^p stability (k=2)"),
        ("table4_poincare_nonconvex",      "Table 4 — Poincaré on non-convex domains"),
    ]

    total_start = time.perf_counter()
    for module_name, description in tables:
        run_module(module_name, description)

    total = time.perf_counter() - total_start
    print(f"\n{'=' * 70}")
    print(f"All tables complete.  Total time: {total / 60:.1f} min")
    print("=" * 70)


if __name__ == "__main__":
    main()
