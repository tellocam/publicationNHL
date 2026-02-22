import Mathlib.Tactic

/-!
# Bounded Domains

This module formalizes bounded open domains in ℝ^n with Lipschitz boundary
as an opaque type, together with geometric predicates used throughout the
Nitsche Hodge-Laplacian error analysis.

## Mathematical objects

- `BoundedDomain n` — opaque type for a bounded Lipschitz domain in ℝ^n
- `BoundedDomain.IsConvex` — convexity predicate (needed for p-independent
  Poincaré-Friedrichs inequalities, see Costabel-McIntosh and [CLV25])
- `BoundedDomain.IsContractible` — contractibility predicate (needed for
  triviality of de Rham cohomology, i.e. exactness of the continuous complex)

## References

- Adams & Fournier, "Sobolev Spaces" (2nd ed.), Academic Press, 2003.
  §4 for Lipschitz domains and extension operators.
- Grisvard, "Elliptic Problems in Nonsmooth Domains", SIAM, 2011.
  Chapter 1 for domain geometry.
- Costabel, McIntosh [CLV25] for p-independent Poincaré on convex domains.

## Design

All fields are opaque or axiomatized. No concrete construction of ℝ^n domains
is given; the axiom system is consistent with standard PDE theory.
-/

noncomputable section

namespace NitscheHodgeLaplacian.Foundations

/-! ## Domains -/

/-- A bounded open domain in ℝ^n with Lipschitz boundary.
This is an opaque type; properties are accessed through axioms.

Reference: Adams-Fournier §4, Grisvard Chapter 1. -/
opaque BoundedDomain (n : ℕ) : Type

/-- The spatial dimension of a domain. -/
axiom BoundedDomain.dim {n : ℕ} (Ω : BoundedDomain n) : ℕ

@[simp]
axiom BoundedDomain.dim_eq {n : ℕ} (Ω : BoundedDomain n) : Ω.dim = n

/-- A bounded domain is nonempty. -/
axiom BoundedDomain.nonempty {n : ℕ} (Ω : BoundedDomain n) : Nonempty (BoundedDomain n)

/-! ## Domain Geometry Predicates -/

/-- Convexity predicate on a bounded domain.
Needed for p-independent Poincaré-Friedrichs inequalities.

Reference: Costabel-McIntosh; also used in [Tello Fachin 2025, §3.5]. -/
opaque BoundedDomain.IsConvex {n : ℕ} (Ω : BoundedDomain n) : Prop

/-- Topological contractibility predicate on a bounded domain.
Needed for triviality of cohomology (kernel of exterior derivative equals
image of exterior derivative in all degrees k ≥ 1).

Reference: de Rham's theorem; [Arnold-Falk-Winther 2006, §2]. -/
opaque BoundedDomain.IsContractible {n : ℕ} (Ω : BoundedDomain n) : Prop

/-- Convex domains are contractible.

A convex domain is star-shaped (hence contractible) via the straight-line
homotopy to any interior point.

Reference: Standard topology; used implicitly throughout [Tello Fachin 2025]. -/
axiom BoundedDomain.IsConvex.isContractible {n : ℕ} {Ω : BoundedDomain n}
    (hconv : Ω.IsConvex) : Ω.IsContractible

end NitscheHodgeLaplacian.Foundations
