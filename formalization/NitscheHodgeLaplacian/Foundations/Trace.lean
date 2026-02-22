import NitscheHodgeLaplacian.Foundations.Sobolev

/-!
# Trace Operator

This module axiomatizes the trace operator T : W^{1,p}(Ω) → L^p(∂Ω)
and its kernel characterization.

## Mathematical objects

- `TraceSpace p Ω` — L^p(∂Ω), the trace space on the boundary
- `trace_operator` — bounded linear map W^{1,p} → TraceSpace
- `trace_kernel_eq_sobolev_zero` — ker(T) = W_0^{1,p}

## References

- Adams & Fournier, "Sobolev Spaces" (2nd ed.), §4.10–§4.11.
- Grisvard, "Elliptic Problems in Nonsmooth Domains", Chapter 1.
- Evans, "Partial Differential Equations" (2nd ed.), §5.5.

## Design

TraceSpace is an opaque type representing L^p(∂Ω). The kernel characterization
`trace_kernel_eq_sobolev_zero` is the fundamental result connecting Dirichlet
boundary conditions to W_0^{1,p} spaces.
-/

noncomputable section

open scoped NNReal ENNReal

namespace NitscheHodgeLaplacian.Foundations

/-! ## Trace Operator -/

/-- The trace space L^p(∂Ω).
For smooth functions, restriction to the boundary lands in this space.

Reference: Adams-Fournier §4.10; Grisvard Chapter 1. -/
opaque TraceSpace (p : ℝ≥0∞) {n : ℕ} (Ω : BoundedDomain n) : Type

/-- TraceSpace is a normed add comm group. -/
axiom TraceSpace.instNormedAddCommGroup {p : ℝ≥0∞} {n : ℕ}
    {Ω : BoundedDomain n} : NormedAddCommGroup (TraceSpace p Ω)

attribute [instance] TraceSpace.instNormedAddCommGroup

/-- L^p(∂Ω) is a normed space over ℝ. -/
axiom TraceSpace.instNormedSpace {p : ℝ≥0∞} {n : ℕ}
    {Ω : BoundedDomain n} : NormedSpace ℝ (TraceSpace p Ω)

attribute [instance] TraceSpace.instNormedSpace

/-- The trace operator is a bounded linear map W^{1,p}(Ω) → L^p(∂Ω).

Reference: Adams-Fournier Theorem 4.10; Evans §5.5 Theorem 1. -/
axiom trace_operator {p : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n} :
    SobolevSpace 1 p Ω →L[ℝ] TraceSpace p Ω

/-- The kernel of the trace operator is W_0^{1,p}.

Reference: Adams-Fournier Theorem 4.11; Evans §5.5 Theorem 2. -/
axiom trace_kernel_eq_sobolev_zero {p : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n}
    (u : SobolevSpace 1 p Ω) :
    trace_operator u = 0 ↔ ∃ u₀ : SobolevSpaceZero 1 p Ω,
      SobolevSpaceZero.toSobolevSpace u₀ = u

end NitscheHodgeLaplacian.Foundations
