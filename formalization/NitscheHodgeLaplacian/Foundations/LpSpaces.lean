import Mathlib.Tactic
import Mathlib.Analysis.InnerProductSpace.Basic
import NitscheHodgeLaplacian.Foundations.Domains
import NitscheHodgeLaplacian.Foundations.Sobolev

/-!
# L^p Spaces

This module axiomatizes L^p(Ω) spaces as opaque Banach spaces, with
identification to Sobolev spaces W^{0,p}(Ω), continuous embeddings,
and Hölder duality.

## Mathematical objects

- `LpSpace p Ω` — L^p(Ω) as a Banach space
- `SobolevSpace.equivLp` — W^{0,p}(Ω) ≃ L^p(Ω)
- `SobolevSpace.toLp` — W^{k,p}(Ω) ↪ L^p(Ω)
- `LpSpace.dualPairing` — L^p × L^q → ℝ for conjugate exponents
- `LpSpace.holder` — Hölder inequality

## References

- Adams & Fournier, "Sobolev Spaces", Chapter 2 (L^p theory).
- Brezis, "Functional Analysis, Sobolev Spaces and PDE", Springer 2011.
  §4 for L^p duality and Hölder.
- Folland, "Real Analysis" (2nd ed.), Wiley 1999. Chapter 6.

## Design

The exponent p : ℝ≥0∞ matches Mathlib's `MeasureTheory.Lp` convention,
allowing p = ∞ to represent L^∞. Conjugate exponents satisfy 1/p + 1/q = 1
in ℝ≥0∞ arithmetic.
-/

noncomputable section

open scoped NNReal ENNReal

namespace NitscheHodgeLaplacian.Foundations

/-! ## L^p Spaces -/

/-- L^p(Ω) as an opaque Banach space. Parameter p : ℝ≥0∞ matches Mathlib convention.

Reference: Adams-Fournier Chapter 2; Brezis §4. -/
opaque LpSpace (p : ℝ≥0∞) {n : ℕ} (Ω : BoundedDomain n) : Type

/-- L^p(Ω) is a normed add comm group. -/
axiom LpSpace.instNormedAddCommGroup {p : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n} :
    NormedAddCommGroup (LpSpace p Ω)

attribute [instance] LpSpace.instNormedAddCommGroup

/-- L^p(Ω) is a normed space over ℝ. -/
axiom LpSpace.instNormedSpace {p : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n} :
    NormedSpace ℝ (LpSpace p Ω)

attribute [instance] LpSpace.instNormedSpace

/-- L^p(Ω) is complete (Banach).

Reference: Riesz-Fischer theorem; Adams-Fournier Theorem 2.7. -/
axiom LpSpace.instCompleteSpace {p : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n} :
    CompleteSpace (LpSpace p Ω)

attribute [instance] LpSpace.instCompleteSpace

/-- L^2(Ω) is a Hilbert space (inner product space).

Reference: Adams-Fournier §2; the L^2 inner product ⟨f,g⟩ = ∫_Ω fg dx. -/
axiom LpSpace.instInnerProductSpace {n : ℕ} {Ω : BoundedDomain n} :
    InnerProductSpace ℝ (LpSpace 2 Ω)

attribute [instance] LpSpace.instInnerProductSpace

/-! ## Embeddings between Sobolev and L^p -/

/-- W^{0,p}(Ω) ≃ L^p(Ω): identification of zero-order Sobolev space with L^p.

Reference: Adams-Fournier §3; this is essentially definitional. -/
axiom SobolevSpace.equivLp {p : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n} :
    SobolevSpace 0 p Ω ≃L[ℝ] LpSpace p Ω

/-- W^{k,p}(Ω) ↪ L^p(Ω): continuous embedding forgetting derivatives.

Reference: Adams-Fournier §4.12 (trivial direction of Sobolev embedding). -/
axiom SobolevSpace.toLp {k : ℕ} {p : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n} :
    SobolevSpace k p Ω →L[ℝ] LpSpace p Ω

/-! ## Duality and Hölder -/

/-- Dual pairing: L^p × L^q → ℝ for conjugate exponents (1/p + 1/q = 1).

Reference: Adams-Fournier Theorem 2.3; Brezis Theorem 4.1. -/
axiom LpSpace.dualPairing {p q : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n}
    (hpq : p⁻¹ + q⁻¹ = 1) :
    LpSpace p Ω →L[ℝ] LpSpace q Ω →L[ℝ] ℝ

/-- Hölder inequality via the dual pairing.

Reference: Adams-Fournier Theorem 2.4; Brezis Theorem 4.6. -/
axiom LpSpace.holder {p q : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n}
    (hpq : p⁻¹ + q⁻¹ = 1) (f : LpSpace p Ω) (g : LpSpace q Ω) :
    |LpSpace.dualPairing hpq f g| ≤ ‖f‖ * ‖g‖

end NitscheHodgeLaplacian.Foundations
