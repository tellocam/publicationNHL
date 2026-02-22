import Mathlib.Tactic
import Mathlib.Analysis.InnerProductSpace.Basic
import NitscheHodgeLaplacian.Foundations.Domains

/-!
# Sobolev Spaces

This module axiomatizes Sobolev spaces W^{k,p}(Ω) and W_0^{k,p}(Ω),
their algebraic and topological instances, and weak derivatives.

## Mathematical objects

- `SobolevSpace k p Ω` — W^{k,p}(Ω): functions with weak derivatives
  up to order k in L^p
- `SobolevSpaceZero k p Ω` — W_0^{k,p}(Ω): closure of C_c^∞(Ω) in W^{k,p}
- `HasWeakDerivative` — predicate for weak partial derivatives

## References

- Adams & Fournier, "Sobolev Spaces" (2nd ed.), Academic Press, 2003.
  The standard reference for W^{k,p} spaces, §3–§5.
- Evans, "Partial Differential Equations" (2nd ed.), AMS, 2010. §5.

## Design

All spaces are opaque types equipped with normed space instances via axioms.
The exponent p : ℝ≥0∞ follows Mathlib's MeasureTheory.Lp convention.
Homogeneous boundary conditions are captured by SobolevSpaceZero
(closure of compactly supported smooth functions), matching the definition
used in [Tello Fachin 2025, §2].
-/

noncomputable section

open scoped NNReal ENNReal

namespace NitscheHodgeLaplacian.Foundations

/-! ## Sobolev Spaces -/

/-- The Sobolev space W^{k,p}(Ω) for k : ℕ, p : ℝ≥0∞, over a bounded domain.
Elements are Lp functions whose weak derivatives up to order k are also in Lp.

Reference: Adams-Fournier, Chapter 3. -/
opaque SobolevSpace (k : ℕ) (p : ℝ≥0∞) {n : ℕ} (Ω : BoundedDomain n) : Type

/-- The Sobolev space W_0^{k,p}(Ω): closure of C_c^∞(Ω) in W^{k,p}(Ω).
Functions in this space satisfy homogeneous boundary conditions.

Reference: Adams-Fournier, §5.5. -/
opaque SobolevSpaceZero (k : ℕ) (p : ℝ≥0∞) {n : ℕ} (Ω : BoundedDomain n) : Type

/-- W_0^{k,p} embeds into W^{k,p}. -/
axiom SobolevSpaceZero.toSobolevSpace {k : ℕ} {p : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n} :
    SobolevSpaceZero k p Ω → SobolevSpace k p Ω

/-- W^{k,p}(Ω) is a normed add comm group. -/
axiom SobolevSpace.instNormedAddCommGroup {k : ℕ} {p : ℝ≥0∞} {n : ℕ}
    {Ω : BoundedDomain n} : NormedAddCommGroup (SobolevSpace k p Ω)

attribute [instance] SobolevSpace.instNormedAddCommGroup

/-- W^{k,p}(Ω) is a normed space over ℝ. -/
axiom SobolevSpace.instNormedSpace {k : ℕ} {p : ℝ≥0∞} {n : ℕ}
    {Ω : BoundedDomain n} : NormedSpace ℝ (SobolevSpace k p Ω)

attribute [instance] SobolevSpace.instNormedSpace

/-- W_0^{k,p}(Ω) is a normed add comm group. -/
axiom SobolevSpaceZero.instNormedAddCommGroup {k : ℕ} {p : ℝ≥0∞} {n : ℕ}
    {Ω : BoundedDomain n} : NormedAddCommGroup (SobolevSpaceZero k p Ω)

attribute [instance] SobolevSpaceZero.instNormedAddCommGroup

/-- W_0^{k,p}(Ω) is a normed space over ℝ. -/
axiom SobolevSpaceZero.instNormedSpace {k : ℕ} {p : ℝ≥0∞} {n : ℕ}
    {Ω : BoundedDomain n} : NormedSpace ℝ (SobolevSpaceZero k p Ω)

attribute [instance] SobolevSpaceZero.instNormedSpace

/-- W_0^{1,2}(Ω) is a Hilbert space (inner product space).

Reference: Adams-Fournier §7; the H^1_0 Hilbert structure is fundamental
to the Lax-Milgram theorem and weak formulations. -/
axiom SobolevSpaceZero.instInnerProductSpace {n : ℕ} {Ω : BoundedDomain n} :
    InnerProductSpace ℝ (SobolevSpaceZero 1 2 Ω)

attribute [instance] SobolevSpaceZero.instInnerProductSpace

/-- W^{k,p}(Ω) is complete. -/
axiom SobolevSpace.instCompleteSpace {k : ℕ} {p : ℝ≥0∞} {n : ℕ}
    {Ω : BoundedDomain n} : CompleteSpace (SobolevSpace k p Ω)

attribute [instance] SobolevSpace.instCompleteSpace

/-- W_0^{k,p}(Ω) is complete. -/
axiom SobolevSpaceZero.instCompleteSpace {k : ℕ} {p : ℝ≥0∞} {n : ℕ}
    {Ω : BoundedDomain n} : CompleteSpace (SobolevSpaceZero k p Ω)

attribute [instance] SobolevSpaceZero.instCompleteSpace

/-! ## Weak Derivatives -/

/-- `HasWeakDerivative Ω u α v` states that `v` is the weak derivative
D^α u on the domain Ω: for all φ ∈ C_c^∞(Ω),
∫_Ω u D^α φ dx = (-1)^|α| ∫_Ω v φ dx.

Reference: Adams-Fournier §1.57; Evans §5.2. -/
axiom HasWeakDerivative {n : ℕ} (Ω : BoundedDomain n)
    (u : SobolevSpace 0 1 Ω) (α : Fin n → ℕ) (v : SobolevSpace 0 1 Ω) : Prop

/-- Weak derivatives are unique (when they exist).

Reference: Adams-Fournier Theorem 1.56. -/
axiom HasWeakDerivative.unique {n : ℕ} {Ω : BoundedDomain n}
    {u : SobolevSpace 0 1 Ω} {α : Fin n → ℕ}
    {v₁ v₂ : SobolevSpace 0 1 Ω}
    (h₁ : HasWeakDerivative Ω u α v₁) (h₂ : HasWeakDerivative Ω u α v₂) :
    v₁ = v₂

end NitscheHodgeLaplacian.Foundations
