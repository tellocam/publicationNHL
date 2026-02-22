import Mathlib.Tactic
import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.LinearAlgebra.Dimension.Finrank
import NitscheHodgeLaplacian.Foundations.Domains
import NitscheHodgeLaplacian.Foundations.DifferentialForms

/-!
# FEEC Discretization

This module axiomatizes finite element exterior calculus (FEEC) discretization:
triangulations, discrete form spaces, discrete exterior derivatives, canonical
interpolation with commuting diagrams, discrete Hodge decomposition, and inverse
estimates.

## Mathematical objects

- `Triangulation Ω` — simplicial triangulation of Ω
- `DiscreteFormSpace k T` — V_h^k: discrete k-forms on T
- `DiscreteFormSpaceZero k T` — V_{h,0}^k: discrete k-forms with zero BC
- `canonicalInterpolation` — Falk-Winther smoothed interpolation Π_h
- `commuting_diagram` — d_h ∘ Π_h = Π_h ∘ d (fundamental FEEC property)
- `discreteHodgeProjectionExact/Coclosed` — orthogonal Hodge decomposition
- `inverse_estimate` — ‖ω_h‖_{HΛᵏ} ≤ C h⁻¹ ‖ω_h‖_{L²Λᵏ}

## References

- Arnold, Falk, Winther, "Finite element exterior calculus: from Hodge theory
  to numerical stability", Bull. AMS 2010. §4–§5 for discrete spaces and estimates.
- Falk, Winther, "Local bounded cochain projections", Math. Comp. 2014.
  The smoothed interpolation operator.
- Christiansen, Winther, "Smoothed projections in finite element exterior calculus",
  Math. Comp. 2008.
- Tello Fachin 2025, §2.2–§2.3 for the specific Nitsche framework.

## Design

The canonical interpolation `canonicalInterpolation` models the *smoothed*
(regularized) Falk-Winther operator, not the classical de Rham map. The classical
de Rham interpolation requires extra Sobolev regularity (W^{s,p} with sp > n)
and is not bounded on HΛᵏ. The smoothed variant is bounded on HΛᵏ while
preserving the commuting diagram and projection properties.
-/

noncomputable section

open scoped NNReal ENNReal

namespace NitscheHodgeLaplacian.Foundations

/-! ## Triangulation -/

/-- A triangulation of a bounded domain Ω.

Reference: Ern-Guermond §1.3; Ciarlet "The Finite Element Method" Chapter 2. -/
opaque Triangulation {n : ℕ} (Ω : BoundedDomain n) : Type

/-- The mesh size (maximum element diameter) of a triangulation. -/
axiom Triangulation.meshSize {n : ℕ} {Ω : BoundedDomain n} (T : Triangulation Ω) : ℝ

/-- Mesh size is positive. -/
axiom Triangulation.meshSize_pos {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) : T.meshSize > 0

/-- Shape regularity constant of a triangulation.

Reference: Ern-Guermond §1.3 (shape regularity = ratio of inradius to circumradius). -/
axiom Triangulation.shapeRegularity {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) : ℝ

/-- Shape regularity is positive. -/
axiom Triangulation.shapeRegularity_pos {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) : T.shapeRegularity > 0

/-- A refinement of the triangulation. -/
axiom Triangulation.refine {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) : Triangulation Ω

/-- Refinement strictly decreases the mesh size. -/
axiom Triangulation.refine_meshSize_lt {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) : T.refine.meshSize < T.meshSize

/-- Quasi-uniform mesh property for a triangulation.
A family of meshes is quasi-uniform if the mesh size h and the minimum element
diameter are bounded below by c*h for some fixed c > 0.

Reference: Brenner-Scott "The Mathematical Theory of Finite Element Methods" §4.4.
Used in inverse estimates and Lp approximation bounds. -/
opaque Triangulation.IsQuasiUniform {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) : Prop

/-- Polynomial degree r of the finite element space.

Reference: Tello Fachin 2025, §2.2; the paper requires r ≥ 2. -/
axiom Triangulation.polynomialDegree {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) : ℕ

/-- W^{l,p} Sobolev regularity space for k-forms (higher regularity than H^1).

Reference: Adams-Fournier Chapter 3; used in interpolation error estimates. -/
opaque WlpFormSpace (l : ℕ) (p : ℝ≥0∞) (k : ℕ) {n : ℕ} (Ω : BoundedDomain n) : Type

axiom WlpFormSpace.instNormedAddCommGroup {l : ℕ} {p : ℝ≥0∞} {k n : ℕ}
    {Ω : BoundedDomain n} : NormedAddCommGroup (WlpFormSpace l p k Ω)

attribute [instance] WlpFormSpace.instNormedAddCommGroup

axiom WlpFormSpace.instNormedSpace {l : ℕ} {p : ℝ≥0∞} {k n : ℕ}
    {Ω : BoundedDomain n} : NormedSpace ℝ (WlpFormSpace l p k Ω)

attribute [instance] WlpFormSpace.instNormedSpace

/-! ## Discrete Form Spaces -/

/-- Discrete k-form space V_h^k associated with a triangulation.

Reference: Arnold-Falk-Winther (2010) §4.2; these are the Raviart-Thomas,
Nédélec, and Whitney form spaces for k=0,1,2,3. -/
opaque DiscreteFormSpace (k : ℕ) {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) : Type

/-- Discrete k-form space V_{h,0}^k with zero boundary conditions.

Reference: Arnold-Falk-Winther (2010) §4.3. -/
opaque DiscreteFormSpaceZero (k : ℕ) {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) : Type

/-! ## Instances: DiscreteFormSpace -/

axiom DiscreteFormSpace.instNormedAddCommGroup {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : NormedAddCommGroup (DiscreteFormSpace k T)

attribute [instance] DiscreteFormSpace.instNormedAddCommGroup

axiom DiscreteFormSpace.instNormedSpace {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : NormedSpace ℝ (DiscreteFormSpace k T)

attribute [instance] DiscreteFormSpace.instNormedSpace

axiom DiscreteFormSpace.instInnerProductSpace {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : InnerProductSpace ℝ (DiscreteFormSpace k T)

attribute [instance] DiscreteFormSpace.instInnerProductSpace

axiom DiscreteFormSpace.instCompleteSpace {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : CompleteSpace (DiscreteFormSpace k T)

attribute [instance] DiscreteFormSpace.instCompleteSpace

axiom DiscreteFormSpace.instFiniteDimensional {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : FiniteDimensional ℝ (DiscreteFormSpace k T)

attribute [instance] DiscreteFormSpace.instFiniteDimensional

/-! ## Instances: DiscreteFormSpaceZero -/

axiom DiscreteFormSpaceZero.instNormedAddCommGroup {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : NormedAddCommGroup (DiscreteFormSpaceZero k T)

attribute [instance] DiscreteFormSpaceZero.instNormedAddCommGroup

axiom DiscreteFormSpaceZero.instNormedSpace {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : NormedSpace ℝ (DiscreteFormSpaceZero k T)

attribute [instance] DiscreteFormSpaceZero.instNormedSpace

axiom DiscreteFormSpaceZero.instInnerProductSpace {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : InnerProductSpace ℝ (DiscreteFormSpaceZero k T)

attribute [instance] DiscreteFormSpaceZero.instInnerProductSpace

axiom DiscreteFormSpaceZero.instCompleteSpace {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : CompleteSpace (DiscreteFormSpaceZero k T)

attribute [instance] DiscreteFormSpaceZero.instCompleteSpace

axiom DiscreteFormSpaceZero.instFiniteDimensional {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} : FiniteDimensional ℝ (DiscreteFormSpaceZero k T)

attribute [instance] DiscreteFormSpaceZero.instFiniteDimensional

/-! ## Embeddings into Continuous Spaces -/

/-- V_h^k ↪ HΛᵏ(Ω): conforming discrete space embeds into continuous space.

Reference: Arnold-Falk-Winther (2010) §4.3 (conforming = subcomplex). -/
axiom DiscreteFormSpace.toSobolevFormSpace {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    DiscreteFormSpace k T →L[ℝ] SobolevFormSpace k Ω

/-- V_{h,0}^k ↪ H₀Λᵏ(Ω): conforming discrete space with zero BC.

Reference: Arnold-Falk-Winther (2010) §4.3. -/
axiom DiscreteFormSpaceZero.toSobolevFormSpaceZero {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    DiscreteFormSpaceZero k T →L[ℝ] SobolevFormSpaceZero k Ω

/-- V_{h,0}^k ↪ V_h^k: inclusion of zero-boundary into full discrete space. -/
axiom DiscreteFormSpaceZero.toDiscreteFormSpace {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    DiscreteFormSpaceZero k T →L[ℝ] DiscreteFormSpace k T

/-- The embedding V_h^k ↪ HΛᵏ(Ω) is injective. -/
axiom DiscreteFormSpace.toSobolevFormSpace_injective {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} {ω₁ ω₂ : DiscreteFormSpace k T}
    (h : DiscreteFormSpace.toSobolevFormSpace ω₁ =
         DiscreteFormSpace.toSobolevFormSpace ω₂) :
    ω₁ = ω₂

/-- The embedding V_{h,0}^k ↪ H₀Λᵏ(Ω) is injective. -/
axiom DiscreteFormSpaceZero.toSobolevFormSpaceZero_injective {k n : ℕ}
    {Ω : BoundedDomain n} {T : Triangulation Ω}
    {ω₁ ω₂ : DiscreteFormSpaceZero k T}
    (h : DiscreteFormSpaceZero.toSobolevFormSpaceZero ω₁ =
         DiscreteFormSpaceZero.toSobolevFormSpaceZero ω₂) :
    ω₁ = ω₂

/-! ## Discrete Exterior Derivative -/

/-- Discrete exterior derivative: d_h : V_{h,0}^k →L[ℝ] V_{h,0}^{k+1}.

Reference: Arnold-Falk-Winther (2010) §4.3; the discrete complex
mirrors the continuous de Rham complex exactly. -/
axiom discreteExtDerivativeZero {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    DiscreteFormSpaceZero k T →L[ℝ] DiscreteFormSpaceZero (k + 1) T

/-- The discrete exterior derivative commutes with the inclusion into continuous spaces:
incl ∘ d_h = d ∘ incl.

Reference: Arnold-Falk-Winther (2010) §4.3 (subcomplex property). -/
axiom discreteExtDerivativeZero_comm {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    DiscreteFormSpaceZero.toSobolevFormSpaceZero.comp discreteExtDerivativeZero =
    extDerivativeZero.comp
      (DiscreteFormSpaceZero.toSobolevFormSpaceZero (k := k) (n := n) (Ω := Ω) (T := T))

/-- Discrete chain complex: d^h_{k+1} ∘ d^h_k = 0.

Reference: Arnold-Falk-Winther (2010) §4.3. -/
axiom discreteExtDerivativeZero_comp {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    (discreteExtDerivativeZero (k := k + 1) (n := n) (Ω := Ω) (T := T)).comp
      (discreteExtDerivativeZero (k := k) (n := n) (Ω := Ω) (T := T)) = 0

/-! ## Canonical Interpolation and Commuting Diagram -/

/-- Canonical interpolation operator: Π_h : H₀Λᵏ(Ω) →L[ℝ] V_{h,0}^k.

**Mathematical note**: This models the *smoothed* (regularized) canonical interpolation
of Falk-Winther (2014) / Christiansen-Winther (2008), NOT the classical de Rham map.
The classical de Rham interpolation requires extra Sobolev regularity (W^{s,p} with sp > n)
and is not bounded on HΛᵏ. The smoothed variant is bounded on HΛᵏ while preserving the
commuting diagram and projection properties.

Reference: Falk-Winther, Math. Comp. 2014; Christiansen-Winther, Math. Comp. 2008. -/
axiom canonicalInterpolation {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    SobolevFormSpaceZero k Ω →L[ℝ] DiscreteFormSpaceZero k T

/-- Projection property: Π_h ∘ incl = id on V_{h,0}^k.

Reference: Falk-Winther (2014) Theorem 1; Christiansen-Winther (2008) §3. -/
axiom canonicalInterpolation_projection {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} (ω_h : DiscreteFormSpaceZero k T) :
    canonicalInterpolation (DiscreteFormSpaceZero.toSobolevFormSpaceZero ω_h) = ω_h

/-- **THE COMMUTING DIAGRAM**: d_k^h ∘ Π_h^k = Π_h^{k+1} ∘ d_k.
This is the fundamental property of FEEC that ensures the discrete complex
inherits the cohomological structure of the continuous complex.

Reference: Arnold-Falk-Winther (2010) Theorem 5.1;
Falk-Winther (2014) Theorem 2. -/
axiom commuting_diagram {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    discreteExtDerivativeZero.comp (canonicalInterpolation (k := k) (n := n) (Ω := Ω) (T := T)) =
    (canonicalInterpolation (k := k + 1) (n := n) (Ω := Ω) (T := T)).comp extDerivativeZero

/-! ## Interpolation Error Estimate -/

/-- Interpolation error: ‖Π_h ω - ω‖ ≤ C h ‖ω‖.

First-order estimate for the smoothed canonical interpolation. Higher-order estimates
(O(h^l) for l ≥ 2) would require tracking higher Sobolev regularity of ω, which is
not available in H₀Λᵏ alone.

Reference: Arnold-Falk-Winther (2010) §5; Falk-Winther (2014) Theorem 3. -/
axiom interpolation_error {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    ∃ C : ℝ, C > 0 ∧ ∀ ω : SobolevFormSpaceZero k Ω,
      ‖DiscreteFormSpaceZero.toSobolevFormSpaceZero
        (canonicalInterpolation (T := T) ω) - ω‖ ≤
        C * T.meshSize * ‖ω‖

/-! ## Discrete Hodge Decomposition -/

/-- Orthogonal projection onto exact forms (image of d_{k-1}^h).

Reference: Arnold-Falk-Winther (2010) §3.1 (Hodge decomposition for discrete complex). -/
axiom discreteHodgeProjectionExact {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    DiscreteFormSpaceZero k T →L[ℝ] DiscreteFormSpaceZero k T

/-- Orthogonal projection onto coclosed forms (complement of exact forms).

Reference: Arnold-Falk-Winther (2010) §3.1. -/
axiom discreteHodgeProjectionCoclosed {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    DiscreteFormSpaceZero k T →L[ℝ] DiscreteFormSpaceZero k T

/-- P_d + P_z = id: the two projections sum to identity.

Reference: Arnold-Falk-Winther (2010) §3.1 Theorem 3.2. -/
axiom discreteHodge_sum_eq_id {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} (ω_h : DiscreteFormSpaceZero k T) :
    discreteHodgeProjectionExact ω_h + discreteHodgeProjectionCoclosed ω_h = ω_h

/-- P_d ∘ P_z = 0: exact projection annihilates coclosed component. -/
axiom discreteHodge_exact_coclosed_zero {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} (ω_h : DiscreteFormSpaceZero k T) :
    discreteHodgeProjectionExact (discreteHodgeProjectionCoclosed ω_h) = 0

/-- P_z ∘ P_d = 0: coclosed projection annihilates exact component. -/
axiom discreteHodge_coclosed_exact_zero {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} (ω_h : DiscreteFormSpaceZero k T) :
    discreteHodgeProjectionCoclosed (discreteHodgeProjectionExact ω_h) = 0

/-- P_d range = image of d_k^h: the exact projection on (k+1)-forms lands in the range of d_k.

Reference: Arnold-Falk-Winther (2010) §3.1. -/
axiom discreteHodgeProjectionExact_range {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} (ω_h : DiscreteFormSpaceZero (k + 1) T) :
    ∃ τ_h : DiscreteFormSpaceZero k T,
      discreteExtDerivativeZero τ_h = discreteHodgeProjectionExact ω_h

/-- P_z maps into ker(d_k^h): the coclosed component is in the kernel of d.

Reference: Arnold-Falk-Winther (2010) §3.1. -/
axiom discreteHodgeProjectionCoclosed_kernel {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} (ω_h : DiscreteFormSpaceZero k T) :
    discreteExtDerivativeZero (discreteHodgeProjectionCoclosed ω_h) = 0

/-! ## Inverse Estimate -/

/-- Inverse estimate: ‖ω_h‖_{HΛᵏ} ≤ C h⁻¹ ‖ω_h‖_{L²Λᵏ}.
Relates the stronger Sobolev norm to the weaker L² norm on finite element spaces.

Reference: Brenner-Scott §4.5 Theorem 4.5.11; Ern-Guermond §1.4.
The inverse estimate holds on quasi-uniform meshes; see also
Arnold-Falk-Winther (2010) Lemma 5.3. -/
axiom inverse_estimate {k n : ℕ} {Ω : BoundedDomain n}
    {T : Triangulation Ω} :
    ∃ C : ℝ, C > 0 ∧ ∀ ω_h : DiscreteFormSpaceZero k T,
      ‖DiscreteFormSpaceZero.toSobolevFormSpaceZero ω_h‖ ≤
        C * T.meshSize⁻¹ *
          ‖SobolevFormSpaceZero.toLpForm (p := 2)
            (DiscreteFormSpaceZero.toSobolevFormSpaceZero ω_h)‖

end NitscheHodgeLaplacian.Foundations
