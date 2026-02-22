import Mathlib.Tactic
import Mathlib.Analysis.InnerProductSpace.Basic
import NitscheHodgeLaplacian.Foundations.Domains
import NitscheHodgeLaplacian.Foundations.Sobolev
import NitscheHodgeLaplacian.Foundations.LpSpaces

/-!
# Differential Forms and the L²-de Rham Complex

This module axiomatizes FEEC (Finite Element Exterior Calculus) foundations:
Sobolev spaces of differential k-forms HΛᵏ(Ω), exterior derivative operators,
the de Rham chain complex property d ∘ d = 0, linking axioms to scalar Sobolev
spaces, trace operators on forms, and abbreviations for the classical n=3 spaces
(H(grad), H(curl), H(div)).

## Mathematical objects

- `SobolevFormSpace k Ω` — HΛᵏ(Ω) with graph norm ‖ω‖² = ‖ω‖²_{L²} + ‖dω‖²_{L²}
- `SobolevFormSpaceZero k Ω` — H₀Λᵏ(Ω) with zero boundary conditions
- `LpFormSpace p k Ω` — L^pΛᵏ(Ω)
- `extDerivativeZero`, `extDerivative` — exterior derivative operators
- Abbreviations: `HCurlSpace`, `HDivSpace`, `HGradSpace` for n=3

## References

- Arnold, Falk, Winther, "Finite element exterior calculus: from Hodge theory
  to numerical stability", Bull. AMS 2010.
- Arnold, Falk, Winther, "Finite element exterior calculus, homological techniques,
  and applications", Acta Numerica 2006.
- Tello Fachin, Tonnon, Zampa 2025, §2 for the specific setup.

## Design

The graph norm characterization `SobolevFormSpace.norm_sq_eq` and the
chain complex axiom `extDerivativeZero_comp_extDerivativeZero` encode the
fundamental structure of the de Rham complex. Forms above the ambient dimension
are trivially zero (`instUnique_above`), which is a topological fact.
-/

noncomputable section

open scoped NNReal ENNReal

namespace NitscheHodgeLaplacian.Foundations

/-! ## Opaque Types -/

/-- HΛᵏ(Ω): Sobolev space of k-forms with graph norm
‖ω‖² = ‖ω‖²_{L²} + ‖dω‖²_{L²}.

Reference: Arnold-Falk-Winther (2010) §3; Tello Fachin 2025, §2.1. -/
opaque SobolevFormSpace (k : ℕ) {n : ℕ} (Ω : BoundedDomain n) : Type

/-- H₀Λᵏ(Ω): Sobolev space of k-forms with zero boundary conditions.
These are the k-forms with vanishing tangential trace on ∂Ω.

Reference: Arnold-Falk-Winther (2010) §3.3. -/
opaque SobolevFormSpaceZero (k : ℕ) {n : ℕ} (Ω : BoundedDomain n) : Type

/-- L^p space of k-forms. Generalizes L²Λᵏ(Ω).

Reference: Arnold-Falk-Winther (2010) §3; used in Lp approximation bounds. -/
opaque LpFormSpace (p : ℝ≥0∞) (k : ℕ) {n : ℕ} (Ω : BoundedDomain n) : Type

/-- Trace space for k-forms on ∂Ω.

Reference: Arnold-Falk-Winther (2010) §3.4; trace of k-forms is the tangential component. -/
opaque FormTraceSpace (k : ℕ) {n : ℕ} (Ω : BoundedDomain n) : Type

/-! ## Instances: SobolevFormSpace -/

axiom SobolevFormSpace.instNormedAddCommGroup {k n : ℕ} {Ω : BoundedDomain n} :
    NormedAddCommGroup (SobolevFormSpace k Ω)

attribute [instance] SobolevFormSpace.instNormedAddCommGroup

axiom SobolevFormSpace.instNormedSpace {k n : ℕ} {Ω : BoundedDomain n} :
    NormedSpace ℝ (SobolevFormSpace k Ω)

attribute [instance] SobolevFormSpace.instNormedSpace

axiom SobolevFormSpace.instCompleteSpace {k n : ℕ} {Ω : BoundedDomain n} :
    CompleteSpace (SobolevFormSpace k Ω)

attribute [instance] SobolevFormSpace.instCompleteSpace

axiom SobolevFormSpace.instInnerProductSpace {k n : ℕ} {Ω : BoundedDomain n} :
    InnerProductSpace ℝ (SobolevFormSpace k Ω)

attribute [instance] SobolevFormSpace.instInnerProductSpace

/-! ## Instances: SobolevFormSpaceZero -/

axiom SobolevFormSpaceZero.instNormedAddCommGroup {k n : ℕ} {Ω : BoundedDomain n} :
    NormedAddCommGroup (SobolevFormSpaceZero k Ω)

attribute [instance] SobolevFormSpaceZero.instNormedAddCommGroup

axiom SobolevFormSpaceZero.instNormedSpace {k n : ℕ} {Ω : BoundedDomain n} :
    NormedSpace ℝ (SobolevFormSpaceZero k Ω)

attribute [instance] SobolevFormSpaceZero.instNormedSpace

axiom SobolevFormSpaceZero.instCompleteSpace {k n : ℕ} {Ω : BoundedDomain n} :
    CompleteSpace (SobolevFormSpaceZero k Ω)

attribute [instance] SobolevFormSpaceZero.instCompleteSpace

axiom SobolevFormSpaceZero.instInnerProductSpace {k n : ℕ} {Ω : BoundedDomain n} :
    InnerProductSpace ℝ (SobolevFormSpaceZero k Ω)

attribute [instance] SobolevFormSpaceZero.instInnerProductSpace

/-! ## Instances: LpFormSpace -/

axiom LpFormSpace.instNormedAddCommGroup {p : ℝ≥0∞} {k n : ℕ} {Ω : BoundedDomain n} :
    NormedAddCommGroup (LpFormSpace p k Ω)

attribute [instance] LpFormSpace.instNormedAddCommGroup

axiom LpFormSpace.instNormedSpace {p : ℝ≥0∞} {k n : ℕ} {Ω : BoundedDomain n} :
    NormedSpace ℝ (LpFormSpace p k Ω)

attribute [instance] LpFormSpace.instNormedSpace

axiom LpFormSpace.instCompleteSpace {p : ℝ≥0∞} {k n : ℕ} {Ω : BoundedDomain n} :
    CompleteSpace (LpFormSpace p k Ω)

attribute [instance] LpFormSpace.instCompleteSpace

axiom LpFormSpace.instInnerProductSpace {k n : ℕ} {Ω : BoundedDomain n} :
    InnerProductSpace ℝ (LpFormSpace 2 k Ω)

attribute [instance] LpFormSpace.instInnerProductSpace

/-! ## Instances: FormTraceSpace -/

axiom FormTraceSpace.instNormedAddCommGroup {k n : ℕ} {Ω : BoundedDomain n} :
    NormedAddCommGroup (FormTraceSpace k Ω)

attribute [instance] FormTraceSpace.instNormedAddCommGroup

axiom FormTraceSpace.instNormedSpace {k n : ℕ} {Ω : BoundedDomain n} :
    NormedSpace ℝ (FormTraceSpace k Ω)

attribute [instance] FormTraceSpace.instNormedSpace

/-! ## Inclusion: H₀Λᵏ ↪ HΛᵏ -/

/-- H₀Λᵏ(Ω) ↪ HΛᵏ(Ω): continuous inclusion of zero-boundary forms.

Reference: Arnold-Falk-Winther (2010) §3.3. -/
axiom SobolevFormSpaceZero.inclusion {k n : ℕ} {Ω : BoundedDomain n} :
    SobolevFormSpaceZero k Ω →L[ℝ] SobolevFormSpace k Ω

/-- The inclusion H₀Λᵏ ↪ HΛᵏ is injective. -/
axiom SobolevFormSpaceZero.inclusion_injective {k n : ℕ} {Ω : BoundedDomain n}
    {ω₁ ω₂ : SobolevFormSpaceZero k Ω}
    (h : SobolevFormSpaceZero.inclusion ω₁ = SobolevFormSpaceZero.inclusion ω₂) :
    ω₁ = ω₂

/-- The inclusion is norm-non-increasing. -/
axiom SobolevFormSpaceZero.inclusion_norm_le {k n : ℕ} {Ω : BoundedDomain n}
    (ω : SobolevFormSpaceZero k Ω) :
    ‖SobolevFormSpaceZero.inclusion ω‖ ≤ ‖ω‖

/-! ## Exterior Derivative -/

/-- Exterior derivative on H₀Λᵏ: d : H₀Λᵏ(Ω) →L[ℝ] H₀Λᵏ⁺¹(Ω).

Reference: Arnold-Falk-Winther (2010) §3.2; Tello Fachin 2025, §2. -/
axiom extDerivativeZero {k n : ℕ} {Ω : BoundedDomain n} :
    SobolevFormSpaceZero k Ω →L[ℝ] SobolevFormSpaceZero (k + 1) Ω

/-- Exterior derivative on HΛᵏ: d : HΛᵏ(Ω) →L[ℝ] HΛᵏ⁺¹(Ω).

Reference: Arnold-Falk-Winther (2010) §3.2. -/
axiom extDerivative {k n : ℕ} {Ω : BoundedDomain n} :
    SobolevFormSpace k Ω →L[ℝ] SobolevFormSpace (k + 1) Ω

/-! ## Chain Complex Property: d ∘ d = 0 -/

/-- d_{k+1} ∘ d_k = 0 on H₀Λᵏ (chain complex).

Reference: Arnold-Falk-Winther (2010) Lemma 2.1; classical d² = 0. -/
axiom extDerivativeZero_comp_extDerivativeZero {k n : ℕ} {Ω : BoundedDomain n} :
    (extDerivativeZero (k := k + 1) (n := n) (Ω := Ω)).comp
      (extDerivativeZero (k := k) (n := n) (Ω := Ω)) = 0

/-- d_{k+1} ∘ d_k = 0 on HΛᵏ (chain complex).

Reference: Arnold-Falk-Winther (2010) Lemma 2.1. -/
axiom extDerivative_comp_extDerivative {k n : ℕ} {Ω : BoundedDomain n} :
    (extDerivative (k := k + 1) (n := n) (Ω := Ω)).comp
      (extDerivative (k := k) (n := n) (Ω := Ω)) = 0

/-- The exterior derivative commutes with inclusion:
d ∘ incl = incl ∘ d₀.

Reference: Arnold-Falk-Winther (2010) §3.3 (naturality of the de Rham complex). -/
axiom extDerivative_inclusion_comm {k n : ℕ} {Ω : BoundedDomain n} :
    (extDerivative (k := k) (n := n) (Ω := Ω)).comp SobolevFormSpaceZero.inclusion =
    SobolevFormSpaceZero.inclusion.comp (extDerivativeZero (k := k) (n := n) (Ω := Ω))

/-! ## Vanishing Above Dimension -/

/-- Forms of degree > n are trivial (unique element).

Reference: Standard differential geometry; Λᵏ(ℝⁿ) = 0 for k > n. -/
axiom SobolevFormSpace.instUnique_above {k n : ℕ} {Ω : BoundedDomain n} (hk : n < k) :
    Unique (SobolevFormSpace k Ω)

/-- Zero-boundary forms of degree > n are trivial.

Reference: Same as above. -/
axiom SobolevFormSpaceZero.instUnique_above {k n : ℕ} {Ω : BoundedDomain n} (hk : n < k) :
    Unique (SobolevFormSpaceZero k Ω)

/-! ## Embeddings to LpFormSpace -/

/-- HΛᵏ(Ω) ↪ L^p Λᵏ(Ω): continuous embedding forgetting derivative regularity.

Reference: Arnold-Falk-Winther (2010) §3. -/
axiom SobolevFormSpace.toLpForm {p : ℝ≥0∞} {k n : ℕ} {Ω : BoundedDomain n} :
    SobolevFormSpace k Ω →L[ℝ] LpFormSpace p k Ω

/-- H₀Λᵏ(Ω) ↪ L^p Λᵏ(Ω): continuous embedding.

Reference: Arnold-Falk-Winther (2010) §3. -/
axiom SobolevFormSpaceZero.toLpForm {p : ℝ≥0∞} {k n : ℕ} {Ω : BoundedDomain n} :
    SobolevFormSpaceZero k Ω →L[ℝ] LpFormSpace p k Ω

/-! ## Graph Norm Characterization -/

/-- ‖ω‖²_{HΛᵏ} = ‖ω‖²_{L²Λᵏ} + ‖dω‖²_{L²Λᵏ⁺¹}.

Reference: Arnold-Falk-Winther (2010) §3.1, Definition 3.1. -/
axiom SobolevFormSpace.norm_sq_eq {k n : ℕ} {Ω : BoundedDomain n}
    (ω : SobolevFormSpace k Ω) :
    ‖ω‖ ^ 2 = ‖SobolevFormSpace.toLpForm (p := 2) ω‖ ^ 2 +
      ‖SobolevFormSpace.toLpForm (p := 2) (extDerivative ω)‖ ^ 2

/-! ## Linking Axioms: Forms ↔ Scalar Spaces -/

/-- HΛ⁰(Ω) ≃ W^{1,2}(Ω): 0-forms with gradient in L² = H¹.

Reference: Arnold-Falk-Winther (2010) §3.2 (identification via gradient). -/
axiom SobolevFormSpace.equivSobolevSpace_zero {n : ℕ} {Ω : BoundedDomain n} :
    SobolevFormSpace 0 Ω ≃L[ℝ] SobolevSpace 1 2 Ω

/-- H₀Λ⁰(Ω) ≃ W₀^{1,2}(Ω).

Reference: Arnold-Falk-Winther (2010) §3.3. -/
axiom SobolevFormSpaceZero.equivSobolevSpaceZero_zero {n : ℕ} {Ω : BoundedDomain n} :
    SobolevFormSpaceZero 0 Ω ≃L[ℝ] SobolevSpaceZero 1 2 Ω

/-- L^p Λ⁰(Ω) ≃ L^p(Ω): scalar functions as 0-forms.

Reference: Standard identification of 0-forms with scalar functions. -/
axiom LpFormSpace.equivLpSpace_zero {p : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n} :
    LpFormSpace p 0 Ω ≃L[ℝ] LpSpace p Ω

/-- HΛⁿ(Ω) ≃ L²(Ω): top-degree forms (trivial d, graph norm = L² norm).

Reference: Arnold-Falk-Winther (2010) §3.2. -/
axiom SobolevFormSpace.equivLpSpace_top {n : ℕ} {Ω : BoundedDomain n} :
    SobolevFormSpace n Ω ≃L[ℝ] LpSpace 2 Ω

/-! ## Trace on Forms -/

/-- Trace operator on k-forms: HΛᵏ(Ω) →L[ℝ] FormTraceSpace k Ω.

Reference: Arnold-Falk-Winther (2010) §3.4; trace = tangential component on ∂Ω. -/
axiom formTraceOperator {k n : ℕ} {Ω : BoundedDomain n} :
    SobolevFormSpace k Ω →L[ℝ] FormTraceSpace k Ω

/-- Trace kernel characterization for forms: Tr(ω) = 0 iff ω ∈ H₀Λᵏ.

Reference: Arnold-Falk-Winther (2010) §3.4, Theorem 3.1. -/
axiom formTrace_kernel_eq_zero {k n : ℕ} {Ω : BoundedDomain n}
    (ω : SobolevFormSpace k Ω) :
    formTraceOperator ω = 0 ↔
    ∃ ω₀ : SobolevFormSpaceZero k Ω, SobolevFormSpaceZero.inclusion ω₀ = ω

/-! ## Abbreviations for n = 3 -/

/-- H(curl, Ω) = HΛ¹(Ω) for n = 3.

Reference: Monk, "Finite Element Methods for Maxwell's Equations", Chapter 3. -/
abbrev HCurlSpace (Ω : BoundedDomain 3) := SobolevFormSpace 1 Ω

/-- H₀(curl, Ω) = H₀Λ¹(Ω) for n = 3. -/
abbrev HCurlSpaceZero (Ω : BoundedDomain 3) := SobolevFormSpaceZero 1 Ω

/-- H(div, Ω) = HΛ²(Ω) for n = 3.

Reference: Monk, "Finite Element Methods for Maxwell's Equations", Chapter 3. -/
abbrev HDivSpace (Ω : BoundedDomain 3) := SobolevFormSpace 2 Ω

/-- H₀(div, Ω) = H₀Λ²(Ω) for n = 3. -/
abbrev HDivSpaceZero (Ω : BoundedDomain 3) := SobolevFormSpaceZero 2 Ω

/-- H(grad, Ω) = HΛ⁰(Ω) for n = 3. -/
abbrev HGradSpace (Ω : BoundedDomain 3) := SobolevFormSpace 0 Ω

/-- H₀(grad, Ω) = H₀Λ⁰(Ω) for n = 3. -/
abbrev HGradSpaceZero (Ω : BoundedDomain 3) := SobolevFormSpaceZero 0 Ω

/-- Gradient operator: d₀ : H₀Λ⁰ → H₀Λ¹ for n = 3. -/
abbrev grad {Ω : BoundedDomain 3} := @extDerivativeZero 0 3 Ω

/-- Curl operator: d₁ : H₀Λ¹ → H₀Λ² for n = 3. -/
abbrev curl_ {Ω : BoundedDomain 3} := @extDerivativeZero 1 3 Ω

/-- Divergence operator: d₂ : H₀Λ² → H₀Λ³ for n = 3. -/
abbrev div_ {Ω : BoundedDomain 3} := @extDerivativeZero 2 3 Ω

/-- L²Λᵏ(Ω) as abbreviation for LpFormSpace 2 k Ω. -/
abbrev L2FormSpace (k : ℕ) {n : ℕ} (Ω : BoundedDomain n) := LpFormSpace 2 k Ω

end NitscheHodgeLaplacian.Foundations
