import Mathlib.Analysis.SpecialFunctions.Pow.Real
import NitscheHodgeLaplacian.Foundations.Discretization

/-!
# Nitsche Hodge-Laplacian Error Estimates

Formalization of Theorem 4.1 and Corollary 4.2 from:

  L. Tello Fachin, G. Tonnon, F. Zampa,
  "Improved A Priori Error Estimates for the Nitsche Hodge-Laplacian",
  2025.

## Architecture

The proof chain decomposes the parametric error bound into granular axioms so that the
paper's novel contributions (Lemma 3.3, Corollary 3.7) are explicitly visible:

- `bnb_with_noncritical` [Established: Tello Fachin 2025, Lem. 3.5.2 + Thm. 3.3.7 + Thm. 3.5.3]
  — BNB quasi-optimality + non-critical terms
- `key_bound` [Novel: Tello Fachin 2025, Lemma 3.3] — adapts [AFG12, Thm 3.6] to 3D
- `lp_approximation` [Novel: Tello Fachin 2025, Corollary 3.7] — via Lp stability of P
  + discrete Poincaré

These chain into proved theorems:
- `critical_term_lp_bound` — chains key_bound + lp_approximation
- `parametric_error_bound` — chains bnb_with_noncritical + critical_term_lp_bound

Main results:
- `nitsche_hodge_laplacian_k1`: k=1 (H(curl)) rate h^{r-1/2} + h^{l-1/2}
- `nitsche_hodge_laplacian_k2`: k=2 (H(div)) rate h^{r-1/2} + h^{l-2/3}
- `nitsche_hodge_laplacian_k2_conditional`: Corollary 4.2 under discrete Hodge Lp stability

## Classification

Every declaration is tagged [Established: ...] or [Novel: ...] in its docstring.
Open problems are tagged [Open: ...].
-/

noncomputable section

open scoped NNReal ENNReal

namespace NitscheHodgeLaplacian

open NitscheHodgeLaplacian.Foundations

/-! ## Section 1: Paper-Specific Setup -/

/-- [Established: Tello Fachin 2025, Def. 3.1.7]
X-norm: mesh-dependent product norm on the error pair (ω_h, σ_h):
  ‖(ω_h, σ_h)‖²_X = ‖ω_h‖²_# + ‖σ_h‖²_{L²} + h² ‖d_{k-1}σ_h‖²_{L²}
where ‖ω_h‖²_# = ‖ω_h‖²_{L²} + ‖d_k ω_h‖²_{L²} + h⁻¹ ‖ω_h^tan‖²_{L²(∂Ω)}
                   + h ‖ω_h^nor‖²_{L²(∂Ω)}.
This is a real-valued function rather than a NormedAddCommGroup instance
because it depends on the mesh. -/
axiom XNorm {k n : ℕ} {Ω : BoundedDomain n} (T : Triangulation Ω)
    (e_ω : SobolevFormSpaceZero k Ω) (e_σ : SobolevFormSpaceZero (k - 1) Ω) : ℝ

/-- [Established: Tello Fachin 2025, Def. 3.1.7]
Non-negativity of the X-norm. -/
axiom XNorm_nonneg {k n : ℕ} {Ω : BoundedDomain n} {T : Triangulation Ω}
    {e_ω : SobolevFormSpaceZero k Ω} {e_σ : SobolevFormSpaceZero (k - 1) Ω} :
    0 ≤ XNorm T e_ω e_σ

/-! ### Intermediate Opaque Values -/

/-- [Novel: Tello Fachin 2025, proof of Thm 4.1]
Critical term supremum from BNB analysis with choice v_h = Pω:
  sup_{y_h ≠ 0} |⟨ω - Pω, d_{k-1}τ_h⟩| / ‖y_h‖_X
where P is the hybrid projection (Definition 2.1 of Tello Fachin 2025). -/
opaque criticalTermSup {k : ℕ} {Ω : BoundedDomain 3}
    (T : Triangulation Ω) (ω : SobolevFormSpaceZero k Ω) : ℝ

/-- [Novel: Tello Fachin 2025, proof of Thm 4.1]
Non-negativity of the critical term supremum. -/
axiom criticalTermSup_nonneg {k : ℕ} {Ω : BoundedDomain 3}
    {T : Triangulation Ω} {ω : SobolevFormSpaceZero k Ω} :
    0 ≤ criticalTermSup T ω

/-- [Novel: Tello Fachin 2025, Corollary 3.7]
Lp approximation error ‖ω - Pω‖_{Lp} for the hybrid projection P. -/
opaque lpApproxError {k : ℕ} {Ω : BoundedDomain 3}
    (T : Triangulation Ω) (p : ℝ) (ω : SobolevFormSpaceZero k Ω) : ℝ

/-- [Novel: Tello Fachin 2025, Corollary 3.7]
Non-negativity of the Lp approximation error. -/
axiom lpApproxError_nonneg {k : ℕ} {Ω : BoundedDomain 3}
    {T : Triangulation Ω} {p : ℝ} {ω : SobolevFormSpaceZero k Ω} :
    0 ≤ lpApproxError T p ω

/-! ## Section 2: Proof Chain Axioms -/

/-- [Open: Issue #2 in Tello Fachin 2025]
Opaque hypothesis for Corollary 4.2: discrete Hodge decomposition is Lp-stable.
Open problem for k=2, p>6. -/
opaque DiscreteHodgeLpStable {n : ℕ} {Ω : BoundedDomain n}
    (T : Triangulation Ω) (p : ℝ) : Prop

/-- [Established: Tello Fachin 2025, Lem. 3.5.2 + Thm. 3.3.7 + Thm. 3.5.3]
BNB quasi-optimality combined with non-critical error terms.
The critical term (involving d_{k-1} test functions) is separated out as `criticalTermSup`. -/
axiom bnb_with_noncritical (k : ℕ) {Ω : BoundedDomain 3}
    (T : Triangulation Ω) (hr : T.polynomialDegree ≥ 2)
    (ω : SobolevFormSpaceZero k Ω)
    (e_ω : SobolevFormSpaceZero k Ω) (e_σ : SobolevFormSpaceZero (k - 1) Ω) :
    ∃ C : ℝ, C > 0 ∧
      XNorm T e_ω e_σ ≤
        C * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) + criticalTermSup T ω)

/-- [Novel: Tello Fachin 2025, Lemma 3.3]
Key bound on the critical term: adapts [AFG12, Thm. 3.6] to 3D.
Bounds the critical supremum by h^{-1/2-1/p} times the Lp approximation error.

Reference: Arnold, Falk, Winther (2006/2010) FEEC papers for the original Thm 3.6. -/
axiom key_bound (k : ℕ) {Ω : BoundedDomain 3}
    (T : Triangulation Ω) (hqu : T.IsQuasiUniform)
    (ω : SobolevFormSpaceZero k Ω) (p : ℝ) (hp : 2 ≤ p) :
    ∃ C : ℝ, C > 0 ∧
      criticalTermSup T ω ≤
        C * T.meshSize ^ (-(1 : ℝ)/2 - 1/p) * lpApproxError T p ω

/-- [Novel: Tello Fachin 2025, Corollary 3.7]
Lp approximation bound via Theorem 3.6 (Lp stability of P) and Lemma 3.5
(discrete Poincaré). Bounds ‖ω - Pω‖_{Lp} by h^l ‖ω‖_{W^{l,p}}.

The discrete Poincaré inequality (Lemma 3.5) follows from
Christiansen-Winther (2008) and the p-independent estimate of Costabel-McIntosh. -/
axiom lp_approximation (k : ℕ) {Ω : BoundedDomain 3}
    (T : Triangulation Ω) (hconv : Ω.IsConvex) (hcontr : Ω.IsContractible)
    (hqu : T.IsQuasiUniform)
    (ω : SobolevFormSpaceZero k Ω)
    (l : ℕ) (hl : 1 ≤ l) (p : ℝ) (hp : 2 ≤ p)
    (ω_Wlp_norm : ℝ) (h_nn : 0 ≤ ω_Wlp_norm) :
    ∃ C : ℝ, C > 0 ∧
      lpApproxError T p ω ≤ C * T.meshSize ^ (l : ℝ) * ω_Wlp_norm

/-- [Novel: Tello Fachin 2025, Thm 4.1(a), limit step]
k=1 final bound via limit p → ∞. Uses [SW82, Thm. 5.1] for C(p) = O(1) when k=1.
The limit step is axiomatized because h^{l-1/2-1/p} → h^{l-1/2} (increasing for h<1)
while C(p)‖ω‖_{W^{l,p}} → ‖ω‖_{W^{l,∞}} (decreasing), requiring genuine analysis.

Reference: Schatz-Wahlbin (1982) for the L∞ stability result. -/
axiom k1_limit_bound {Ω : BoundedDomain 3}
    (T : Triangulation Ω)
    (hconv : Ω.IsConvex) (hcontr : Ω.IsContractible)
    (hqu : T.IsQuasiUniform) (hr : T.polynomialDegree ≥ 2)
    (ω : SobolevFormSpaceZero 1 Ω)
    (l : ℕ) (hl : 1 ≤ l)
    (e_ω : SobolevFormSpaceZero 1 Ω) (e_σ : SobolevFormSpaceZero 0 Ω)
    (ω_Wlinfty_norm : ℝ) (h_nn : 0 ≤ ω_Wlinfty_norm) :
    ∃ C : ℝ, C > 0 ∧
      XNorm T e_ω e_σ ≤
        C * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) +
             T.meshSize ^ ((l : ℝ) - 1/2) * ω_Wlinfty_norm)

/-- [Novel: Tello Fachin 2025, Corollary 4.2]
Under the assumption that discrete Hodge decomposition is Lp-stable for all p ≥ 2,
the H(div) error achieves the same rate as H(curl): h^{r-1/2} + h^{l-1/2}. -/
axiom k2_conditional_lp_stability {Ω : BoundedDomain 3}
    (T : Triangulation Ω)
    (hconv : Ω.IsConvex) (hcontr : Ω.IsContractible)
    (hqu : T.IsQuasiUniform) (hr : T.polynomialDegree ≥ 2)
    (ω : SobolevFormSpaceZero 2 Ω)
    (l : ℕ) (hl : 1 ≤ l)
    (e_ω : SobolevFormSpaceZero 2 Ω) (e_σ : SobolevFormSpaceZero 1 Ω)
    (ω_Wlinfty_norm : ℝ) (h_nn : 0 ≤ ω_Wlinfty_norm)
    (h_lp_stable : ∀ p : ℝ, 2 ≤ p → DiscreteHodgeLpStable T p) :
    ∃ C : ℝ, C > 0 ∧
      XNorm T e_ω e_σ ≤
        C * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) +
             T.meshSize ^ ((l : ℝ) - 1/2) * ω_Wlinfty_norm)

/-! ## Section 3: Proved Theorems -/

/-- Exponent arithmetic for k=2 case: l - 1/2 - 1/6 = l - 2/3. -/
lemma exponent_k2 (l : ℕ) : (l : ℝ) - 1/2 - 1/6 = (l : ℝ) - 2/3 := by ring

/-- Helper: absorb two constants into one via max.
If a ≤ C₁ * (x + C₂ * y), then a ≤ (C₁ * max 1 C₂) * (x + y). -/
private lemma absorb_constant {C₁ C₂ a x y : ℝ}
    (hC₁ : 0 < C₁) (_hC₂ : 0 < C₂) (_ha : 0 ≤ a) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (hbound : a ≤ C₁ * (x + C₂ * y)) :
    a ≤ C₁ * max 1 C₂ * (x + y) := by
  have hm : 1 ≤ max 1 C₂ := le_max_left 1 C₂
  have hC₂m : C₂ ≤ max 1 C₂ := le_max_right 1 C₂
  calc a ≤ C₁ * (x + C₂ * y) := hbound
    _ ≤ C₁ * (max 1 C₂ * x + max 1 C₂ * y) := by
        apply mul_le_mul_of_nonneg_left _ (le_of_lt hC₁)
        apply add_le_add
        · exact le_mul_of_one_le_left hx hm
        · exact mul_le_mul_of_nonneg_right hC₂m hy
    _ = C₁ * max 1 C₂ * (x + y) := by ring

/-- [Novel: Tello Fachin 2025, chains Lemma 3.3 + Corollary 3.7]
Bound on the critical term in terms of mesh size and W^{l,p} norm:
  criticalTermSup T ω ≤ C * h^{l - 1/2 - 1/p} * ‖ω‖_{W^{l,p}}

Proof: chain `key_bound` (h^{-1/2-1/p} * lpApproxError) with
`lp_approximation` (lpApproxError ≤ h^l * ‖ω‖), combine exponents. -/
theorem critical_term_lp_bound (k : ℕ) {Ω : BoundedDomain 3}
    (T : Triangulation Ω) (hconv : Ω.IsConvex) (hcontr : Ω.IsContractible)
    (hqu : T.IsQuasiUniform)
    (ω : SobolevFormSpaceZero k Ω)
    (l : ℕ) (hl : 1 ≤ l) (p : ℝ) (hp : 2 ≤ p)
    (ω_Wlp_norm : ℝ) (h_nn : 0 ≤ ω_Wlp_norm) :
    ∃ C : ℝ, C > 0 ∧
      criticalTermSup T ω ≤
        C * T.meshSize ^ ((l : ℝ) - 1/2 - 1/p) * ω_Wlp_norm := by
  obtain ⟨C₁, hC₁, h_key⟩ := key_bound k T hqu ω p hp
  obtain ⟨C₂, hC₂, h_lp⟩ := lp_approximation k T hconv hcontr hqu ω l hl p hp ω_Wlp_norm h_nn
  have hh_pos := Triangulation.meshSize_pos T
  have hh := le_of_lt hh_pos
  -- Combine rpow exponents: h^{-1/2-1/p} * h^l = h^{l-1/2-1/p}
  have hrpow : T.meshSize ^ (-(1 : ℝ)/2 - 1/p) * T.meshSize ^ (l : ℝ) =
      T.meshSize ^ ((l : ℝ) - 1/2 - 1/p) := by
    rw [← Real.rpow_add hh_pos]; congr 1; ring
  refine ⟨C₁ * C₂, mul_pos hC₁ hC₂, ?_⟩
  -- Chain key_bound with lp_approximation
  have h1 : criticalTermSup T ω ≤
      C₁ * T.meshSize ^ (-(1 : ℝ)/2 - 1/p) * (C₂ * T.meshSize ^ (l : ℝ) * ω_Wlp_norm) :=
    le_trans h_key (mul_le_mul_of_nonneg_left h_lp
      (mul_nonneg (le_of_lt hC₁) (Real.rpow_nonneg hh _)))
  -- Rearrange: a*b*(c*d*e) → a*c*(b*d)*e, then combine b*d via hrpow
  have rearrange : ∀ (a b c d e : ℝ), a * b * (c * d * e) = a * c * (b * d) * e :=
    fun _ _ _ _ _ => by ring
  rw [rearrange, hrpow] at h1
  exact h1

/-- [Proved: chains bnb_with_noncritical + critical_term_lp_bound]
Parametric error bound: for form degree k, with mesh polynomial degree r ≥ 2,
regularity l ≥ 1, and integrability p ≥ 2, the X-norm error satisfies:

  ‖error‖_X ≤ C (h^{r-1/2} + h^{l-1/2-1/p} ‖ω‖_{W^{l,p}})

Previously a monolithic axiom (A1), now proved by chaining:
1. `bnb_with_noncritical` — BNB framework separating critical term
   [Tello Fachin 2025, Lem. 3.5.2 + Thm. 3.3.7 + Thm. 3.5.3]
2. `critical_term_lp_bound` — bounding critical term via key_bound + lp_approximation
   [Tello Fachin 2025, Lemma 3.3 + Corollary 3.7] -/
theorem parametric_error_bound (k : ℕ) {Ω : BoundedDomain 3}
    (T : Triangulation Ω)
    (hconv : Ω.IsConvex) (hcontr : Ω.IsContractible)
    (hqu : T.IsQuasiUniform) (hr : T.polynomialDegree ≥ 2)
    (ω : SobolevFormSpaceZero k Ω)
    (l : ℕ) (hl : 1 ≤ l) (p : ℝ) (hp : 2 ≤ p)
    (e_ω : SobolevFormSpaceZero k Ω) (e_σ : SobolevFormSpaceZero (k - 1) Ω)
    (ω_Wlp_norm : ℝ) (h_nn : 0 ≤ ω_Wlp_norm) :
    ∃ C : ℝ, C > 0 ∧
      XNorm T e_ω e_σ ≤
        C * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) +
             T.meshSize ^ ((l : ℝ) - 1/2 - 1/p) * ω_Wlp_norm) := by
  obtain ⟨C₁, hC₁, h_bnb⟩ := bnb_with_noncritical k T hr ω e_ω e_σ
  obtain ⟨C₂, hC₂, h_crit⟩ := critical_term_lp_bound k T hconv hcontr hqu ω l hl p hp
    ω_Wlp_norm h_nn
  refine ⟨C₁ * max 1 C₂, mul_pos hC₁ (lt_of_lt_of_le one_pos (le_max_left 1 C₂)), ?_⟩
  apply absorb_constant hC₁ hC₂ XNorm_nonneg
    (Real.rpow_nonneg (le_of_lt (Triangulation.meshSize_pos T)) _)
    (mul_nonneg (Real.rpow_nonneg (le_of_lt (Triangulation.meshSize_pos T)) _) h_nn)
  calc XNorm T e_ω e_σ
      ≤ C₁ * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) + criticalTermSup T ω) := h_bnb
    _ ≤ C₁ * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) +
          C₂ * T.meshSize ^ ((l : ℝ) - 1/2 - 1/p) * ω_Wlp_norm) := by
        apply mul_le_mul_of_nonneg_left _ (le_of_lt hC₁)
        linarith
    _ = C₁ * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) +
          C₂ * (T.meshSize ^ ((l : ℝ) - 1/2 - 1/p) * ω_Wlp_norm)) := by ring

/-- [Proved: applies k1_limit_bound]
**Theorem 4.1a (k=1, H(curl))**: For the Nitsche discretization of the
Hodge-Laplacian on 1-forms (Tello Fachin, Tonnon, Zampa 2025, Theorem 4.1a),
the X-norm error satisfies

  ‖error‖_X ≤ C (h^{r-1/2} + h^{l-1/2} ‖ω‖_{W^{l,∞}})

where r is the polynomial degree and l is the Sobolev regularity index.
The improved rate h^{l-1/2} (vs. the generic h^{l-2/3}) is achieved via a p → ∞
limit argument exploiting L∞ stability for k=1. -/
theorem nitsche_hodge_laplacian_k1 {Ω : BoundedDomain 3}
    (T : Triangulation Ω)
    (hconv : Ω.IsConvex) (hcontr : Ω.IsContractible)
    (hqu : T.IsQuasiUniform) (hr : T.polynomialDegree ≥ 2)
    (ω : SobolevFormSpaceZero 1 Ω)
    (l : ℕ) (hl : 1 ≤ l)
    (e_ω : SobolevFormSpaceZero 1 Ω) (e_σ : SobolevFormSpaceZero 0 Ω)
    (ω_Wlinfty_norm : ℝ) (h_nn : 0 ≤ ω_Wlinfty_norm) :
    ∃ C : ℝ, C > 0 ∧
      XNorm T e_ω e_σ ≤
        C * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) +
             T.meshSize ^ ((l : ℝ) - 1/2) * ω_Wlinfty_norm) :=
  k1_limit_bound T hconv hcontr hqu hr ω l hl e_ω e_σ ω_Wlinfty_norm h_nn

/-- [Proved: applies parametric_error_bound at p=6]
**Theorem 4.1b (k=2, H(div))**: For the Nitsche discretization of the
Hodge-Laplacian on 2-forms (Tello Fachin, Tonnon, Zampa 2025, Theorem 4.1b),
the X-norm error satisfies

  ‖error‖_X ≤ C (h^{r-1/2} + h^{l-2/3} ‖ω‖_{W^{l,6}})

where r is the polynomial degree, l is the regularity index, and the exponent
l-2/3 arises from substituting p=6 in the parametric bound (l-1/2-1/6 = l-2/3).
The choice p=6 is natural in 3D from the Sobolev embedding W^{1,2} ↪ L^6. -/
theorem nitsche_hodge_laplacian_k2 {Ω : BoundedDomain 3}
    (T : Triangulation Ω)
    (hconv : Ω.IsConvex) (hcontr : Ω.IsContractible)
    (hqu : T.IsQuasiUniform) (hr : T.polynomialDegree ≥ 2)
    (ω : SobolevFormSpaceZero 2 Ω)
    (l : ℕ) (hl : 1 ≤ l)
    (e_ω : SobolevFormSpaceZero 2 Ω) (e_σ : SobolevFormSpaceZero 1 Ω)
    (ω_Wl6_norm : ℝ) (h_nn : 0 ≤ ω_Wl6_norm) :
    ∃ C : ℝ, C > 0 ∧
      XNorm T e_ω e_σ ≤
        C * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) +
             T.meshSize ^ ((l : ℝ) - 2/3) * ω_Wl6_norm) := by
  obtain ⟨C, hC, hbound⟩ :=
    parametric_error_bound 2 T hconv hcontr hqu hr ω l hl 6 (by norm_num) e_ω e_σ
      ω_Wl6_norm h_nn
  exact ⟨C, hC, by rwa [show (l : ℝ) - 1/2 - 1/6 = (l : ℝ) - 2/3 from by ring] at hbound⟩

/-- [Proved: applies k2_conditional_lp_stability]
**Corollary 4.2 (k=2, conditional)**: Under the additional hypothesis that discrete
Hodge decomposition is Lp-stable for all p ≥ 2 (Tello Fachin, Tonnon, Zampa 2025,
Corollary 4.2), the H(div) error achieves the same improved rate as H(curl):

  ‖error‖_X ≤ C (h^{r-1/2} + h^{l-1/2} ‖ω‖_{W^{l,∞}})

The discrete Hodge Lp stability is an open problem for k=2 and p > 6. -/
theorem nitsche_hodge_laplacian_k2_conditional {Ω : BoundedDomain 3}
    (T : Triangulation Ω)
    (hconv : Ω.IsConvex) (hcontr : Ω.IsContractible)
    (hqu : T.IsQuasiUniform) (hr : T.polynomialDegree ≥ 2)
    (ω : SobolevFormSpaceZero 2 Ω)
    (l : ℕ) (hl : 1 ≤ l)
    (e_ω : SobolevFormSpaceZero 2 Ω) (e_σ : SobolevFormSpaceZero 1 Ω)
    (ω_Wlinfty_norm : ℝ) (h_nn : 0 ≤ ω_Wlinfty_norm)
    (h_lp_stable : ∀ p : ℝ, 2 ≤ p → DiscreteHodgeLpStable T p) :
    ∃ C : ℝ, C > 0 ∧
      XNorm T e_ω e_σ ≤
        C * (T.meshSize ^ ((T.polynomialDegree : ℝ) - 1/2) +
             T.meshSize ^ ((l : ℝ) - 1/2) * ω_Wlinfty_norm) :=
  k2_conditional_lp_stability T hconv hcontr hqu hr ω l hl e_ω e_σ ω_Wlinfty_norm h_nn
    h_lp_stable

end NitscheHodgeLaplacian
