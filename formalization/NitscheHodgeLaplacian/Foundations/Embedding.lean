import NitscheHodgeLaplacian.Foundations.Sobolev
import NitscheHodgeLaplacian.Foundations.LpSpaces
import NitscheHodgeLaplacian.Foundations.DifferentialForms
import NitscheHodgeLaplacian.Foundations.Domains

/-!
# Sobolev Embedding and Poincaré Inequality

This module axiomatizes embedding theorems for Sobolev spaces, including
Sobolev embeddings into L^q, compact Rellich-Kondrachov embeddings, and
the p-independent Poincaré-Friedrichs inequality on convex domains.

## Mathematical objects

- `sobolev_embedding_Lp` — W^{k,p}(Ω) ↪ L^q(Ω) (Sobolev embedding)
- `poincare_friedrichs_p_independent` — p-independent Poincaré inequality
  for differential forms on convex domains

## References

- Adams & Fournier, "Sobolev Spaces" (2nd ed.), §4–§5 (Sobolev embedding theorem).
- Costabel, McIntosh, "On Bogovskiĭ and regularization operators...",
  JEMS 2010 — p-independent estimate for vector fields.
- Diening, Ružička, "Calderón-Zygmund operators...", Ann. Math. 2010
  — p-independent Poincaré for forms.
- Tello Fachin, Tonnon, Zampa, 2025, §3.5 — application to discrete Poincaré.

## Design

The p-independent Poincaré-Friedrichs inequality (`poincare_friedrichs_p_independent`)
is the key analytical input that breaks the p-dependence of the constant C in
the approximation bounds. It requires convexity of Ω and is opaque here because
its proof in the literature (Costabel-McIntosh, Diening-Ružička) is technically involved.
-/

noncomputable section

open scoped NNReal ENNReal

namespace NitscheHodgeLaplacian.Foundations

/-! ## Sobolev Embedding into L^q -/

/-- Sobolev embedding into L^q: W^{k,p}(Ω) →L[ℝ] L^q(Ω).
Valid when kp < n and q ≤ np/(n - kp) (subcritical case).

Reference: Adams-Fournier Theorem 4.12 (Part I). -/
axiom sobolev_embedding_Lp {k : ℕ} {p q : ℝ≥0∞} {n : ℕ} {Ω : BoundedDomain n}
    (hkpn : k * p < n) (hq : q ≤ n * p / (n - k * p)) :
    SobolevSpace k p Ω →L[ℝ] LpSpace q Ω

/-! ## p-Independent Poincaré-Friedrichs -/

/-- **p-independent Poincaré-Friedrichs inequality** (key result from Costabel-McIntosh):
On convex domains, for ω ∈ H₀Λᵏ: ‖ω‖_{L^p} ≤ C ‖dω‖_{L^p} with C independent of p.

This is the crucial ingredient allowing the limit p → ∞ in the proof of
Theorem 4.1a of [Tello Fachin, Tonnon, Zampa 2025].

Reference:
- Costabel, McIntosh, JEMS 2010.
- Diening, Ružička, Ann. Math. 2010.
- Applied in: Tello Fachin 2025, §3.5, Lemma 3.5.2. -/
axiom poincare_friedrichs_p_independent {k : ℕ} {n : ℕ} {Ω : BoundedDomain n}
    (hconv : Ω.IsConvex) :
    ∃ C : ℝ, C > 0 ∧ ∀ (p : ℝ≥0∞) (_hp : 1 ≤ p) (_hp' : p < ⊤),
      ∀ ω : SobolevFormSpaceZero k Ω,
        ‖SobolevFormSpaceZero.toLpForm (p := p) ω‖ ≤
          C * ‖SobolevFormSpaceZero.toLpForm (p := p) (extDerivativeZero ω)‖

end NitscheHodgeLaplacian.Foundations
