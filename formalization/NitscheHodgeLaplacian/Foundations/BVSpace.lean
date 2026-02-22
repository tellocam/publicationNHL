import Mathlib.Tactic
import NitscheHodgeLaplacian.Foundations.Domains
import NitscheHodgeLaplacian.Foundations.LpSpaces
import NitscheHodgeLaplacian.Foundations.Sobolev

/-!
# BV Spaces (Functions of Bounded Variation)

This module axiomatizes BV(Ω) — the space of functions of bounded variation —
with embeddings into L^p spaces and from W^{1,1}.

## Mathematical objects

- `BVSpace Ω` — space of functions of bounded variation
- `BVSpace.toLp` — BV(Ω) ↪ L^1(Ω)
- `BVSpace.toLp_sharp` — BV(Ω) ↪ L^q(Ω) for q ≤ n/(n-1) (sharp Sobolev embedding)
- `SobolevSpace.toBV` — W^{1,1}(Ω) ↪ BV(Ω)

## References

- Evans, Gariepy, "Measure Theory and Fine Properties of Functions",
  CRC Press, 2015. Chapter 5 for BV spaces.
- Ambrosio, Fusco, Pallara, "Functions of Bounded Variation and Free Discontinuity
  Problems", Oxford 2000. Part I.
- Adams-Fournier §4.12 for the sharp Sobolev embedding.

## Design

BVSpace is not used directly in the Nitsche Hodge-Laplacian error bounds but
is included for completeness of the PDE foundation library. The sharp embedding
exponent n/(n-1) is the critical Sobolev exponent for W^{1,1} ↪ L^{n/(n-1)}.
-/

noncomputable section

open scoped NNReal ENNReal

namespace NitscheHodgeLaplacian.Foundations

/-! ## BV Space -/

/-- BV(Ω): space of functions of bounded variation on a bounded domain.

Reference: Evans-Gariepy Chapter 5; Ambrosio-Fusco-Pallara Part I. -/
opaque BVSpace {n : ℕ} (Ω : BoundedDomain n) : Type

/-- BVSpace is a normed add comm group. -/
axiom BVSpace.instNormedAddCommGroup {n : ℕ} {Ω : BoundedDomain n} :
    NormedAddCommGroup (BVSpace Ω)

attribute [instance] BVSpace.instNormedAddCommGroup

/-- BVSpace is a normed space over ℝ. -/
axiom BVSpace.instNormedSpace {n : ℕ} {Ω : BoundedDomain n} :
    NormedSpace ℝ (BVSpace Ω)

attribute [instance] BVSpace.instNormedSpace

/-- BVSpace is complete (Banach).

Reference: Ambrosio-Fusco-Pallara Theorem 3.23. -/
axiom BVSpace.instCompleteSpace {n : ℕ} {Ω : BoundedDomain n} :
    CompleteSpace (BVSpace Ω)

attribute [instance] BVSpace.instCompleteSpace

/-! ## Embeddings -/

/-- BV(Ω) ↪ L^1(Ω): continuous embedding.

Reference: Evans-Gariepy §5.1 Theorem 1. -/
axiom BVSpace.toLp {n : ℕ} {Ω : BoundedDomain n} :
    BVSpace Ω →L[ℝ] LpSpace 1 Ω

/-- BV(Ω) ↪ L^q(Ω) for q ≤ n/(n-1) when n ≥ 2 (sharp Sobolev-type embedding).

Reference: Evans-Gariepy §5.6 Theorem 1 (Gagliardo-Nirenberg-Sobolev for BV). -/
axiom BVSpace.toLp_sharp {n : ℕ} {Ω : BoundedDomain n} (hn : n ≥ 2)
    (q : ℝ≥0∞) (hq : q ≤ n / (n - 1)) :
    BVSpace Ω →L[ℝ] LpSpace q Ω

/-- W^{1,1}(Ω) ↪ BV(Ω): Sobolev functions are BV.

Reference: Evans-Gariepy §5.1; the distributional gradient of W^{1,1} functions
is absolutely continuous, hence has bounded variation. -/
axiom SobolevSpace.toBV {n : ℕ} {Ω : BoundedDomain n} :
    SobolevSpace 1 1 Ω →L[ℝ] BVSpace Ω

end NitscheHodgeLaplacian.Foundations
