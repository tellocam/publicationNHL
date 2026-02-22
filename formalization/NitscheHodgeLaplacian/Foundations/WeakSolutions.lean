import Mathlib.Tactic
import Mathlib.Analysis.InnerProductSpace.Basic

/-!
# Weak Solution Framework

This module provides continuous bilinear forms on arbitrary Banach spaces,
coercivity, weak solutions, Lax-Milgram theorem, inf-sup conditions, and
Brezzi's theorem for saddle-point problems.

## Mathematical objects

- `ContinuousBilinearForm V W` — bounded bilinear map V × W → ℝ
- `IsWeakSolution` — a(u,v) = f(v) for all v
- `InfSup` — Ladyženskaya-Babuška-Brezzi inf-sup condition
- `lax_milgram` — existence and uniqueness for coercive problems
- `brezzi_theorem` — existence and uniqueness for saddle-point problems

## References

- Brezzi, Fortin, "Mixed and Hybrid Finite Element Methods", Springer 1991.
  Chapter 1 for saddle-point theory.
- Ern, Guermond, "Theory and Practice of Finite Elements", Springer 2004.
  §2.1–§2.2 for abstract variational framework.
- Lax, Milgram (1954) for the original theorem.
- Tello Fachin 2025, §3.3 (BNB = Babuška-Nečas-Brezzi framework).

## Design

`ContinuousBilinearForm` is a structure (not an axiom) because the axiom system
needs to reason about specific instances (e.g., the Nitsche bilinear form).
The Lax-Milgram and Brezzi theorems are axiomatized as they require deep
functional analysis (closed range theorem, Banach open mapping).
-/

noncomputable section

namespace NitscheHodgeLaplacian.Foundations

/-! ## Continuous Bilinear Forms -/

/-- A continuous bilinear form on Banach spaces V × W → ℝ.
Generalizes the previous BilinearForm which was hardcoded to W₀^{1,2}.

Reference: Ern-Guermond §2.1; Brezzi-Fortin Chapter 1. -/
structure ContinuousBilinearForm (V : Type*) (W : Type*)
    [NormedAddCommGroup V] [NormedSpace ℝ V]
    [NormedAddCommGroup W] [NormedSpace ℝ W] where
  /-- The bilinear map. -/
  toFun : V → W → ℝ
  /-- Linearity in the first argument (additivity). -/
  map_add_left : ∀ u v w, toFun (u + v) w = toFun u w + toFun v w
  /-- Linearity in the second argument (additivity). -/
  map_add_right : ∀ u v w, toFun u (v + w) = toFun u v + toFun u w
  /-- Scalar multiplication in the first argument. -/
  map_smul_left : ∀ (c : ℝ) u v, toFun (c • u) v = c * toFun u v
  /-- Scalar multiplication in the second argument. -/
  map_smul_right : ∀ (c : ℝ) u v, toFun u (c • v) = c * toFun u v

/-! ## Properties of Bilinear Forms -/

/-- Boundedness: |a(u,v)| ≤ M ‖u‖ ‖v‖.

Reference: Ern-Guermond Definition 2.1. -/
def ContinuousBilinearForm.IsBounded {V W : Type*}
    [NormedAddCommGroup V] [NormedSpace ℝ V]
    [NormedAddCommGroup W] [NormedSpace ℝ W]
    (a : ContinuousBilinearForm V W) : Prop :=
  ∃ M : ℝ, M > 0 ∧ ∀ u v, |a.toFun u v| ≤ M * ‖u‖ * ‖v‖

/-- Coercivity: a(u,u) ≥ α ‖u‖² (for V = W).

Reference: Ern-Guermond Definition 2.2; Lax-Milgram hypothesis. -/
def ContinuousBilinearForm.IsCoercive {V : Type*}
    [NormedAddCommGroup V] [NormedSpace ℝ V]
    (a : ContinuousBilinearForm V V) : Prop :=
  ∃ α : ℝ, α > 0 ∧ ∀ u : V, a.toFun u u ≥ α * ‖u‖ ^ 2

/-! ## Weak Solutions -/

/-- Weak solution: a(u,v) = f(v) for all test functions v.

Reference: Evans §6.1; Ern-Guermond §2.1. -/
def IsWeakSolution {V : Type*}
    [NormedAddCommGroup V] [NormedSpace ℝ V]
    (a : ContinuousBilinearForm V V) (f : V →L[ℝ] ℝ) (u : V) : Prop :=
  ∀ v : V, a.toFun u v = f v

/-! ## Inf-Sup and Saddle-Point Theory -/

/-- Inf-sup condition (LBB): ∃ β > 0, ∀ q ∈ Q, ∃ v ∈ V \ {0},
b(v,q) ≥ β ‖v‖ ‖q‖.

Reference: Brezzi-Fortin §II.1; Ern-Guermond §2.2.
The LBB condition (Ladyženskaya-Babuška-Brezzi) is the discrete analogue of
the surjectivity condition for well-posedness of saddle-point problems. -/
def InfSup {V Q : Type*}
    [NormedAddCommGroup V] [NormedSpace ℝ V]
    [NormedAddCommGroup Q] [NormedSpace ℝ Q]
    (b : ContinuousBilinearForm V Q) : Prop :=
  ∃ β : ℝ, β > 0 ∧ ∀ q : Q,
    ∃ v : V, v ≠ 0 ∧ b.toFun v q ≥ β * ‖v‖ * ‖q‖

/-- Ellipticity on the kernel of b: a(v,v) ≥ α ‖v‖² for all v with b(v,q) = 0 ∀ q.

Reference: Brezzi-Fortin §II.1 (Brezzi's ellipticity-on-kernel condition). -/
def IsEllipticOnKernel {V Q : Type*}
    [NormedAddCommGroup V] [NormedSpace ℝ V]
    [NormedAddCommGroup Q] [NormedSpace ℝ Q]
    (a : ContinuousBilinearForm V V) (b : ContinuousBilinearForm V Q) : Prop :=
  ∃ α : ℝ, α > 0 ∧ ∀ v : V,
    (∀ q : Q, b.toFun v q = 0) → a.toFun v v ≥ α * ‖v‖ ^ 2

/-- Saddle-point solution: (u, p) ∈ V × Q satisfies
a(u,v) + b(v,p) = f(v) ∀ v,  and  b(u,q) = g(q) ∀ q.

Reference: Brezzi-Fortin §II.1; Ern-Guermond §2.2. -/
def IsSaddlePointSolution {V Q : Type*}
    [NormedAddCommGroup V] [NormedSpace ℝ V]
    [NormedAddCommGroup Q] [NormedSpace ℝ Q]
    (a : ContinuousBilinearForm V V) (b : ContinuousBilinearForm V Q)
    (f : V →L[ℝ] ℝ) (g : Q →L[ℝ] ℝ) (u : V) (p : Q) : Prop :=
  (∀ v : V, a.toFun u v + b.toFun v p = f v) ∧
  (∀ q : Q, b.toFun u q = g q)

/-! ## Existence and Uniqueness Theorems -/

/-- **Lax-Milgram theorem**: If a is a coercive, bounded bilinear form on a
Hilbert space and f is a bounded linear functional, then there exists a unique
weak solution.

Reference: Lax-Milgram (1954); Evans §6.2 Theorem 1; Brezis §6.2. -/
axiom lax_milgram {V : Type*}
    [NormedAddCommGroup V] [InnerProductSpace ℝ V] [CompleteSpace V]
    (a : ContinuousBilinearForm V V)
    (ha_coercive : a.IsCoercive) (ha_bounded : a.IsBounded)
    (f : V →L[ℝ] ℝ) :
    ∃! u : V, IsWeakSolution a f u

/-- **Brezzi theorem**: Existence and uniqueness for saddle-point problems.
Given Hilbert spaces V, Q, bounded bilinear forms a : V×V → ℝ and b : V×Q → ℝ,
if a is elliptic on ker(b) and b satisfies the inf-sup condition, then for any
f ∈ V* and g ∈ Q*, there exists a unique (u, p) ∈ V × Q satisfying the
saddle-point system.

Reference: Brezzi (1974); Brezzi-Fortin §II.1 Theorem 1.1;
Ern-Guermond §2.2 Theorem 2.21. -/
axiom brezzi_theorem {V Q : Type*}
    [NormedAddCommGroup V] [InnerProductSpace ℝ V] [CompleteSpace V]
    [NormedAddCommGroup Q] [InnerProductSpace ℝ Q] [CompleteSpace Q]
    (a : ContinuousBilinearForm V V) (b : ContinuousBilinearForm V Q)
    (ha_bounded : a.IsBounded) (hb_bounded : b.IsBounded)
    (ha_elliptic : IsEllipticOnKernel a b) (hb_infsup : InfSup b)
    (f : V →L[ℝ] ℝ) (g : Q →L[ℝ] ℝ) :
    ∃! up : V × Q, IsSaddlePointSolution a b f g up.1 up.2

end NitscheHodgeLaplacian.Foundations
