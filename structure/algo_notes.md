# Track "Algorithme" — upper bounds and structure for P3 : x a x² = c in S_N

Conventions: image tuples, product right to left, (pq)(j) = p(q(j)); x a x² applied to j is x(a(x(x(j)))).
Code: `algo_track.py` (numpy, exhaustive over all of S_N for N ≤ 9). Data: `algo_results.json`.
Every statement below is tagged **[Proved]**, **[Verified numerically]** (with the range) or **[Conjectured / open]**.

## 0. Two identities used throughout

**[Proved] (conjugation form).** x a x² = x (a x³) x⁻¹, hence c = x a x² ⟺ x⁻¹ c x = a x³. In particular
ctype(c) = ctype(a x³) for every solution x.

**[Proved] (first-moment identity).** Fix a of cycle type α and a class γ. Then
Σ_{c ∈ class γ} #{x : x a x² = c} = #{x ∈ S_N : ctype(a x³) = γ}.
Proof: the map x ↦ c(x) = x a x² sends S_N onto the multiset of targets; c(x) ∈ class γ iff a x³ ∈ class γ by the
conjugation form. ∎  The right-hand side is a class function of a and is computable in polynomial time from the
class algebra of S_N (number of cube roots per class × connection coefficients), so the *average* number of solutions
over a (ctype a, ctype c) block is polynomial-time computable. Checked exhaustively against the tables at N = 6, 7, 8
(121 + 225 + 484 blocks, all equal). It does *not* decide solvability (the solution count is ≈ Poisson(1) inside a block).

## 1. Special families of a

### 1.1 a an N-cycle
**[Proved] The "b then x" combination gives no gain.** Write b = x a x⁻¹. Since C(a) = ⟨a⟩ has order N, each N-cycle b
determines x up to N choices (x ∈ x₀⟨a⟩), and b ranges over the (N−1)! N-cycles: N · (N−1)! = N! candidates, i.e.
exactly all of S_N. The second constraint x³ = b⁻¹c is then the original equation. Equivalently, with a = τ (shift j ↦ j+1),
x is the "discrete exponential" j ↦ b^j(x₀) and x³ = b⁻¹c is a functional equation with no visible algebraic structure.

**[Proved, small remark] Lemma A5 is nearly empty here.** Solutions commuting with a are cube roots of a⁻¹c in ⟨a⟩ — at most
N candidates, checked in O(N²); all other solutions have b ≠ a.

**[Verified numerically, N ≤ 9] No cycle-type invariant decides the N-cycle case.** With a fixed as the N-cycle, the
solvable c form a union of ⟨a⟩-orbits, not a union of conjugacy classes: solvability of (a, c) is *not* determined by any
of I1 = ctype c, I2 = I1 + ctype(ac), I3 = I2 + ctype(a⁻¹c), I4 = cycle types of all words of length ≤ 4 in a^{±1}, c^{±1}
(length ≤ 3 at N = 8, 9), I5 = I4 + orbit structure of ⟨a, c⟩. Smallest witness (N = 6, a = (0 1 2 3 4 5)):
 - c_u = [5,3,2,1,4,0] = (0 5)(1 3): **0 solutions**;
 - c_s = [5,1,4,3,2,0] = (0 5)(2 4): **1 solution**;
 both have ctype c = (2,2,1,1), ctype ac = ctype a⁻¹c = (3,2,1), ctype a²c = (5,1), ctype ac² = (6), ctype acac = (3,1,1,1),
 ctype a²c² = (3,3), ctype [a,c] = (4,2), and ⟨a,c⟩ transitive in both cases (a is a 6-cycle). Verified by brute force
 over S_6 independently of the numpy pipeline.
Fraction of solvable c for a = N-cycle: 0.667, 0.639, 0.601, 0.597, 0.597 (N = 5…9) — indistinguishable from the
all-pairs fraction (0.708, 0.610, 0.587, 0.591, 0.591) and from 1 − 1/e = 0.632 at these N.

**[Verified numerically] Propagation-solver cost is the same as for random a.** Median nodes to exhaust the search on a
random target (8 trials, budget 10⁵): N = 10: 110 / 115 / 139; N = 14: 1911 / 2188 / 1936; N = 18: 38 224 / 39 608 / 38 829
for a = N-cycle / distinct-coprime cycle lengths / random a. Growth ≈ ×2.1 per unit of N (≈ 2^N over this range);
at N = 20 no run exhausts within 10⁵ nodes for any family. Planted instances are found much faster (median 4 469 / 10 752 /
5 547 nodes at N = 18) but without a family effect. (Data: `partD_solver_cost` in `algo_results.json`.)

**[Open]** Whether the N-cycle case is polynomial, or whether general P3 reduces to it. Neither direction was found.

### 1.2 a with distinct cycle lengths, a with pairwise-coprime cycle lengths
C(a) is then abelian (a product of cyclic groups), so Lemma A5 (b = a) is a tiny search, but the same negative results
hold: none of I0–I5 decides solvability for a ∈ {(5,1), (4,2), (3,2,1)} at N = 6, {(7),(6,1),(5,2),(4,3),(4,2,1)} at N = 7,
{(8),(7,1),(6,2),(5,3),(5,2,1),(4,3,1)} at N = 8 and all such types at N = 9 (`per_a.*.invariants_decide`). Solver cost
(family "coprime" above) is the same as for random a. **[Open]** polynomiality of these families.

**[Proved] The only a-families made trivial by centraliser size are the bounded-support ones (Lemma A3, N^{O(s)}) and
their mirror on x (below): the size of C(a) is irrelevant to the b ≠ a solutions.** (Small centraliser ⇒ few
Lemma-A5 solutions, but the solvable fraction stays ≈ 1 − 1/e, so almost all solutions are of the b ≠ a kind.)

## 2. The cycle-type enumeration route (fix λ = ctype x)

Write P3_λ for "is there a solution x of type λ". P3 = ∨_λ P3_λ over p(N) = e^{O(√N)} types, so a poly(N) algorithm for
P3_λ would give a 2^{O(√N)} algorithm for P3.

**[Proved] Reformulation.** P3_λ ⟺ ∃ g ∈ S_N : (g⁻¹ a g, g⁻¹ c g) = (a′, x_λ a′ x_λ²) for some a′, i.e. the simultaneous
conjugacy orbit of the pair (a, c) meets the graph {(a′, x_λ a′ x_λ²) : a′ ∈ class(a)} of the fixed bijection
a′ ↦ x_λ a′ x_λ². (Substitute x = g x_λ g⁻¹ in the conjugation form.) Simultaneous conjugacy of a *given* pair of pairs is
polynomial, but a′ is free; the trivial enumeration costs N!/z_λ per λ, and Σ_λ N!/z_λ = N!: fixing λ only divides the
search by z_λ, which is ≤ N for λ = (N). So the route needs a genuinely new per-λ idea; it is not a reduction in itself.

**[Proved] Polynomial λ's.**
 (i) λ with all parts ≤ 2 (x² = 1): the equation becomes x a = c, x = c a⁻¹; P3_λ ⟺ ctype(c a⁻¹) = λ. Decided by I3.
 (ii) λ with all parts in {1, 3} (x³ = 1): the equation becomes x a x⁻¹ = c; P3_λ ⟺ a and c are conjugate by an element of
     cycle type λ — a "conjugator of prescribed type" problem on the coset g C(a).
 (iii) λ with support s = N − m₁(λ) bounded: guess supp(x) and x on it (N^{O(s)}), the rest is x = id where x is trivial,
     i.e. c = a off the support — the mirror of Lemma A3.
**[Verified numerically, N = 7, 8, 9]** exactly the λ of type (i) are decided by I3 for every a; λ = (3,1^{N−3}) is
decided by I5 at N = 7, 8 but **not** at N = 9 (`first_failure_per_lambda`); every λ with a part ≥ 4, or with ≥ 2 parts
≥ 3, or of the form (3,2,…) fails I5 for some a already at N = 7 (e.g. λ = (3,2,1,1), (4,1,1,1), (3,2,2), (3,3,1)).
So the per-λ problem is not decided by short-word cycle types + orbit structure for any λ of large support.

**[Verified numerically]** Number of distinct λ realised by the solutions of a solvable pair (a = N-cycle, N = 8):
15 072 pairs with one λ, 6 644 with 2, 1 832 with 3, 452 with 4, 64 with 5, 8 with 6 — the per-λ subproblems are far from
disjoint, and a solvable instance typically has solutions of a single type (so guessing λ must be exhaustive).

**[Conjectured]** P3_λ is not polynomial for λ of unbounded support (in particular λ = (N)): the candidate count N!/z_λ and
the failure of all tested invariants are consistent with the per-λ problem retaining the full difficulty. No
2^{O(√N log N)} algorithm follows from this route; the m = 36 parameters are not threatened by it.

## 3. Exact structure of the solution set, N ≤ 9

Exhaustive over all x ∈ S_N for each cycle-type representative a (this covers every pair (a, c) up to simultaneous
conjugation): N = 3…9, 5 040 + 40 320 + 362 880 × (#types) products, cross-checked against pure-Python brute force at N = 5.

**[Verified numerically] Global statistics** (all pairs (a, c), weighted by class sizes):
| N | solvable fraction | P(1 sol) | P(2 sol) |
|---|---|---|---|
| 4 | 0.5833 | 0.266 | 0.250 |
| 5 | 0.7083 | 0.488 | 0.175 |
| 6 | 0.6097 | 0.339 | 0.179 |
| 7 | 0.5874 | 0.311 | 0.185 |
| 8 | 0.5911 | 0.338 | 0.159 |
| 9 | 0.5907 | 0.331 | 0.163 |
Poisson(1) would give 0.632 / 0.368 / 0.184; the excess of ≥ 3 solutions comes from pairs with large common centraliser
(e.g. a = c = id has N! solutions of the cube-root type — hist keys 9, 81, … at N = 4, 6).

**[Verified numerically] No polynomial cycle-type invariant of the pair decides solvability.** For each N the first
failure of each invariant (smallest witness) is recorded in `first_failure_per_invariant`:
 - I1 = ctype c fails from N = 3; I0 = orbit sizes of ⟨a,c⟩ fails from N = 3;
 - I2 = I1 + ctype(ac) fails from N = 5: a = (0 1 2 3 4), c_u = [4,3,2,0,1] (0 solutions) vs c_s = [4,1,3,0,2] (2 solutions),
   both ctype c = (4,1), same ctype(ac);
 - I3, I4 (all words of length ≤ 4), I5 (I4 + orbit structure) all fail from N = 6 with the a = 6-cycle witness of §1.1.
 At N = 7, 8, 9 all six invariants fail for every a except the degenerate a = id (I1 suffices: c must be a cube) and
 a = a transposition (I2 suffices; consistent with Lemma A3 at s = 2). At N = 6 the type (2,2,2) is also decided by I2.
The complete invariant is, tautologically, the simultaneous-conjugacy class of (a, c) (equivalently the isomorphism type of
the 2-coloured functional graph), which is an N!-sized family — the tables show that no coarsening of it by
short-word cycle types or orbit data survives already at N = 6, and the witness is inside the N-cycle family.

## 4. What this track establishes for the project

1. The "cheap" upper-bound ideas — centraliser size of a, cycle-type enumeration of x, cycle-type invariants of the pair —
   all fail provably or numerically; the only polynomial cases found are the known ones (bounded support of a or x,
   x² = 1, x³ = 1 ⇒ conjugacy with prescribed conjugator type) plus the exact first-moment formula of §0.
2. The special families a = N-cycle / distinct / coprime cycle lengths behave exactly like random a in solvable fraction,
   in solution-count distribution and in propagation-solver cost (N ≤ 20), and contain the smallest invariant-breaking
   witness (N = 6). If a hardness reduction is found, restricting a to an N-cycle is a natural target, not an obstacle.
3. Open: polynomiality of P3_λ for λ of unbounded support; polynomiality of the N-cycle case; any reduction between them.

## Deviations / scope
Exhaustive tabulation stopped at N = 9 (N = 10 feasible with the same code, ≈ 40 × 3 min, not run within the budget);
word length for I4 reduced from 4 to 3 at N = 8, 9; solver benchmark reached N = 20 only (loop planned to N = 22, stopped by the 200 s time guard), 8 trials per family, node budget 10⁵
(N = 20 already exceeds the budget for exhaustion; planted only). The first-moment identity was checked at N = 6–8 only.
