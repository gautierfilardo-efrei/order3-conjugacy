# CPWI single-round system — structural notes (voie A, Algèbre)

Convention. Permutations of {0,…,m−1} as image tuples, products right to left, (pq)(j)=p(q(j)).
Round: U' = U A V U, V' = V B U V with A,B ∈ A_m public. Single-round system
R(U',V';A,B) : unknowns (U,V) ∈ A_m×A_m. All statements below were checked by
`round_algebra.py` (own toolkit, stdlib only): exhaustively at m=5 (all 3600 (U,V) for 63
(A,B) pairs, all 3600 targets for the reductions), and on 3000 random instances each at
m=6 and m=7 (0 failures everywhere; numbers in `round_algebra_results.json`).

## A0. Primary normal form: one-variable equation x a x² = c in the wreath product

Let K = (A_m × A_m) ⋊ C₂ = A_m wr C₂ ≤ S_{2m}: (U,V) acts by U on {0,…,m−1} and by V on
{m,…,2m−1}; σ : j ↦ j+m mod 2m swaps the blocks. Put
  x = (U,V)·σ,   a = (B,A) (base group),   c = (U',V')·σ.

**Theorem A0.** (U,V) solves R(U',V';A,B) ⇔ x a x² = c in K. Every solution x ∈ K lies in the
σ-coset (the σ-exponent of x a x² is that of x), and x ↦ (U,V) is a bijection between the
solutions of the two problems.
*Proof.* Evaluate on j < m: σ sends j to block 1, x(j) = V(j)+m, x²(j) = UV(j), a x²(j) = BUV(j),
x a x²(j) = VBUV(j)+m = c(j) ⇔ V' = VBUV. On j = k+m: x(j) = U(k), x²(j) = VU(k)+m,
a x²(j) = AVU(k)+m, x a x²(j) = UAVU(k) = c(j) ⇔ U' = UAVU. ∎
Verified: all 3600 (U,V) × 21 (A,B) at m=5 (453 600 checks), 3000 random instances at m=6 and
m=7, zero failures; for all 7200 targets of two m=5 systems the K-solution set equals the image
of the round solutions and no base-group element solves (3600 × 50 targets checked).

Consequences (all verified in `part_K`):
- **Three positive occurrences of the unknown**, versus five with both signs in Lemma A2. The
  round problem is the fixed-shape equation x a x² = c (equivalently, with y = x⁻¹: y² a⁻¹ y = c⁻¹,
  the shape found for A = id in Part B) over a permutation group given as input.
- a = 1 ⇔ A = B = id: cube root x³ = c (Proposition A5 is exactly this case).
- x³ = (UVU, VUV)·σ and a x³ = x⁻¹ c x, hence the necessary condition ctype(a x³) = ctype(c).
  Since ctype((P,Q)σ) in S_{2m} is the doubled cycle type of PQ, this reads
  ctype(B·UVU·A·VUV) = ctype(U'V') on the pair (U,V): it keeps 32.8 % (m=5) / 27.4 % (m=6) of
  the pairs — same order as Σ_λ p(λ)² (Lemma A7), a polynomial factor.
- Base-group conjugation g = (g₁,g₂): x ↦ gxg⁻¹ realises (U,V) ↦ (g₁Ug₂⁻¹, g₂Vg₁⁻¹),
  a ↦ (g₁Bg₁⁻¹, g₂Ag₂⁻¹): this is Lemma A4 (difficulty depends only on the classes of A and B).
- **Relaxing to S_{2m} is not sound:** on 3 planted m=5 instances, x a x² = c had 2, 1, 1 solutions in
  the σ-coset of K (= the round solutions) but 2, 0, 2 additional solutions in S₁₀ outside
  S₅ wr C₂ (none in S₅ wr C₂ ∖ K). A solver working in S_{2m} must impose the block structure or
  filter its output.
- The twist obstruction of Proposition A6 is the constant a: for a ≠ 1 the equation is the
  "twisted cube root" x³ = a⁻¹·x⁻¹cx, not a power equation.

## A. Normal forms and invariants (elimination form)

**Lemma A1 (elimination).** V = A⁻¹U⁻¹U'U⁻¹ and, symmetrically, U = B⁻¹V⁻¹V'V⁻¹.
*Proof.* U' = UAVU ⇒ AVU = U⁻¹U' ⇒ V = A⁻¹U⁻¹U'U⁻¹. ∎

**Lemma A2 (single-unknown equation).** (U,V) solves R iff V is given by A1 and
  B U A⁻¹ U⁻¹ U' U⁻¹ = U U'⁻¹ U A V'.
U occurs 3 times, U⁻¹ twice: five occurrences, both signs, with the four constants
B, A⁻¹, U', U'⁻¹, A V' interleaved. *Proof.* Substitute A1 in V' = VBUV and
left-multiply by U U'⁻¹ U A. The converse holds because A1 recovers a V satisfying the first
equation and the second one is then equivalent to the displayed identity. ∎

**Lemma A3 (symmetric form).** (Secondary to Theorem A0, which is the primary normal form.) With Y = AV, Z = BU:
  Y (Z A⁻¹) Y = A V'  and  Z (Y B⁻¹) Z = B U'.
Also U'A = (UA) V (UA) and V'B = (VB) U (VB). (All verified; 226 800 exhaustive checks at m=5.)

**Lemma A4 (conjugation normal form — new).** For any g,h ∈ A_m,
  (U,V) solves R(U',V';A,B) ⇔ (Ug, g⁻¹V) solves R(U'g, g⁻¹V'; g⁻¹Ag, B)
                             ⇔ (h⁻¹U, Vh) solves R(h⁻¹U', V'h; A, h⁻¹Bh).
*Proof.* UAVU = (Ug)(g⁻¹Ag)(g⁻¹V)(Ug)g⁻¹ and VBUV = g(g⁻¹V)B(Ug)(g⁻¹V). ∎
Consequently the difficulty of R depends on (A,B) only through the conjugacy classes of A
and of B (taken in A_m; with g or h odd the statement holds in S_m, the unknowns moving to the
odd coset). Checked: solution counts coincide for all 3600 targets on 6 random (A,B,g) at
m=5 and 300 random instances at m=6. In particular: A ~ B reduces to A = B; A = B⁻¹ reduces to
A = B whenever A is A_m-conjugate to A⁻¹ (true for every A ∈ A_5 and A_6; false e.g. for a 7-cycle
in A_7, whose two A_7-classes are inverse to each other); (A,B) → (B,A) with (U,V) → (V,U) is the
swap symmetry.

**Proposition A5 (untwisted system reduces to cube roots; = case a = 1 of Theorem A0).** Let Y Z Y = c₁, Z Y Z = c₂ in a
group G. Then (YZ)³ = (YZY)(ZYZ) = c₁c₂ and (ZY)³ = c₂c₁. The solution set is exactly
  { (Y,Z) = (W⁻¹c₁, c₁⁻¹W²) : W ∈ G, W³ = c₁c₂ },
the map W ↦ (Y,Z) being injective; solutions exist iff c₁c₂ is a cube in G, and their number
equals the number of cube roots of c₁c₂ (no further consistency condition).
*Proof.* Given a solution put W = YZ; then WY = c₁ gives Y = W⁻¹c₁ and Z = Y⁻¹W = c₁⁻¹W².
Conversely, for W³ = c₁c₂: YZY = W⁻¹c₁c₁⁻¹W²W⁻¹c₁ = c₁ and ZYZ = c₁⁻¹W²W⁻¹c₁c₁⁻¹W² = c₁⁻¹W³ = c₂. ∎
Checked at m=5 for all 3600 targets (count = #cube roots in A_5, 0 mismatches; explicit set
checked on 200 targets): 2400/3600 targets have a solution; cube-root histogram of A_5 =
{0: 20 elements (the 3-cycles), 1: 39, 21: identity}. Existence/counting of k-th roots of a
permutation is polynomial (cycle-type arithmetic) and enumeration is output-polynomial, so the
untwisted system — which is exactly CPWI with A = B = id — is polynomial.

**Proposition A6 (the obstruction — why the CPWI shape does not reduce).** In the wreath form the obstruction is simply a ≠ 1
(x a x² = c is a cube root iff a = 1). In the (Y,Z) coordinates: put Z̃ = ZA⁻¹. Then
  Y Z̃ Y = c₁ := A V'   and   Z̃ (A Y B⁻¹) Z̃ = c₂' := B U' A⁻¹.
The second equation has inner letter AYB⁻¹, not Y. The identity behind A5 gives
(Y Z̃)³ = c₁ · (Z̃ Y Z̃), and Z̃ Y Z̃ is *not* determined by the data: it differs from c₂' by the
"twist defect" δ = A Y B⁻¹ Y⁻¹ (inserted inside the word). Precisely, let F(U',V') be the
cube-root family of A5 built from (c₁, c₂') and pulled back via V = A⁻¹Y, U = B⁻¹Z̃A. Then
  F ∩ Sol = Sol ∩ {V⁻¹AV = B} = F ∩ {V⁻¹AV = B}.
*Proof.* If (U,V) ∈ F ∩ Sol then Z̃(AYB⁻¹)Z̃ = c₂' = Z̃YZ̃, so AYB⁻¹ = Y, i.e. V⁻¹AV = B.
Conversely on the locus the two systems coincide. ∎ Symmetrically, untwisting with Ỹ = YB⁻¹
yields the locus {U⁻¹BU = A}. Numerically (m=5, 97 484 targets over 41 systems): the three sets
agree on every target; of 147 600 solutions only 10 080 (6.8 %) lie on the V-locus, 10 080 on the
U-locus, 4 054 on both; of 107 244 candidates produced by the cube-root construction only those
10 080 are solutions. When A is not conjugate to B the locus is empty: *no* solution is
reachable by the reduction. When A ~ B the reduction finds exactly the solutions with
V ∈ g·C_{A_m}(B) (g conjugating), a fraction |C(A)|/|A_m| of planted solutions — exponentially
small in m for generic A.
Non-invariance: for a random (A,B) at m=5 the solution count is not a function of the cycle
type of c₁c₂' (all 4 classes carry several counts) nor of (ctype U', ctype V', ctype U'V')
(36 of 37 classes carry several counts). Contrast with A5, where the count is a class function.

**Lemma A7 (cycle-type invariant).** With P = U'U⁻¹ = UAV and Q = V'V⁻¹ = VBU,
  PB = U (AQ) U⁻¹ and QA = V (BP) V⁻¹, hence ctype(U'U⁻¹B) = ctype(A V'V⁻¹).
(Both conditions are the same one, since PB ~ BP and AQ ~ QA.) Pruning power: as U ranges over
A_m, U'U⁻¹B is uniform on A_m, so the fraction of (U,V) pairs passing is Σ_λ p(λ)² over even
cycle types: 0.334 (m=5), 0.263 (6), 0.206 (7), 0.161 (8), 0.133 (9), 0.107 (10), 0.075 (12),
0.042 (16), 0.026 (20) — measured 0.334 / 0.263 at m=5,6. After elimination (V = V(U)), 34.7 %
(m=5) resp. 26.8 % (m=6) of the candidates U pass the test while 1/|A_m| are solutions. The
invariant is a polynomial (≈ m^{-1.8}) factor, useless against 2^{Θ(m log m)} search.

## B. Special cases (m=5 exhaustive, all 3600 targets per (A,B); frequencies exact)

| case | reduction | polynomial? | freq. (m=5 / 6 / 7) | frac. targets solvable | planted mean | max |
|---|---|---|---|---|---|---|
| generic (40 samples) | Lemma A2 (5 occurrences) | unknown | — | 0.665 | 1.92 | 9 |
| A = B = id | cube roots (A5), 0 mismatches | **yes** | 1/|A_m|² | 0.667 | 8.0 | 21 |
| A = id (or B = id) | X = U'⁻¹U, **3 occurrences**: X² (V'U') X = U'⁻¹BU' ⇔ X³(V'U') = X(U'⁻¹BU')X⁻¹; 0 mismatches on 28 800 targets (m=5) and 300 (m=6) | unknown ("twisted cube root"); yes iff V'U' = id or B = id | 1/60, 1/360, 1/2520 | 0.583 | 2.32 | 21 |
| A = B | swap symmetry Sol(U',V') ↔ Sol(V',U'); locus V ∈ C(A) (cube roots) covers 6.8 % of solutions; if moreover U'=V', diagonal solutions U=V ⇔ U A U² = U' (checked) | no reduction found | 1/60, 1/360, 1/2520 | 0.637 | 2.19 | 21 |
| A = B⁻¹ | = case A = B by A4 when A ~_{A_m} A⁻¹ (all A at m=5,6); stats identical to A = B | as A = B | same | 0.637 | 2.19 | 21 |
| A ~ B, A ≠ B | = case A = B by A4 (0 mismatches) | as A = B | 0.254, 0.183 (A_m); 0.334, 0.263, 0.206 (S_m) | 0.636 | 2.09 | 9 |
| A,B commute, non-trivial | none found | no | k(A_m)/|A_m| = 0.083, 0.019 | 0.628 | 2.08 | 9 |
| ord A = 2 | A = A⁻¹ only (A = B⁻¹ ⇔ A = B) | no | 0.25, 0.125, 0.042 | 0.708 | 1.74 | 4 |
| U' = V' (generic A,B) | none (swap needs A = B) | no | 1/|A_m| | 0.703 | 1.03 (random-target mean) | 9 |
| U',V' commute | none | no | k(A_m)/|A_m| | 0.684 | 1.02 (random-target mean) | 9 |

Reference solver (`round_attack.solve_round`, full enumeration, 25 planted instances per cell),
mean nodes m=6 / m=7: generic 64 / 171; A=B 66 / 196; A=id 49 / 115; A=B=id 79 / 326 (median
27 / 77, max 334 / 1340 — the tail is the enumeration of many cube roots); A,B commuting 58 / 165.
No special case is materially easier for the reference solver; A = id gives a ~30 % saving
consistent with the 5→3 occurrence drop.

**General obstruction (statement).** Write the system as Y Z̃ Y = c₁, Z̃ φ(Y) Z̃ = c₂' with
φ(Y) = A Y B⁻¹ (the two twists are the left action of A and the right action of B). Power-type
reductions (roots, conjugacy) require the two equations to be words in the *same* letters; this
holds iff φ(Y) = Y on the solution, i.e. on the locus V⁻¹AV = B (equivalently, after the other
untwisting, U⁻¹BU = A). Constants can be conjugated freely (Lemma A4) but the *relative* twist
— which sits between two occurrences of the unknown — cannot be removed by any substitution
Y ↦ gYh, Z ↦ g'Zh' with constant g,h,g',h' unless A and B are both trivial (up to the normal form
of A4: the only class pair for which the twist vanishes is (id,id)). Every polynomial case found
is therefore either (id,id) or lies on a sub-locus of exponentially small density; the
5-occurrence (generic) or 3-occurrence (A = id) single-unknown equations with interleaved
constants remain, and for these no root/conjugacy reduction exists.

## C. Solution statistics (number of preimages of a target)

Exact fact: Σ_targets #Sol = |A_m|², so the mean over uniformly random targets is exactly 1
for every (A,B); the planted (size-biased) mean equals 1 + Var(random).

| m | systems | random: Var, TV to Poisson(1) | P(no solution) (Poisson: 0.368) | planted: mean, Var, TV to 1+Poisson(1) | P(planted unique) (Poisson: 0.368) | max |
|---|---|---|---|---|---|---|
| 5 | all 3600 (A,B), 3600 targets each | 0.911, 0.057 | 0.335 | 1.911, 1.10, 0.064 | 0.423 | 21 (only (A,B,U',V') = (id,id,id,id), the 21 cube roots of id) |
| 6 | 50 random (A,B), 129 600 targets each | 1.016, 0.0013 | 0.369 | 2.016, 1.14, 0.0031 | 0.368 | 21 |
| 7 | 20 random (A,B), 6.35·10⁶ targets each (exhaustive via A_7 table) | 1.059, 0.011 | 0.376 | 2.059, 1.27, 0.016 | 0.361 | 72 (one system; P ≈ 1.4·10⁻⁷) |
| 7 | 20 planted, reference solver | — | — | histogram {1:6, 2:8, 3:4, 4:2} | 0.30 | 4 |

The random-target law is Poisson(1) up to TV ≈ 10⁻² already at m=6 and the planted law is
size-biased Poisson (1 + Poisson(1)); the excess variance at m=7 comes from a rare heavy tail
(targets with 10–72 preimages, total probability ~10⁻⁵) attached to structured (A,B). The
m=5 excess of unique planted solutions (0.423 vs 0.368) is a small-group effect. For the
security model: a planted target has ≥2 preimages with probability ≈ 0.63, so backward
inversion of n rounds is a tree with mean branching 1 (random targets) but ~2 along the
planted path; the expected tree size stays O(n).

## D. Checks of the identities stated in the task

All identities in the task statement are correct as written and were verified: elimination,
the 5-occurrence equation (3×U, 2×U⁻¹), the symmetric form, U'A = (UA)V(UA), V'B = (VB)U(VB),
(YZ)³ = c₁c₂ for the untwisted system, PB = U(AQ)U⁻¹, QA = V(BP)V⁻¹ and the cycle-type
invariant. Two nuances rather than errors: (1) "inverting n rounds is exactly n sequential
instances" presumes unique preimages — planted targets have ≈2 preimages on average, so the
backward search is a tree (mean total size O(n)); (2) the case list treats A = B, A = B⁻¹ and
A ~ B as distinct — by Lemma A4 they are the same problem (up to relabelling the targets)
whenever A ~_{A_m} A⁻¹, which is always the case at m = 5, 6 but fails e.g. for 7-cycles in A_7.

## Open questions raised

1. Complexity of the 3-occurrence equation x a x² = c (⇔ x³ = a⁻¹x⁻¹cx) in A_m wr C₂ — by Theorem A0 this IS the
   round problem — and of its A_m analogue X² C X = D (the A = id sub-case); both are "twisted cube root" problems; a polynomial algorithm for it
   would break the (id,B) family, a hardness proof would be the first for a fixed-shape
   equation over a growing group.
2. Whether the 5-occurrence equation of Lemma A2 admits a reduction to at most 3 occurrences
   for some other conjugacy-class pair (A,B) (Lemma A4 makes this a question about class pairs).
