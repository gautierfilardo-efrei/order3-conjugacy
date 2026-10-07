# Track « quatre couches » — the M = ord_L(−2) layered padding Conj₃ → P₃

Files: `fourlayer_track.py` (instance generator, planted lifts, exhaustive bitmask-domain solver `P3Solver`, solution analysis),
`fourlayer_results.json` (all counts, witnesses, scan over exponent vectors, first-moment values, capped N = 40 runs),
`fourlayer_statement.tex` (proposition + proof + remark, ready for the paper), `scan_all_e.py`, `run_instance.py`,
`first_moment_mc.py` (drivers). Conventions as in the paper: image tuples, products right to left, x a x x applied to j is x(a(x(x(j)))).
Padded point (j, r) ∈ Ω × Z_M is encoded as r·|Ω| + j.

**Bottom line.** The construction is NOT a reduction Conj₃ → P₃: it already fails at k = 1 on the Conj₃‑NO instance
(a, c) = (5‑cycle, its inverse), m = 4, for the exponent vector e = (1,2,3,3) of padding_notes.md **and for every one of the 160
admissible exponent vectors**; the failure propagates to every k ≥ 1 by block products. The zeros observed at m = 2, 3 for
e = (1,2,3,3) are accidental, not structural: the exact first moment of the solution count over the class of 20‑cycles is
12337775/11639628 ≈ 1.06, and other admissible e give 5 or 20 solutions at m = 2, 3. No NP‑completeness theorem for P₃ follows
from this track; the proposition in `fourlayer_statement.tex` records the negative result.

## 1. PROVED

**P1 (counter‑example, k = 1).** L = 5, M = 4, ρ = (r ↦ r+1), e = (1,2,3,3), a = (0 1 2 3 4), c = a⁴ = a⁻¹.
Conj₃(a, c) is NO: every conjugator of a to a⁻¹ is a reflection z ↦ −z + t (block lemma), of order ≤ 2; equivalently the affine criterion
4³ = 64 ≡ 4 ≠ 1 (mod 5) (also brute‑forced over S₅: 0 order‑3 conjugators). P₃(a′, c′) is YES: with the encoding (j, r) ↦ 5r + j,
```
a' = (6,7,8,9,5, 12,13,14,10,11, 18,19,15,16,17, 3,4,0,1,2)
c' = (9,5,6,7,8, 13,14,10,11,12, 17,18,19,15,16, 2,3,4,0,1)
x  = (0,13,1,7,10, 14,9,15,5,4, 8,3,17,2,18, 16,6,11,12,19)
```
satisfies x a′ x² = c′ (20 direct evaluations; re‑checked with `verify_campaign.comp`). Cycle type of x: (15,3,1,1), x³ of type (5³,1⁵);
x does not respect the layers B × {r}. Exhaustive count: exactly 10 solutions (P3Solver, tree fully explored; word_solver.py with
5·10⁶ node budget, 208 694 nodes, budget not hit: 10 solutions, same set), forming two orbits of size 5 under conjugation by the common
centraliser ⟨a⁴ ⊗ id⟩ of a′ and c′ (a′ᵗ commutes with c′ iff 4 | t).

**P2 (every admissible exponent vector fails at k = 1).** The 160 vectors e ∈ {1,2,3,4}⁴ with Σe ≢ 0 and ⟨v, e⟩ ≢ 0 (mod 5),
v = (1,2,4,3) the left null vector of 2I + P_ρ — exactly the vectors for which Lemma P3(iii) of padding_notes.md has no blockwise‑power
solution — fall into 10 orbits of size 16 under layer rotation × unit scaling (instances in an orbit are simultaneously conjugate, so the
counts agree). Representatives and exhaustive counts (N = 20, all trees fully explored):

| e (rep.) | #sol m=2 | #sol m=3 | #sol m=4 |
|---|---|---|---|
| (1,1,1,3) | 0 | 0 | 30 |
| (1,1,1,4) | 0 | 0 | 10 |
| (1,1,2,2) | 0 | 0 | 5 |
| (1,1,2,4) ∋ (1,2,3,3) | 0 | 0 | 10 |
| (1,1,3,2) | 0 | 0 | 30 |
| (1,1,3,4) | 5 | 5 | 40 |
| (1,1,4,3) | 20 | 20 | 15 |
| (1,2,1,3) | 0 | 0 | 35 |
| (1,2,1,4) | 0 | 0 | 25 |
| (1,2,3,2) | 5 | 5 | 5 |

So m = 4 is spurious for every admissible e, and m = 2, 3 are spurious for 3 of the 10 orbits.

**P3 (failure for every k).** a = k disjoint 5‑cycles, c = a⁻¹ on B₁ and c = a on B₂,…,B_k (or c = a⁻¹ on every block). Conj₃ is NO by
the affine criterion (B₁ on a σ‑3‑cycle gives slope product −1; B₁ fixed gives (−1)³ = −1). Since a′ and c′ preserve every B_i × Z₄, the
disjoint union of the P1 witness on B₁ × Z₄ with the identity (c′ = a′ there) or with further P1 witnesses on the other blocks solves
x a′ x² = c′. For k = 2, m = (4,4) the 100 product solutions were verified explicitly (N = 40).

**P4 (exact first moment, k = 1).** With the paper's identity Σ_{c∈γ} #{x : x a x² = c} = #{x : ctype(a x³) = γ}, γ = (20), a a 20‑cycle,
whose character vanishes off the 20 hook partitions, and the hook generating function Σ_i χ^{(n−i,1^i)}(σ) tⁱ = Π_j (1 − (−t)^{ℓ_j})/(1+t):
mean number of solutions over the 19! twenty‑cycles c′ = (1/20!) Σ_h χ_h(1)⁻¹ Σ_x χ_h(x³) = **12337775/11639628 ≈ 1.0600**.
The formula was checked against brute force at N = 4, 5, 6 (5/3, 17/12, 7/5). Hence a padded NO instance is expected to have about one
solution; zeros are not structurally enforced.

**P5 (unpadded).** Already in S₅, (a, a⁻¹) is Conj₃‑NO and P₃‑YES (exactly 1 solution, brute force) — the padding is not needed to
separate the two problems on this instance.

## 2. VERIFIED (numerically)

- k = 1, e = (1,2,3,3): m = 2 → 0, m = 3 → 0 (reproduces padding_notes.md), m = 4 → 10, m = 1 → 5 (the lifts aᵗ ⊗ id, t = 0..4, all
  layer‑affine; the planted lift id is among them). All four trees fully explored by P3Solver and by word_solver (budgets not hit).
- P3Solver agrees with brute force on 40 random instances in S₅–S₇ (all solution sets identical, all trees complete).
- Monte‑Carlo first moment: k = 1: 1.051 ± 0.007 (4·10⁵ samples; exact 1.060); k = 2, class (20,20): 1.03 ± 0.02 (1.5·10⁶ samples).
- The solution x of P1 and the 100 product solutions at k = 2 were re‑evaluated with the independent helpers of verify_campaign.py.

## 3. OBSTRUCTIONS / FAILED

1. **The reduction is refuted** (P1–P3): non‑layered solutions exist on Conj₃‑NO instances for all k, for every admissible exponent vector,
   with L = 5 (the only prime L ≤ 50 with ord_L(−2) = 4; L = 15 is the only other L ≤ 50). The mechanism is not an artefact of a special e:
   the spurious solutions have assorted cycle types ((15,3,1,1), (8,7,3,2), (17,3), (5,5,5,5), …) and are not affine on layers. The
   heuristic reason is P4: over a conjugacy class the mean count is ≈ 1, so correctness on all NO instances would need a structural
   vanishing that the layered padding does not provide (contrast Conj₃, where x³ = 1 forces the block/affine structure).
2. **N = 40 (k = 2) exhaustion failed within budget**: six instances m ∈ {(2,3),(2,2),(4,2),(1,2),(4,4),(1,1)}, e = (1,2,3,3), 1500 s each
   in parallel: 38k–65k nodes, no tree completed (~38 ms/node). Only (1,1) produced solutions in the explored part (10, all layer‑affine,
   planted lift found); for (4,4) the search did not even reach the 100 product solutions that provably exist. **All k = 2 zero counts are
   inconclusive** (lower bounds). The question "are the (2,3)/(2,2) padded instances solution‑free?" remains open but is moot for the
   reduction, which fails anyway by P3.
3. **Not run (time):** the next M: L = 11 (M = ord₁₁(−2) = 5, N = 55 at k = 1), L = 7 (M = 6, N = 42; L ≡ 1 mod 3 is outside the family),
   L = 13 (M = 12), L = 17 (M = 8). By P4's argument there is no reason to expect them to behave differently; the k = 1, L = 5 scan over all
   admissible e is the decisive negative evidence for the whole family of M = ord_L(−2) layered paddings.

## 4. OPEN

- Are the k = 2 padded instances m = (2,3), (2,2) solution‑free? (needs a solver ~50× faster or symmetry reduction; moot for the reduction.)
- Is there any padding (a, c) ↦ (a′, c′) with a′ ∼ c′ whose spurious‑solution count vanishes structurally on the hard family? P4 suggests that
  any construction must break the "random‑target" behaviour, i.e. must make a′x³ ∼ c′ impossible for all x not conjugating a′ to c′; the
  layered paddings do not.
- Whether a single P₃ equation can simulate the system {x a x² = c, x a² x² = c²} (Lemma P2 of padding_notes.md) — unchanged.

## 5. Solver notes

`P3Solver(a, c).solve(max_solutions=None, node_budget=None, time_cap=None)` → (solutions, stats); `stats["complete"]` is True iff the whole
search tree was explored (exhaustive list). Propagation: all‑different (naked + hidden singles) and generalised arc consistency on the three
positions of the chain p → x(p) → x²(p) → a′x²(p) → x(a′x²(p)) = c′(p), with bitmask domains; branching on the smallest domain.
Throughput: N = 20 ≈ 2.7 ms/node (≈ 2000 nodes per instance), N = 40 ≈ 38 ms/node. Caveat on word_solver.py: its `st["exhausted"]` flag
means the node budget was hit (search incomplete), the opposite of "tree exhausted".
