# Track "Ordre 3" — the order-3 core D0 of P3, and rigidity

Conventions: permutations of [N] = range(N) as image tuples, product right to left,
(pq)(j) = p(q(j)). P3(a,c): ∃x, x a x² = c.  D0(a,c): ∃x, x³ = 1 and x a x⁻¹ = c.
Since x⁻¹ = x² when x³ = 1, **D0 = P3 restricted to solutions of order dividing 3**,
and D0(a,c) ⇔ the conjugator coset {x : x a x⁻¹ = c} = g·C(a) contains an element of order | 3.

Status legend: **[PROVED]** = complete proof below; **[VERIFIED]** = checked by code, see
`order3_results.json`; **[CONJ]** = conjecture / open.

---

## 1. Literature (exact statements, DOIs)

| Ref | Statement used | Identifier | status |
|---|---|---|---|
| Lubiw 1981, *Some NP-complete problems similar to graph isomorphism*, SIAM J. Comput. 10(1) 11–21 | deciding whether a graph has a fixed-point-free automorphism is NP-complete; equivalently FPF element in a permutation group given by generators. In her reduction every FPF automorphism is an involution. | doi:10.1137/0210002 | DOI verified |
| Buchheim–Cameron–Wu 2009, *On the subgroup distance problem*, Discrete Math. 309(4) 962–968 | Hamming subgroup-distance problem (∃h∈H with d(π,h) ≤ K) is NP-complete even for H abelian of exponent 2 (Thm 1) | doi:10.1016/j.disc.2008.01.036 | DOI verified |
| Cameron–Wu 2010, *The complexity of the weight problem for permutation and matrix groups*, Discrete Math. 310(3) 408–416 | weight problem (∃g∈G with w(g)=k) NP-complete for Hamming and most metrics; FPF-element problem NP-complete even for elementary abelian 2-groups with small orbits (Sect. 4) | doi:10.1016/j.disc.2009.03.005 | DOI verified |
| Lohrey–Rosowski–Zetzsche 2022/2025, *Membership problems in finite groups*, MFCS 2022 / J. Algebra | membership in a product of three cyclic permutation groups is NP-complete (knapsack, n ≥ 3 fixed) | doi:10.4230/LIPIcs.MFCS.2022.71 ; 10.1016/j.jalgebra.2025.03.011 | verified |
| Lohrey–Rosowski 2025, *Finding cycle types in permutation groups with few generators* | "∃ element of prescribed cycle type in ⟨S⟩" is NP-complete even for 2-generated abelian groups; logspace for cyclic groups | arXiv:2503.02864 ; doi:10.1007/978-981-95-0218-9_27 | verified (abstract) |
| Garey–Johnson 1979, [SP16] Numerical 3-Dimensional Matching | strongly NP-complete | book | standard, not re-checked online |
| Breusch 1932, Math. Z. 34, 505–526 | for n ≥ 7 there is a prime ≡ 2 (mod 3) in (n, 2n) | — | **UNVERIFIED** (from memory; only used to bound the prime L in the reduction) |

None of these results applies directly to D0: they all concern groups **given by generators**
with tailored actions, whereas C(a) is a fixed, very structured group (∏_L C_L wr S_{k_L}).
No paper on "element of order 3 in a coset of a centralizer" was found (consistent with the
87-reference survey of the memo).

---

## 2. The block lemma  [PROVED, VERIFIED]

**Setting (block instances).** Fix L ≥ 2 and k ≥ 1, N = kL, blocks A_i = {iL, …, iL+L−1} ≅ Z/L.
Let a act on each block as the translation s : z ↦ z+1 (a is a product of k disjoint L-cycles),
and let c act on each block A_i as an L-cycle c_i = π_i s π_i⁻¹ with π_i ∈ Sym(Z/L).
Write C = ⟨s⟩ ≤ Sym(Z/L).

**Lemma B.** D0(a,c) holds iff there is σ ∈ S_k with σ³ = 1 such that

* for every fixed point i of σ: Fix(i) :⇔ ∃r, (π_i s^r)³ = 1;
* for every 3-cycle i → j → k → i of σ: Adm(i,j,k) :⇔ ∃ r_i, r_j, r_k with
  π_i s^{r_k} π_k s^{r_j} π_j s^{r_i} = 1 in Sym(Z/L)  (equivalently 1 ∈ C π_k C π_j C π_i).

*Proof.* Let x a x⁻¹ = c. x maps a-cycles onto c-cycles; both families are the blocks, so
x(A_i) = A_{σ(i)} for some σ ∈ S_k. In block coordinates x_i := x|A_i : Z/L → Z/L satisfies
x_i s x_i⁻¹ = c_{σ(i)} = π_{σ(i)} s π_{σ(i)}⁻¹, hence π_{σ(i)}⁻¹ x_i ∈ C(s) = C, i.e.
x_i = π_{σ(i)} s^{r_i}. Now x³ = 1 ⇔ σ³ = 1 and, for every i, x_{σ²(i)} x_{σ(i)} x_i = id_{Z/L}.
For a fixed point this reads (π_i s^{r})³ = 1; for a 3-cycle i→j→k it reads
π_i s^{r_k} π_k s^{r_j} π_j s^{r_i} = 1 (the three conditions around a 3-cycle are conjugate to
each other, so one suffices). Conversely any σ, (r_i) satisfying the conditions define, block by
block, a permutation x with x a x⁻¹ = c and x³ = 1. ∎

*Verification.* Exhaustive over all π-tuples for (L,k) ∈ {(2,3),(3,2),(3,3),(2,4)} (N ≤ 9,
brute force over all order-3 elements of S_N) and, with the propagation solver, 300 random
π-tuples for L = 5, k ∈ {3,4} — criterion and brute force agree on every instance
(`block_lemma_exhaustive`, `block_lemma_L5_random`).

---

## 3. Affine blocks  [PROVED, VERIFIED]

Take L prime and π_i(z) = m_i z with m_i ∈ (Z/L)^×. Then c_i = π_i s π_i⁻¹ is the translation
z ↦ z + m_i, i.e. **c|A_i = a^{m_i}** — both a and c are block-wise translations and c ∈ C(a).

**Lemma C.** For affine blocks, Adm(i,j,k) ⇔ m_i m_j m_k = 1 (mod L) and Fix(i) ⇔ m_i³ = 1 (mod L).

*Proof.* C π C = {z ↦ m z + b : b ∈ Z/L} is the full set of affine maps with multiplier m, so
C π_k C π_j C π_i C is the set of affine maps with multiplier m_k m_j m_i; it contains the
identity iff m_i m_j m_k = 1. For f(z) = m z + b, f³(z) = m³ z + b(1+m+m²); if m³ = 1 then f³ = id
(b = 0 works, and any b if m ≠ 1); if m³ ≠ 1, f³ ≠ id for every b. ∎

Hence, with L prime and **3 ∤ L−1** (L ≡ 2 mod 3, L ≥ 5), the only unit with m³ = 1 is m = 1, and:

**Corollary.** For L prime ≡ 2 (mod 3), m_i ∈ (Z/L)^× \ {1}:
D0(a,c) ⇔ [k] can be partitioned into triples {i,j,k} with m_i m_j m_k = 1 (mod L).

*Verification.* Exhaustive over all m-tuples for (L,k) ∈ {(5,3),(7,3),(5,4),(7,4),(11,3)}
(N up to 33, 2 832 instances), criterion vs. propagation solver, all agree
(`affine_criterion_vs_solver`). (L = 7 has 3 | L−1 and exercises Fix(i) with m ∈ {2,4}.)

---

## 4. Theorem: D0 is NP-complete  [PROVED; reduction VERIFIED on 250 instances]

**Theorem 1.** D0 = {(a,c) ∈ S_N² : ∃x, x³ = 1, x a x⁻¹ = c} is NP-complete. It remains
NP-complete when a is a product of k disjoint L-cycles (L prime) and c = a^{m_i} on the i-th
block, i.e. c ∈ C(a); and the certificate x may be required to be fixed-point-free of cycle type
3^{N/3}.

*Proof.* Membership in NP is clear. Hardness: reduce from **Numerical 3-Dimensional Matching**
(N3DM, strongly NP-complete [Garey–Johnson SP16]): disjoint W, X, Y of size q, sizes
s(·) ∈ Z_{>0}, bound B; question: partition into q triples (w,x,y) with s(w)+s(x)+s(y) = B.
W.l.o.g. all sizes lie in [1, B−2] (otherwise the instance is trivially NO; output a fixed NO
instance of D0 such as a = (0 1), c = (2 3) in S_4, which has no order-3 conjugator — checked).

Choose a prime L ≡ 2 (mod 3) with R := L − 1 ≥ 10B (exists in [10B, 20B] by Breusch's
theorem — see caveat in §1 — and is found by trial division in poly(B) = poly(input) time
because N3DM is strongly NP-complete, so B may be assumed polynomially bounded).
Let ω be a primitive root mod L and K := 2B. Define exponents
  e_w = s(w),  e_x = s(x) + K,  e_y = s(y) + R − K − B   (mod R),
m_i := ω^{e_i}, k := 3q blocks, N := 3qL, a := block translation by 1, c := block translation by m_i.

*Correctness.* By the Corollary, D0(a,c) ⇔ [k] splits into triples with e_i + e_j + e_k ≡ 0 (mod R).
Let S = sum of the three sizes ∈ [3, 3B−6]. Going through the 10 class patterns:
WWW: S; WWX: S+K; WXX: S+2K; XXX: S+3K; XXY: S+K−B; WWY: S−K−B; WYY: S−2K−2B; XYY: S−K−2B;
YYY: S−3K−3B; WXY: S − B. With K = 2B and R ≥ 10B every pattern except WXY lies in an open
interval avoiding all multiples of R (e.g. WWY ∈ (−3B, 0), YYY ∈ (−9B, −6B), XXX ∈ (6B, 9B)), while
WXY ≡ 0 ⇔ S = B because |S − B| < 3B < R. So the admissible triples are exactly the W×X×Y triples
of size-sum B, and D0(a,c) ⇔ N3DM is YES. Moreover e_i ≢ 0 (mod R) for every i
(e_w ∈ [1,B), e_x ∈ [2B+1, 3B), e_y ∈ (R−3B, R−2B)), so m_i ≠ 1, Fix(i) never holds, σ is fixed-point-free
and every certificate x moves every point: ctype(x) = 3^{N/3}. The map is polynomial (N = 3qL = O(qB)). ∎

*Verification.* `reduce_n3dm` implemented; 250 random/planted N3DM instances with q ∈ {2,3,4},
B ∈ {6,7,8,9}, (q,B,L) ∈ {(2,6,71),(2,8,83),(3,7,71),(3,9,101),(4,8,83)} (N ∈ {426,498,639,909,996}, max N = 996): brute-force N3DM answer equals the D0 criterion on every
instance (88 YES / 162 NO). Ten planted YES instances (N = 426) were additionally solved by the
generic point-level propagation solver, which found an order-3 conjugator in 76 nodes each.
(Deviation: the point-level solver could not decide the NO instances at N > 300 within 10 min — the
first attempt was aborted; NO instances are certified through the *proved* Lemmas B–C only.)

**Remarks.**
1. Two polynomial one-variable equations (conjugacy x a x⁻¹ = c and the root x³ = 1) intersect in an
   NP-complete system; equivalently, P3 with the side condition x³ = 1 — or with ctype(x) = 3^{N/3}
   prescribed — is NP-complete. This is the "coset with prescribed cycle type" phenomenon
   (Lohrey–Rosowski) transported to *centralizer* cosets of the very rigid form g·(C_L wr S_k).
2. The hard instances have c ∈ C(a) (c = a^{m_i} block-wise) — the problem is hard even when a and c
   commute.
3. The block structure could be replaced by any a whose cycles all have the same length L; the
   reduction only uses the S_k part of C(a) = C_L wr S_k and the affine normaliser of C_L.

---

## 5. Rigidity: does D0-hardness transfer to P3?  [VERIFIED negative for the hard family; census]

P3 solutions satisfy x a x⁻¹ = b := c x⁻³ with b ~ a. Rigidity of (a,c) := every P3 solution has x³ = 1
(then P3(a,c) = D0(a,c)).

**Fact R1 [PROVED].** The affine block family of Theorem 1 is *never* rigid, and P3 is trivially YES on
it: since L is prime ≠ 3, t_i := (m_i − 1)/3 mod L exists and x := block-wise translation by t_i
satisfies x ∈ C(a), x³ = a⁻¹c, hence x a x² = a x³ = c (Lemma A5 solution) while x³ ≠ 1.
Checked on (L,m) = (5,(2,2,4)) (`affine_L5_commuting_solution`). So **the chain "D0 hard + rigidity ⇒
P3 hard" does not close with this family.** Killing the commuting solutions requires a⁻¹c to have no
cube root in C(a) = C_L wr S_k, i.e. 3 | L and m_i ≢ 1 (mod 3) and no three blocks with equal m_i —
but then L is not prime, and non-commuting solutions (b ≠ a, x³ ≠ 1) are not controlled either.

**Census [VERIFIED], N ≤ 8** (all a up to conjugacy, all c, all x ∈ S_N; `rigidity_census`):

| N | solvable (a,c) | rigid | rigid with c ≠ a | rigid with unique solution | pairs with x³ constant |
|---|---|---|---|---|---|
| 6 | 4 594 | 242 | 232 | — | 2 793 |
| 7 | 43 138 | 1 140 | 1 128 | 579 | 25 211 |
| 8 | 505 999 | 4 038 | 4 020 | 2 647 | 302 573 |

Rigid pairs are ≈ 0.8 % of solvable pairs at N = 8 and about two thirds of them are "accidental"
(a unique solution which happens to have order 3); no cycle type of a has rigid fraction above 5.2 %
(N=7 best: (3,2,1,1) 5.2 %, (4,3) 4.8 %, (7) 4.7 %). "x³ constant across all solutions" is common
(≈ 60 %) but this is again the unique-solution effect (Poisson(1) solution counts).

**Observation R2 [VERIFIED, small sample] → [CONJ].** For block instances with **L = 3** (a of cycle
type 3^k, blocks Z/3, arbitrary π_i ∈ S_3), a random sample of 30 instances at k = 3 (N = 9, exhaustive
P3) gave 19 solvable, of which 17 rigid and 2 with *no* order-3 solution — and **no mixed instance**
(P3 solutions either all of order 3 or none of order 3). Conjecture R2: for a of type 3^k and c an
arbitrary 3-cycle on each block, the set of P3 solutions is either contained in {x³=1} or disjoint from
it. If true, P3 restricted to the L = 3 block family would be D0 on that family plus an extra
decision; but Lemma C needs L prime with (Z/L)^× ≠ {1,2}, so L = 3 gives no hardness by itself. This is
the natural next experiment (k = 4, N = 12 with the word solver; then a structural proof).

---

## 6. Summary of what is established

* **Proved:** Lemma B (block lemma), Lemma C (affine blocks), Theorem 1: **D0 is NP-complete**, already for
  a = k disjoint L-cycles and c = a^{m_i} block-wise (c ∈ C(a)), with certificates of cycle type 3^{N/3}.
  Hence "P3 with x³ = 1" / "P3 with ctype(x) prescribed" is NP-complete, and the conjugator coset
  g·C(a) — a coset of C_L wr S_k — can be NP-hard to search for an element of order 3. Fact R1.
* **Verified numerically:** solver ≡ exhaustive (360 pairs, N ≤ 8); Lemma B (276 exhaustive + 300 solver
  instances); Lemma C (2 832 instances, N ≤ 33); reduction (250 N3DM instances, N ≤ 996); rigidity
  census N ≤ 8; Observation R2 (30 instances).
* **Conjectured / failed:** rigidity for the hard family fails (R1), so P3 NP-hardness is **not**
  obtained on this track; Conjecture R2 (L = 3 dichotomy) is the surviving lead. Breusch 1932 citation
  unverified (not load-bearing: the prime L can be searched for directly).
