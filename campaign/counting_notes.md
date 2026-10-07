# Counting track — P3 : x a x² = c in Sym(N)

Conventions: image tuples on {0..N-1}, (pq)(j) = p(q(j)); the word x a x x applied to j is x(a(x(x(j)))).
Code: `counting_track.py` (numpy, exhaustive over all x ∈ S_N for every cycle-type representative a, N ≤ 9).
Data: `counting_results.json` (phases A–D). Every statement is tagged PROVED / VERIFIED / FAILED / OPEN.
Notation: #Sol(a,c) = #{x : x a x² = c}; α = ctype(a), γ = ctype(c); C(a) the centraliser; ρ₃(ν) = number of cube
roots of a permutation of cycle type ν; z_λ = |C(x)| for x of type λ.

## A. First moment — PROVED (closed formula) and VERIFIED (N = 3..9, all 1813 (α,γ) blocks)

**Theorem A.1 (first-moment identity, PROVED).** For a of type α and any class γ,
  Σ_{c ∈ γ} #Sol(a,c) = #{x ∈ S_N : ctype(a x³) = γ} = Σ_ν ρ₃(ν) · K(γ, ν; α),
where K(γ,ν;α) = #{(z,w) ∈ γ × ν : z w = a} = (|γ||ν|/N!) Σ_χ χ(γ)χ(ν)χ(α)/χ(1) (Frobenius; S_N-characters are real).

*Proof.* x a x² = x (a x³) x⁻¹, so c = x a x² lies in γ iff a x³ ∈ γ; summing #Sol(a,c) over c ∈ γ counts x with
ctype(a x³) = γ. Group these x by y = x³: #{x : x³ = y} = ρ₃(ctype y) is a class function, and for fixed a the number of
y of type ν with a y ∈ γ equals #{(z,w) ∈ γ×ν : z w = a} via z = a y, w = y⁻¹ (ctype(y⁻¹) = ctype(y)). ∎

**Lemma A.2 (cube roots, PROVED; VERIFIED by brute force N ≤ 7).** If ν has m_L cycles of length L,
  ρ₃(ν) = Π_L Σ_{k} m_L! / ((m_L − 3k)! k! 6^k) · (2L²)^k,
the sum over 0 ≤ k ≤ m_L/3, restricted to 3k = m_L when 3 | L.
*Proof.* x³ = y: an x-cycle of length M with 3 ∤ M cubes to one M-cycle, one with 3 | M cubes to three M/3-cycles. Hence
the L-cycles of y are partitioned into singletons (3 ∤ L only; on an L-cycle y₀ the root is the unique power y₀^t,
3t ≡ 1 mod L, since x commutes with y and preserves the cycle) and unordered triples merged into a 3L-cycle of x
(2 cyclic orders × L × L choices of the interleaving). k triples among m_L cycles: m_L!/((m_L−3k)! k! 6^k) ways. ∎

The right-hand side of A.1 is computable in time polynomial in the number of classes used (here by Murnaghan–Nakayama;
p(N) = e^{O(√N)} terms). Agreement with the exhaustive tables: N = 3,…,9: 9, 25, 49, 121, 225, 484, 900 blocks, all equal
(`phaseA_first_moment`). The N = 9 tables were also compared block-by-block (900 solution-count histograms) with the
independent pipeline `algo_results_N9.json`: 0 mismatches.

**Consequence (PROVED).** The *average* of #Sol over an (α,γ) block is polynomial-time computable (up to p(N) character
values), but the identity says nothing about an individual pair: inside a block the distribution is ≈ Poisson(mean).

## B. Second moment — identity PROVED (reformulation only), no closed formula found; VERIFIED numerics

**Proposition B.1 (PROVED, elementary).** M₂(α,γ) := Σ_{c∈γ} #Sol(a,c)² = #{(x,y) ∈ S_N² : x a x² = y a y², ctype(a x³) = γ}
= Σ_{x : ctype(a x³)=γ} #Sol(a, x a x²), a class function of a. (Σ over γ: #{(x,y) : x a x² = y a y²} = #solutions of a
two-variable equation with constant a.)
Substituting y = x u turns x a x² = y a y² into a x² = u a (x u)², which is not a product of conjugacy/power constraints;
no class-algebra expression was found (FAILED/OPEN).

**Numerics (VERIFIED, N = 4..9, all α).** `phaseB_second_moment` stores per block n = |γ|, M₁, M₂. Prediction under
"Poisson(M₁/n) inside each block": M₂^pred = M₁ + M₁²/n.
- Total over γ, a = N-cycle: M₂/N! = 1.967, 2.139, 2.224, 2.217 (N = 6..9); ratio M₂/Σ_γ M₂^pred = 0.955, 1.049, 1.099, 1.099.
  Over all α the ratio ranges in [0.89,1.07] (N=6), [0.97,1.33] (7), [0.98,1.48] (8), [0.99,1.75] (9).
- Block level: only ≈ 48–50 % of the non-trivial blocks are within 10 % of the Poisson prediction; the excess is
  concentrated on small blocks where c is close to a (large common centraliser — the commuting cube-root solutions of
  Lemma A5 come in multiples, e.g. N = 9, α = (9), γ = (3,1⁶): n = 168, M₁ = 81, M₂ = 405, ratio 3.37).
So the second moment is systematically above Poisson and the excess grows with N; there is no evidence of an exact
Poisson or other simple law per block.

## C. Is #Sol a function of polynomially computable pair data? — FAILED, with a PROVED obstruction at N = 7

Trivial remark (PROVED). #Sol(a,c) is invariant under simultaneous conjugation, hence a function of the pair class; the
pair class is the isomorphism class of a 2-coloured functional graph of valence ≤ 4, so a canonical form *is* computable in
polynomial time (Luks, bounded valence). "Function of the pair class" therefore carries no complexity content; the
question is whether #Sol is an efficiently computable function of that canonical form. Number of pair classes with a
fixed = number of C(a)-orbits on S_N = (1/|C(a)|) Σ_{g∈C(a)} z_{ctype g} (Burnside): 136 (N=6, 6-cycle), 726 (7-cycle),
5100 (8-cycle), 40362 (9-cycle) — i.e. N!/|C(a)| up to lower order, exponential.

Invariant tested: J_L = (multiset of orbit sizes of ⟨a,c⟩) + (ctype w(a,c) for every cyclically reduced word w of
length ≤ L in a^{±1}, c^{±1}, one representative per rotation/inversion class, no proper powers: 2,2,4,9,24,58,156,405
words for L = 1..8). Computed for all c ∈ S_N at once (hashing two independent 61/59-bit rolling hashes; collisions can
only merge classes, i.e. only make the test *more* likely to fail, never produce a false "determined").
Tested: all α at N = 6, 7 with L ≤ 8; all α at N = 8 with L ≤ 7; α ∈ {(9),(8,1),(7,2),(5,4),(3,3,3)} at N = 9 with L ≤ 6.

**Findings (VERIFIED, `phaseC_word_invariants`).**
1. Whenever #Sol becomes a function of J_L, J_L is (or is within a few classes of) a *complete* invariant of the pair
   class: e.g. N = 6, a = 6-cycle: #Sol determined first at L = 7, where J_7 has 136 classes = number of pair classes;
   N = 7, a of type (5,2),(5,1,1),(4,3),(3,2,1,1): first at L = 7, with 536/536, 536/536, 448/448, 484/484 classes.
   The only cases where #Sol is determined with a strictly coarser invariant are (4,2),(4,1,1),(2,2,1,1) at N=6 and (3,2,2)
   at N=7 (104/108, 105/108 for solvability, 68/70, 252/254 classes), and the trivial α = (1^N), (2,1^{N−2}).
2. For a = N-cycle (and (N−1,1), (N−2,2), …) no tested L decides #Sol or even solvability:
   N = 7: J_8 has 705 classes vs 726 pair classes, not deciding; N = 8: J_7 has 5089 vs 5100; N = 9: J_6 has 40237 vs 40362.
3. **Theorem C.1 (PROVED by finite computation, N = 7).** Let a = (0 1 2 3 4 5 6), c₁ = [0,4,5,1,3,6,2], c₂ = [0,2,5,6,3,1,4]
   (both of cycle type (3,3,1)). Then #Sol(a,c₁) = 0 and #Sol(a,c₂) = 1 (brute force over S_7, pure Python, independent of
   the numpy pipeline), while the assignment w(a,c₁) ↦ w(a,c₂) is a well-defined isomorphism φ : ⟨a,c₁⟩ → ⟨a,c₂⟩ (both of
   order 168, ≅ PSL(2,7), transitive, distinct but conjugate subgroups of S_7) with ctype(φ(g)) = ctype(g) for all 168
   elements. Consequently, for **every** word w ∈ F₂ (not only |w| ≤ 8), ctype(w(a,c₁)) = ctype(w(a,c₂)); the orbit
   structures, |⟨a,c⟩|, the abstract isomorphism type of ⟨a,c⟩ and its permutation character (transported along φ) all
   coincide. No function of these data determines #Sol(a,c), nor even solvability.
   *Proof.* BFS over the Cayley graph of ⟨a,c₁⟩: for each reached g and each generator s ∈ {a,c₁} the image of s·g is
   compared with s'·φ(g); no conflict in 168 vertices ⇒ φ is a homomorphism of the free group factoring through ⟨a,c₁⟩;
   cycle types compared element-wise; the same BFS in the other direction is also conflict-free (so φ is bijective).
   Script: `word_iso_certificate` in `counting_track.py`; record `N7_witness_certificate` in the JSON. ∎
   Remark: the N = 6 witness of the earlier tracks (a = 6-cycle, (0 5)(1 3) vs (0 5)(2 4)) is separated by a word of
   length 7 (certificate fails at group order 242 reached), so N = 7 is the smallest N for which we have a word-proof
   obstruction; the two pairs realise the two 7-point actions of PSL(2,7) interchanged by its outer automorphism, which
   preserves cycle types in S_7 but is not induced by conjugation on the pair.
4. Further witnesses (VERIFIED, hash-based): N = 8, a = 8-cycle, c = [0,3,7,1,6,5,4,2] (0 solutions) vs [0,6,4,3,2,7,1,5]
   (2 solutions), J_7 equal; N = 9, a = 9-cycle, c = [0,6,5,4,3,2,1,8,7] (0) vs [0,2,1,8,7,6,5,4,3] (1), J_6 equal.

**Verdict on Goal 1.** The only exact counting formula we have is the first moment (block average, Theorem A.1). #Sol is
not a function of any word-cycle-type / orbit / group-theoretic invariant of the pair (Theorem C.1, N = 7), and the
complete invariant has exponentially many values. Whether #Sol is in FP remains OPEN; nothing here shows it is #P-hard
either (Theorem B of the project shows the restricted count #{x : x³ = 1, x a x⁻¹ = c} is NP-hard to decide, hence
#Sol restricted to the class 3^{N/3} is #P-hard; the unrestricted count is not covered by that reduction because the
commuting cube-root solutions always exist on the hard family).

## D. Goal 2 — the per-λ route P3^λ (λ = ctype x)

**The obstacle, stated exactly (PROVED).** Fix λ and a representative x_λ. By the conjugation form, x = g x_λ g⁻¹ solves
x a x² = c iff (g⁻¹ a g, g⁻¹ c g) = (a′, x_λ a′ x_λ²) with a′ := g⁻¹ a g.
Hence P3^λ(a,c) ⟺ the diagonal S_N-orbit of the pair (a,c) meets the graph Γ_λ = {(a′, F_λ(a′)) : a′ ∈ class α} of the
bijection F_λ : a′ ↦ x_λ a′ x_λ² of S_N. For a *given* a′ this is simultaneous conjugacy of two pairs (polynomial), but a′
is a free parameter ranging over |α| = N!/z_α elements, and Γ_λ is not a union of diagonal orbits (F_λ does not commute
with conjugation unless the conjugator centralises x_λ: F_λ(h a′ h⁻¹) = h F_λ(a′) h⁻¹ iff h ∈ C(x_λ)). The enumeration
over a′ modulo C(x_λ) costs N!/(z_α z_λ)·O(poly), so fixing λ only divides the N!-search by z_λ (≤ N for λ = (N)).
Equivalently, P3^λ is the orbit-intersection problem "orbit of (a,c) ∩ Γ_λ ≠ ∅" where both sets are exponential and Γ_λ
has no visible group structure.

**Special case (PROVED).** If every part of λ is coprime to 3, put e = lcm(λ) and k = 3⁻¹ mod e; then x = (x³)^k for
x ∈ λ, and with z = x³ (also of type λ, cubing being a bijection on the class) P3^λ ⟺ ∃ z ∈ class λ : c = z^k a z^{1−k},
i.e. c is conjugate to a z by the power z^k of z. This is a one-variable two-sided equation of the same nature as P3; it
gives no reduction. For λ = (N), 3 ∤ N: ∃ N-cycle z with c = z^k (a z) z^{−k}.

**Polynomial λ (PROVED earlier, cited).** parts ≤ 2 (x = c a⁻¹, decided by ctype(c a⁻¹), a word of length 2); parts in
{1,3} ⟺ conjugacy of a and c by an element of type λ (NP-complete for λ = 3^{N/3}, Theorem B); bounded support.

**Numerical test (VERIFIED).** For each λ and each a, P3^λ(a,c) was tested over all c against J_L (same invariant as §C):
`phaseD_per_lambda` records the smallest L such that P3^λ is a function of J_L (None = not decided up to the maximal L).
N = 7 (all a, L ≤ 8), N = 8 (all a, L ≤ 7), N = 9 (five a, L ≤ 6):
- parts ≤ 2: decided at L = 2 for every a (consistent with x = c a⁻¹).
- parts in {1,3} (x³ = 1): decided at small L for every a tested: λ = (3,1^{N−3}) at L ≤ 3 (N = 7, 8), L ≤ 4 (N = 9, five a);
  (3,3,1) at L ≤ 5 (N=7); (3,3,1,1) at L ≤ 5 (N=8); (3,3,1,1,1) and (3,3,3) at L ≤ 4 (N=9, five a). This refines the earlier
  "I5 fails for λ = (3,1⁶) at N = 9" (I5 used words ≤ 3): words of length 4 suffice at N = 9 on the a tested. By Theorem B
  this cannot persist for all N unless P = NP (the hard family needs N ≥ 15), so this is small-N behaviour.
- All other λ (exact lists in `phaseD_per_lambda.summary_lists`):
  * undecided (None) for a = N-cycle within the tested L: N = 7 (L ≤ 8): (7), (5,1,1), (4,2,1), (3,2,2);
    N = 8 (L ≤ 7): (7,1), (6,1,1), (5,1,1,1), (4,2,2), (4,1⁴); N = 9 (L ≤ 6): (9), (8,1), (7,2), (7,1,1), (6,2,1), (5,4),
    (5,3,1), (5,2,2), (5,2,1,1), (5,1⁴), (4,4,1), (4,3,2), (4,3,1,1), (4,2,2,1), (4,2,1³), (3,2,2,1,1).
  * undecided for *some* a (but decided for the N-cycle at the largest L, where J_L is nearly complete): N = 7: (5,2),
    (4,3), (3,2,1,1); N = 8: (8), (6,2), (5,3), (5,2,1), (4,4), (4,3,1), (4,2,1,1), (3,2,2,1), (3,2,1³); N = 9 (five a):
    (6,3), (6,1,1,1), (4,1⁵), (3,3,2,1), (3,2,2,2) — decided for a = 9-cycle at L = 6 but not for all five a.
  * decided for all a tested but only at L equal (or close) to the length where J_L becomes essentially complete, so no
    evidence of a genuine criterion: (6,1), (4,1,1,1) at N = 7 (L = 7); (3,3,2) at N = 8 (L = 7); (3,2,1⁴) at N = 9 (L = 6).
  Summary: among the λ with a part ≥ 4 or of the form (3,2,…), all but four — (6,1), (4,1,1,1) [N=7], (3,3,2) [N=8],
  (3,2,1⁴) [N=9, five a] — fail to be a function of the tested J_L for some (N ≤ 9, a); those four are decided for
  every a tested, but only at L = 7, 7, 7, 6 respectively. Decisions by words of length ≤ 4 occur for these λ only when
  a has small support or is an involution (`phaseD_per_lambda.decidability_le4_by_a`): a = id (L = 0, P3^λ ⟺ c is a
  cube of an element of type λ), a = transposition (all such λ decided at L ≤ 2 — consistent with Lemma A3, support 2),
  a = 3-cycle (3 of 9 such λ at N=7, 5 of 14 at N=8, at L = 2..4), a = (2,2,2,2) (6 of 14 at N=8), a = (2,2,1⁴) (1 of 14).
  For every a with a part ≥ 4, or with ≥ 2 parts ≥ 3, or of type (3,2,…), (2,2,2,1), (2,2,1,1,1) at N = 7, and for all
  five a at N = 9, no such λ is decided below L = 5 (minimum first_len_decided = 5), where J_L already separates
  ≥ 95.5 % (N=7), ≥ 98.2 % (N=8), ≥ 99.0 % (N=9) of the pair classes of the N-cycle.
**What survives:** only the already-known polynomial λ (parts ≤ 2; bounded support) and — at N ≤ 9 only — the order-3
types. No new polynomial criterion for any λ with a part ≥ 4 or of the form (3,2,…) (CONJECTURED: none exists; the
per-λ problem retains the orbit-intersection obstacle above). Caveat: at N ≤ 9 the tested word lengths are close to the
length at which J_L separates all pair classes, so 'decided at the largest L' carries little information.

## E. Does anything here change the complexity status of P3?  — No.
PROVED: Theorem A.1 + Lemma A.2 (first moment / block average in closed form), Proposition B.1 (second-moment
reformulation), the exact statement of the per-λ obstacle and the coprime-to-3 reformulation, Theorem C.1 (N = 7
obstruction to all word/group invariants).
VERIFIED only: second-moment numerics, all J_L tables, per-λ tables.
OPEN: #Sol ∈ FP?; a class-algebra formula for M₂; polynomiality of P3^λ for λ of unbounded support with a part ≥ 4.

## Deviations
Word length capped at L ≤ 8 (N ≤ 7), L ≤ 7 (N = 8), L ≤ 6 and five of the 30 cycle types of a (N = 9); no N = 10 tables
(feasible, not run). Second moment: numerics only. Hash-based grouping (two independent moduli) for the J_L tests;
the N = 7 witness was re-verified exactly in pure Python.
