# Referee audit of `p3_article.tex` v6 (+ `theorems_new.tex`) — Journal of Symbolic Computation

Scope: every lemma/proposition/theorem/corollary of Sections 2–4 and 6 read step by step; Section 5 read for internal consistency;
abstract / results list / open problems / numbering / notation / numerical claims cross-checked. Every doubtful point was tested by
brute force (`audit_checks.py`, results in `audit_checks_results.json`; stdlib + numpy, N ≤ 9). The manuscript itself was not edited;
replacement LaTeX is in `audit_patches.tex`, each snippet preceded by the exact original text it replaces.

Verdict: **the mathematics of the hardness part (Thm 4.2 main, Lemma 4.3 block, Lemma 4.4 affine, Lemma affinep, Thm conjp, Thm conj2,
Lemma twoeq) is correct**; the structural part is correct but two statements are false as written (Cor. lambda, Prop. moment's
complexity claim), one lemma is only sketched (Lemma support), one proof has a confused case split (Thm transposition), and several
numbers are not traceable to the shipped record files. Three BLOCKING items, fifteen MUST FIX, twelve STYLE.

Labels used below: PROVED = I re-derived the proof; VERIFIED = enumeration (what, how many, script); FAILED = claim not reproduced.

---

## A. BLOCKING

### B1. Corollary `cor:lambda`, third clause, contradicts Open Problem 2
*Location.* Sec. 4.4, `cor:lambda`: "and in time $N^{O(s)}$ when $\lambda$ has at most $s$ parts larger than $1$."
*Problem.* $\lambda=(N)$ has one part larger than 1, so the clause asserts $\Pthree^{(N)}\in\mathrm P$ (time $N^{O(1)}$); Open Problem 2
asks whether $\Pthree^{(N)}$ is NP-complete and says no polynomial criterion was found. Nothing in the paper proves the clause (Lemma
`lem:support` bounds $|\supp a|$, not the shape of $\lambda$). What is true and trivial: if the parts of $\lambda$ larger than 1 sum to
at most $s$ then $|\supp x|\le s$ and the $\le N^{s}$ candidates are checked directly.
*Fix.* Patch P-B1 (restate the clause with $|\supp x|\le s$).

### B2. Proposition `prop:moment`: "computable in polynomial time from the class algebra" is unjustified (and very likely false)
*Location.* Sec. 3.2, statement of `prop:moment`, last sentence; proof: "both polynomial-time computable".
*Problem.* The closed form (Thm `thm:firstmoment`) is a sum over all $p(N)=e^{\Theta(\sqrt N)}$ classes $\nu$ of a connection
coefficient $K(\gamma,\nu;\alpha)$, itself a sum over all $p(N)$ irreducible characters. No polynomial-time algorithm for connection
coefficients of $\Sym(N)$ (or for individual character values) is known; the paper gives none. In a complexity paper a false "polynomial
time" is a referee stopper. The identity itself is PROVED and VERIFIED (closed form against brute force: all $(\alpha,\gamma)$ blocks,
$N\le6$, 0 discrepancies; $\rho_3$ formula against cube-root counts, all classes $N\le7$; `audit_checks.py`).
*Fix.* Patch P-B2: drop the complexity claim; state "evaluation polynomial in $p(N)$ given the character table"; and merge `prop:moment`
into `thm:firstmoment` (duplicate statement, see M6).

### B3. Lemma `lem:support` is a "proof sketch" with an incorrect step
*Location.* Sec. 3.1, `lem:support`, "Proof sketch".
*Problem.* (i) The results list (item 2) and the abstract sell bounded support as a *proved* polynomial case; a sketch is not acceptable.
(ii) The sentence "The prescribed images involve at most $2s$ cycles of $t$ and fix their grouping and phases up to a number of choices
bounded by a function of $s$" is wrong: a prescribed image $x(p)=q$ with $p,q$ in different cycles of $t$ forces those two cycles into a
triple whose *third* cycle is a free choice among up to $N$ cycles; the number of choices is $N^{O(s)}$, not $f(s)$. The conclusion
$N^{O(s)}$ survives, the argument does not. (iii) The verification note "(Verified exhaustively for $s=2$, $N=6$: 30/30)" sits inside the
proof environment.
*Fix.* Patch P-B3 gives a complete proof (cube roots of $t$ commute with $t$, hence are length-preserving permutations $\pi$ of the
cycles of $t$ with $\pi$-cycles of length 1 or 3 plus phases; $\le N^{s}$ injections $f$, $\le N^{2s}$ completions on the cycles meeting
$\supp a$, divisibility-by-3 condition on the rest). The algorithm of that proof was implemented and VERIFIED against brute force on
360 instances ($N=5,6,7$; $s\in\{0,2,3,4\}$; 0 discrepancies; `support_decide` in `audit_checks.py`).

---

## B. MUST FIX

### M1. Theorem `thm:transposition`: case split and "one solution per admissible $\ell$"
*Location.* Sec. 3.1, proof of `thm:transposition`.
*Problems.* (a) "The points $0$ and $d$ lie in different arcs" is asserted without reason. Reason: $y$ commutes with $y^{3}=(u\,v)a$ and
has no cycle of length divisible by 3, so $y$ preserves each arc; as $y(0)=u\in C_1$ and $y(d)=v\in C_2$, $0\in C_1$, $d\in C_2$. This
also disposes of the case $\ell=N-\ell$ (nothing special happens: $y$ cannot swap the arcs). (b) "It remains to impose
$\{y(0),y(d)\}=\{u,v\}$" is a *set* condition although $u,v$ were *defined* as $y(0),y(d)$; the ensuing "case $0\in C_2$, $d\in C_1$ …
gives the same $y$" cannot occur under the definitions and confuses the count (a reader may think $\ell$ and $N-\ell$ are both counted
for one $y$). (c) "One solution per admissible $\ell$" needs injectivity of $y\mapsto\ell$: $\ell$ is the length of the $y$-cycle through
$0$, hence determined by $y$. (d) The convention "position" on an arc and the identity $\varepsilon(\ell)\equiv-e_1\pmod\ell$ should be
written. Statement and formula are correct: VERIFIED by brute force for all $d$, $N=2,\dots,8$ (28/28; the paper's record claims
$N\le11$). The admissible sets for $\ell$ and $N-\ell$ were never both non-empty for $N\le12$.
*Fix.* Patch P-M1 (rewritten proof).

### M2. Lemma `lem:affinep`: parenthetical false for composite $L$; non-degeneracy of the translation condition unstated
*Location.* Sec. 4.3, proof of `lem:affinep`.
*Problems.* The lemma is stated for arbitrary $L\ge3$. The parenthetical "(and, when $m_i\ne1$, $t_i$ is then arbitrary; when
$m_i=1$, $t_i=0$)" is false for composite $L$: for $L=9$, $p=3$, $m_i=4$ one has $m_i^{3}\equiv1$ but $x^{3}(z)=z+3t_i$, so only
$t_i\in\{0,3,6\}$ work (checked). Also for $m_i=1$ the condition is $pt_i\equiv0\pmod L$, so $t_i=0$ is *a* solution, unique only when
$p\nmid L$. The lemma's *iff* is unaffected because $t_i=0$ always works. On $p$-cycles, "one linear condition on $p$ free parameters"
should say that $t_{i_p}$ has coefficient $1$ (the same remark is made correctly in `lem:affine`), so the condition is non-degenerate and
has exactly $L^{p-1}$ solutions. VERIFIED: criterion against enumeration of all $L^{k}k!$ conjugators for
$(L,k,p)\in\{(9,3,3),(8,3,3),(7,3,3),(9,2,3),(5,3,3),(4,3,2),(5,2,2),(8,2,2)\}$, 40 random slope vectors each, 0 discrepancies; the
"type $p^{N/p}$" clause held on every positive instance satisfying its hypothesis.
*Fix.* Patch P-M2.

### M3. Choice of the prime $L$ (Thm `thm:main`) and the prime alternative in Thm `thm:conjp`
*Location.* Proof of `thm:main`, first paragraph; statement of `thm:conjp` "(or any prime $L$ with $p\nmid L-1$)" and the proof's
"alternatively the least prime $L$ with $L-1\ge T_0$ and $p\nmid L-1$".
*Problems.* (a) "for $B\le10^{4}$ one checks directly that $L\le2(12B+2)$" adds nothing to the asymptotic statement and suggests the
polynomial bound is only verified numerically. The correct argument: the least prime $\equiv2\pmod3$ above $12B+1$ is $O(B)$ by the
prime number theorem for arithmetic progressions; the implied constant exists and need not be known since $L$ is found by trial division
over successive candidates, each test polynomial in $B$. (b) In `thm:conjp` the prime alternative is asserted with no bound on $L$; it needs
the same remark (primes $L\not\equiv1\pmod p$ have Dirichlet density $(p-2)/(p-1)>0$, so the least one above $T_0$ is $O(T_0)$), plus
$L\ne p$ (automatic since $L>T_0>p$). Either add the sentence or drop the parenthetical; the $2$-power variant needs none of this.
*Fix.* Patch P-M3a, P-M3b.

### M4. Proof of Theorem `thm:main`: wrong interval in the case $\gamma=1,\beta=2$
*Location.* "and $\gamma=1,\beta=2$ gives $S+M-B\in[3B+3,5B]$". With $S\ge3$, $M=3B$: $S+M-B\ge2B+3$, so the interval is $[2B+3,5B]$.
The conclusion (never $\equiv0$) is unaffected. The whole ten-case analysis was VERIFIED by enumeration for $B\le12$, $T\in\{12B+1,12B+2,12B+3\}$,
all $(\alpha,\beta,\gamma)$ and all $S\in[3,3B]$: 0 discrepancies. The analogous analysis of `thm:conjp` ("admissible $p$-sets are exactly
the matching tuples") was VERIFIED arithmetically for $p=3$ ($B'\le8$), $p=5$ ($B'\le4$), $p=7$ ($B'=2$): 0 discrepancies.
*Fix.* Patch P-M4.

### M5. Section 4.3 preamble: tautology and duplicated setup
*Location.* Lines "$\mathrm{Conj}_3$ is the problem $\Dzero$ of Section~\ref{sec:hard}" and the paragraph "The blockwise family and
the affine criterion for $x^{p}=1$ … Fix $L\ge3$ … for $\sigma\in\Sym(k)$ and $t\in\ZZ_L^{k}$."
*Problems.* `\Dzero` is defined as `\mathsf{Conj}_3`, so the sentence prints "Conj$_3$ is the problem Conj$_3$ of Section 4" — a
tautology, and "Section 4" is the section the reader is in. The blockwise family is defined a second time (already in Sec. 4.1); the
notation changes from $\ZZ_L^{\times}$ (4.1) to $\ZZ_L^{*}$ (4.3). Also, as the section now contains all primes, its title "Conjugacy by an
element of order three is NP-complete" is too narrow.
*Fix.* Patch P-M5 (labels `def:dzero`, `sec:block`; one-paragraph preamble; unified $\ZZ_L^{\times}$; new section title).

### M6. Duplicated statements
*Location.* (a) `prop:moment` and `thm:firstmoment` both state $\sum_{c\in\gamma}\#\mathrm{Sol}(a,c)=\#\{x:\ctype(ax^{3})=\gamma\}$ with
separate proofs. (b) `prop:three` ends with "$\Pthree$ is self-dual: $xax^{2}=c$ iff $ycy^{2}=a$ with $y=x^{-1}$" and `prop:dual` restates it.
(c) Lemma `lem:block` is restated in prose in 4.3 (see M5).
*Fix.* P-B2 merges (a); P-M6 removes the duplicate from `prop:three` and lets `prop:dual` carry the duality (it is needed there for the
$N^{O(s)}$ consequence and for `thm:transposition`).

### M7. Introduction: "Length one is trivial or conjugacy"
*Location.* Sec. 1, paragraph "For one-variable words of small $x$-length…".
*Problem.* A word of $x$-length one is $c_0x^{\pm1}c_1=c$, never a conjugacy (conjugacy has $x$-length two), contradicting Prop. `prop:short`.
*Fix.* Patch P-M7.

### M8. Numerical claims not traceable to a shipped record, or not reproduced
(a) **FAILED** — Remark `rem:notproved`: "in an exhaustive census for $N\le8$ they [rigid pairs] are about $0.8\%$ of the solvable pairs,
mostly with a unique solution." Recomputed (`rigid_census2`, `audit_checks.py`; rigid = solvable and every solution has order exactly 3):
$N=5,6,7,8$ give $7.1\%, 6.2\%, 3.3\%, 1.1\%$ of solvable pairs (with "order dividing 3": $7.9, 6.3, 3.3, 1.1\%$), and the fraction of
rigid pairs with a unique solution is $42, 59, 61, 78\%$. Neither convention yields $0.8\%$; the figure must be recomputed or its exact
definition stated (the decreasing trend is the real message). Proposed text in P-M8a is labelled provisional until the author's record agrees.
(b) The cost slopes "1.12, 1.44 and 1.21 bits per point of $N$" (after `lem:inv`) are in none of the five JSON files of `data/` (grep on the
values); they need a record file or must be dropped.
(c) **FAILED** — Sec. 5, "from $m=70$ with the first-free-point fit to $m=89$–$116$": with the fits of Table 1 and the rule
$\TR(m)\ge2^{128}$, the first-free-point fits give $m=74$ (linear, $1.85m-7.53\ge128$) and $m=59$ ($m\log m$: `wreath_summary.json` itself lists
$128.17$ bits at $m=59$); $89$ and $116$ are reproduced for the most-constrained fits. "$m=70$" matches neither; state the rule and recompute.
(d) "all $1813$ pairs of classes $(\alpha,\gamma)$ for $N\le9$": $1813=\sum_{N=3}^{9}p(N)^{2}$, so write "$3\le N\le9$".
(e) Theorems `thm:conj2` and `thm:conjp` carry **no** verification record although the campaign memo has them (Conj$_2$: 1520/1520 pairs,
$N\le7$, brute force over all involutions; affine criterion for $p\in\{3,5\}$: 74/74 instances, 10 positive; the "type $p^{N/p}$" clause
tested on *one* positive instance only). Add them with their honest scope (P-M8e), matching the exact numbers in `conj_p_results.json`.
My independent checks: `thm:conj2` criterion against brute force on all 873 pairs ($N\le6$, one $a$ per class, every $c\sim a$) and 300
random pairs at $N=7$: 0 discrepancies.
(f) The statistics paragraph: "the probabilities of exactly one and exactly two solutions are $0.31$–$0.34$ and $0.16$–$0.19$" holds for
$N=6,\dots,9$ only (reproduced: $P_1=0.339,0.311,0.338,0.331$; $P_2=0.179,0.185,0.159,0.163$); at $N=5$, $P_1=0.49$. Say "for $N=6,\dots,9$".
The solvable fractions $0.583,0.708,0.610,0.587,0.591,0.591$ are reproduced exactly. Table 1 entries are reproduced from `wreath_summary.json`.

### M9. Open Problem 4: "$\#\mathrm P$-hard by the reduction of Theorem main"
*Problem.* The reduction is not stated to be parsimonious and no reference for the hardness of counting N3DM solutions is given. What the
reduction gives: every admissible triple contributes exactly $L^{2}$ translation vectors, so $\#\{x:x^{3}=1,xax^{-1}=c\}=L^{2q}\cdot\#\{\text{N3DM
solutions}\}$; hence the count is at least as hard as counting N3DM solutions. Say that, or supply a verified reference for $\#$N3DM.
*Fix.* Patch P-M9.

### M10. Lemma `lem:inv2`(a): "(a) is Lemma inv1" and the $O(|O|)$ decision procedure
*Problem.* Lemma `lem:inv1` is about a global involution $x$; (a) is about an abstract twisted isomorphism $O\to Gv$ and must be argued
(same computation: $\varphi(gj_0)=g^{\tau}v$). The procedure must check *all* $2|O|$ edges of the Schreier graph (not a spanning tree), that
the map is injective, and that its image is $Gv$ (automatic: a non-empty $a,c$-invariant subset of an orbit). Also note that $\varphi a=c\varphi$
on a bijection implies $\varphi a^{-1}=c^{-1}\varphi$, so inverse edges need no separate treatment.
*Fix.* Patch P-M10. (The complexity bound $O(N^{3})$ in `thm:conj2` is correct but loose: $\sum_{O,O'}|O||O'|=N^{2}$, so the all-pairs
tests cost $O(N^{2})$; optional.)

### M11. Theorem `thm:psl`: "well-defined bijective homomorphism" and "same permutation character"
*Problem.* $\varphi$ is defined on *elements* via one BFS word per element, so it is a function by construction; what the $168\times168$
check establishes is the homomorphism property, which implies independence of the chosen words. "Same permutation character" should be
stated as $\pi_2\circ\varphi=\pi_1$ (fixed-point counts are read off cycle types). VERIFIED: $|G_i|=168$, both transitive, $\varphi$ a
cycle-type-preserving isomorphism, solution counts $0$ and $1$ (`audit_checks.py`).
*Fix.* Patch P-M11.

### M12. Lemma `lem:cent`: sufficiency asserted, verification degenerate
*Problem.* "conversely such a grouping yields a cube root" is asserted; for singletons it needs the invertibility of $I+P+P^{2}$ over $\ZZ_L$
(determinant $\pm3$ when $3\nmid\ell$), for triples the identification of the image of $I+P+P^{2}$ with the vectors of equal rotation sums.
The shipped verification covers only $L=2$, where $\ZZ_2\wr\Sym(n)$ is degenerate. VERIFIED here for $(L,n)\in\{(2,4),(4,3),(5,3),(7,3),(4,4),(5,4),(2,6)\}$
(every element of $\Cent(a)$, $70\,800$ elements in total; 0 discrepancies; `cent_criterion`). The phrase "three ways of reading the rotation sum
… using $\gcd(L,3)=1$ to divide by $3$" is imprecise: each of the three $\rho$-cycle sums equals the total rotation sum along the $\pi$-cycle,
no division needed; $\gcd(L,3)=1$ is needed for the singletons.
*Fix.* Patch P-M12.

### M13. Observation `obs:cost`: "exponent about one fifth of brute force"
*Problem.* Brute force over $\Alt(m)$ costs $\log_2|\Alt(m)|\approx m\log_2m-1.44m$ bits; the fitted coefficients are $0.374$ (first free
point) and $0.223$ (most constrained): between a fifth and two fifths, not "about one fifth".
*Fix.* Patch P-M13.

### M14. Statement of Theorem `thm:main`: dependence of $L$ on the instance
"Hardness holds already when $a$ consists of $k$ disjoint cycles of one prime length $L\equiv2\pmod3$" — $L$ grows with the instance; write
"for a prime $L\equiv2\pmod3$ depending on the instance" (same in `thm:conjp`). Patch P-M14.

### M15. `theorems_new.tex` (if used for further integration)
Its definition of $\rho_3$ ("$\#\{x:\ctype(x^{3})=\nu$ for a fixed $y\in\nu$, $x^{3}=y\}$") is garbled (the article version is correct) and it
references a non-existent label `sec:hardness`. Not compiled, so harmless, but do not copy from it.

---

## C. STYLE

S1. Figures `fig:red`, `fig:stat`, `fig:wreath` are never referenced in the text (only `fig:cost`, `fig:var` are). Add "(Figure~\ref{…})".
S2. Open Problem 1 cites "Remark after Lemma~\ref{lem:twoeq}": remarks are numbered; label it (`rem:padding`) and cite the number.
S3. "The tautological complete invariant … has $N!/|\Cent(a)|$ values in general": the number of $\Cent(a)$-orbits on $\Sym(N)$ is at least
$N!/|\Cent(a)|$ (the paper's own $726>720$ for a $7$-cycle shows the inequality). Write "at least".
S4. Abstract: "so it is the shortest one-variable shape … whose complexity is not classified" — mixed-sign words of $x$-length three are also
unclassified (Remark after `prop:three`); write "a shortest".
S5. Verification notes inside proof environments (`lem:support`, `lem:cent`, `prop:moment`, `thm:psl`): move them after the proof.
S6. Section 4 title (see M5). Subsection 4.4 title "Consequences for $\Pthree$" fine.
S7. Notation: $\ZZ_L^{\times}$ vs $\ZZ_L^{*}$; `\mathrm{Conj}_p` (roman) vs `\Dzero` = `\mathsf{Conj}_3` (sans-serif) — define
`\Conj{p}` once as `\mathsf{Conj}_{p}` and set `\Dzero` to `\Conj{3}`.
S8. "found no other identical solution" (after `lem:inv`) → "no other word in $a^{\pm1},c^{\pm1}$ that is identically a solution".
S9. `thm:conj2`: define "involutive twisted automorphism" (a twisted automorphism $x$ of $O$ with $x^{2}=1$) in the statement.
S10. Reproducibility: the release tag `v0.3.0` must be bumped if new scripts (`verify_round2.py`, `audit_checks.py`) are added; the files named
in the proof of `thm:psl` (`counting_track.py`, `verify_campaign.py`) must be in the release.
S11. `prop:short`: the citation for the square-root criterion is `blum1974`; the criterion is classical — fine, but the phrase "cycles of even
length occur an even number of times for every length" should read "for every even length $m$, the number of $m$-cycles is even".
S12. English: "one checks directly" (M3), "the mechanisms of Theorem main are polynomial here" → "the mechanisms of Theorem~\ref{thm:main} are
polynomial-time here"; "Finite certificate:" → "Proof by finite certificate."

---

## D. Statement-by-statement verdict (for the v7 checklist)

| Statement | Proof status | Numerical check (this audit) | Action |
|---|---|---|---|
| Prop `prop:short` | PROVED (correct) | — | S11 wording |
| Prop `prop:three` | PROVED (bijection re-derived) | — | M6 (remove duplicate duality) |
| Lemma `lem:conj` | PROVED | — | none |
| Lemma `lem:support` | sketch, one wrong step | algorithm of new proof: 360/360 | **B3: replace by full proof** |
| Lemma `lem:cent` | sufficiency asserted | 70 800 elements, $L\in\{2,4,5,7\}$: 0 disc. | M12: add argument |
| Lemma `lem:inv` | PROVED | — | S8 |
| Prop `prop:dual` | PROVED | — | M6 |
| Thm `thm:transposition` | correct, case split confused | 28/28 ($N\le8$, all $d$) | M1: rewrite proof |
| Prop `prop:moment` | identity PROVED; complexity claim unjustified | closed form $N\le6$: 0 disc. | **B2: fix + merge** |
| Thm `thm:firstmoment` | PROVED ($\rho_3$ and Frobenius step re-derived) | $\rho_3$ all classes $N\le7$: 0 disc. | fold `prop:moment` into it |
| Thm `thm:psl` | PROVED (finite certificate re-run) | 168, hom, bijective, ctype-preserving, counts 0/1 | M11 wording |
| Lemma `lem:block` | PROVED (any $L\ge2$) | implicit in affine checks | none |
| Lemma `lem:affine` | PROVED | 320 instances incl. composite $L$ | none (typo-free) |
| Thm `thm:main` | PROVED; one interval typo; prime bound wording | ten cases $B\le12$: 0 disc. | M3a, M4, M14 |
| Lemma `lem:affinep` | PROVED (iff); parenthetical false for composite $L$ | 320 instances: 0 disc. | M2 |
| Thm `thm:conjp` | PROVED (a)(b)(c) re-derived | admissible sets $p=3,5,7$: 0 disc. | M3b, M8e, M14 |
| Lemma `lem:inv1` | PROVED | — | none |
| Lemma `lem:inv2` | (a) under-argued | — | M10 |
| Thm `thm:conj2` | PROVED (necessity, sufficiency, class well-definedness checked) | 873 exhaustive + 300 random: 0 disc. | M8e, S9 |
| Cor `cor:classification` | follows | — | none |
| Cor `cor:lambda` | third clause FALSE as written | — | **B1** |
| Cor `cor:noenum` | follows from `cor:lambda`(i) | — | none |
| Rem `rem:notproved` | commuting solution PROVED; census FAILED | 7.1/6.2/3.3/1.1 % | M8a |
| Lemma `lem:twoeq` | PROVED | 400 random pairs $N=5$: 0 disc. | none |
| Thm `thm:wreath` | PROVED (computation re-done) | — | none |
| Obs `obs:cost`, Table 1 | empirical; Table reproduced from `wreath_summary.json` | $m=70$ FAILED | M8c, M13 |
| Open problems 1–7 | consistent after B1 fix; OP4 over-claims | — | M9 |

Remaining for a human reader (not checkable by enumeration): the asymptotic prime bounds (M3; standard but must be worded as in the
patch), the $\#$N3DM reference (M9), the exact scopes of the campaign verification records to be quoted (M8e), and the recomputation of the
rigid census and of the sizing degree (M8a, M8c).
