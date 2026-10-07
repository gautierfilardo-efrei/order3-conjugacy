# Novelty report — results added in v6 of `p3_article.tex` (JSC submission)

Date: 2026-10-07. Sources queried: Crossref REST API (`api.crossref.org/works`, bibliographic queries + per-DOI resolution) and the arXiv API (`export.arxiv.org/api/query`). **OpenAlex was NOT queried** (see Deviations). Every DOI listed below was resolved through Crossref on 2026-10-07; the corresponding BibTeX entries are in `novelty_refs.bib` (21 entries, keys disjoint from `p3_refs.bib`). Entries already in `p3_refs.bib` (Goldmann–Russell 2002, Idziak–Kawałek–Krzaczkowski 2022, Weiß 2020, Lohrey–Rosowski–Zetzsche MFCS 2022, Luks, Seress, Sims, …) are referred to by their existing keys and not duplicated.

Labels: **KNOWN** = the statement itself is in the literature; **SPECIAL CASE KNOWN** = a special case or the underlying phenomenon is known, the statement as formulated is not; **NEW** = no literature found after the searches recorded in §7 (a negative search result is evidence of absence only to the extent of the queries run).

---

## 1. Theorem `thm:main` / `thm:conjp` — Conj_3 and Conj_p (p ≥ 3 prime) are NP-complete

**Verdict: NEW** (no hit on "conjugating element of prescribed order", "order of the conjugating permutation", "elements of order p in a coset of a centraliser"; queries b1–b3 returned nothing on the problem).

Closest literature, to be cited for positioning:

| ref | what it covers | what it does not cover |
|---|---|---|
| Mattes, Ushakov, Weiß, *Complexity of Spherical Equations in Finite Groups*, LATIN 2024, LNCS, pp. 383–397, doi:10.1007/978-3-031-52113-3_27 (arXiv:2308.12841) | Spherical equations z₁⁻¹g₁z₁⋯z_k⁻¹g_kz_k = 1 (quadratic, several variables) over S_n and A_n **with n part of the input** are NP-complete (their Thm 16, Cor 18); polynomial when the group is fixed or given by its multiplication table. | One-variable equations; non-quadratic words (x appears three times in x a x²); constraints of the form x^p = 1 on the conjugator. |
| Lohrey, Rosowski, Zetzsche, *Membership problems in finite groups*, J. Algebra 675 (2025) 23–58, doi:10.1016/j.jalgebra.2025.03.011 (journal version of the MFCS 2022 paper already cited as `lohrey2022mfcs`) | Subset-sum / knapsack / membership in products of three cyclic permutation groups NP-complete, n in the input. | Equations with a conjugation; order constraints. Suggest updating `lohrey2022mfcs` to the journal version or citing both. |
| Goldmann–Russell 2002 (`goldmann2002`), Idziak et al. (`idziak2022`), Weiß 2020 (`weiss2020`); Aichinger–Grünbacher, arXiv:2503.07285 (2025, S₄; no DOI, not in .bib) | Equation satisfiability over a **fixed** finite group. | The group S_N growing with the input; this is exactly the regime of the paper. |
| Larrauri–Živný, *Solving Promise Equations over Monoids and Groups*, ACM TOCL 26 (2024) 1–24, doi:10.1145/3698106 | Promise/approximation variants of equations over fixed groups. | Group as input. |
| Hulpke, *Conjugacy classes in finite permutation groups via homomorphic images*, Math. Comp. 69 (1999) 1633–1651, doi:10.1090/S0025-5718-99-01157-6 | Practical computation of conjugacy classes and class representatives in permutation groups (GAP). | The decision problem "is there a conjugator of order p" and its complexity. |

Positioning sentence (Section 4 / Introduction): *"Hardness results for equations over symmetric groups with N part of the input are recent: Mattes, Ushakov and Weiß [mattes2024spherical] prove NP-completeness of spherical (multi-variable quadratic) equations over S_N and A_N, and Lohrey, Rosowski and Zetzsche [lohrey2025membership] of membership in products of three cyclic subgroups. Theorem `thm:conjp` concerns a different regime: a single unknown, a non-quadratic word, and an order constraint x^p = 1 on the conjugator; to our knowledge the complexity of 'conjugacy by an element of prescribed order' has not been studied before."*

Note on the proof: the reduction from N3DM / NpDM to the blockwise-translation family is self-contained; no prior reduction of this shape was found.

## 2. Theorem `thm:conj2` — Conj_2 (conjugacy by an involution) is polynomial, O(N³)

**Verdict: NEW as an algorithmic statement; SPECIAL CASE KNOWN** (c = a⁻¹).

- The case c = a⁻¹ is the notion of a *strongly real* element (x a x⁻¹ = a⁻¹ with x² = 1). In S_N every element is strongly real (every permutation is inverted by an involution — classical; every irreducible character of S_N is real-valued and the Frobenius–Schur indicators are +1). The literature on strongly real elements concerns other families: Tiep–Zalesski, *Real conjugacy classes in algebraic groups and finite groups of Lie type*, J. Group Theory 8 (2005), doi:10.1515/jgth.2005.8.3.291; Rämö, *Strongly real elements of orthogonal groups in even characteristic*, J. Group Theory 14 (2011), doi:10.1515/jgt.2010.036. Neither treats a prescribed target c ≠ a⁻¹, nor the algorithmic question. (Guralnick–Montgomery, arXiv:math/0703681, on Frobenius–Schur indicators of Drinfeld doubles of Weyl groups, is related to involutions inverting elements but its DOI could not be resolved during this run — Crossref returned HTTP 429 — so it is not in the .bib.)
- No paper was found on "given a, c, decide whether an involution conjugates a to c" (queries a1–a6). The orbit-pairing criterion (plain vs twisted isomorphism classes of ⟨a,c⟩-orbits, parity condition on self-twisted classes) appears to be new.

Positioning sentence: *"For c = a⁻¹ the question is whether a is strongly real, which in S_N always holds [tiep2005real, §1]; for general c the only invariant available is the cycle type of ca⁻¹ when x² = 1 is imposed on x in x a x² = c, and the problem of deciding whether some involution conjugates a to c does not seem to have been considered. Theorem `thm:conj2` settles it in O(N³) time."* (Point the reader to the distinction: strongly real ⇔ Conj_2(a, a⁻¹).)

## 3. Corollary `cor:classification` — dichotomy over primes (p = 2 polynomial, p ≥ 3 NP-complete)

**Verdict: NEW** (follows from §1–§2; no dichotomy of this kind found). No positioning sentence needed beyond §1–§2.

## 4. Theorem `thm:psl` — no word invariant decides P3 (PSL(2,7) witness at N = 7)

**Verdict: phenomenon KNOWN (classical Gassmann triple), application NEW.**

- The statement constructs an isomorphism φ : ⟨a,c₁⟩ → ⟨a,c₂⟩ ≅ PSL(2,7) with φ(a) = a that preserves cycle types, i.e. the two degree-7 permutation representations have the same permutation character but are not conjugate in S₇. This is the textbook Gassmann triple (G = PSL(2,7) ≅ GL(3,2) acting on the 7 points and the 7 lines of the Fano plane; the point and line stabilisers are non-conjugate subgroups with the same permutation character). Perlis, *On the equation ζ_K(s) = ζ_{K'}(s)*, J. Number Theory 9 (1977) 342–360, doi:10.1016/0022-314X(77)90070-1, shows this is the smallest-degree non-trivial Gassmann triple (degree 7) and uses it for arithmetically equivalent number fields; Sunada, *Riemannian coverings and isospectral manifolds*, Ann. Math. 121 (1985), doi:10.2307/1971195, for isospectral manifolds. Kammeyer–Kionke, *Gassmann triples with special cycle types and applications*, Proc. Edinburgh Math. Soc. 67 (2024) 1115–1124, doi:10.1017/S0013091524000579, is the recent reference relating cycle types in coset actions to almost-conjugacy — directly the vocabulary of `thm:psl`. (Gassmann's original 1926 note in Math. Z. 25 could not be located through Crossref — no DOI record found — and is therefore omitted from the .bib; cite via Perlis.)
- INFERRED, to be confirmed by the author against the recorded isomorphism in `counting_track.py`: the pair (a, c₁) ↦ (a, c₂) should correspond to the point/line duality (inverse–transpose automorphism of GL(3,2)) composed with an inner automorphism fixing a. If confirmed, the theorem can be stated conceptually ("the two Gassmann-equivalent degree-7 actions of PSL(2,7) distinguish P3") instead of as a bare computational certificate.
- What is new: the use of a Gassmann pair to show that **no** function of the word cycle types, orbit structure, group order, isomorphism type or permutation character of ⟨a,c⟩ decides solvability of x a x² = c. No prior use of Gassmann equivalence as an obstruction to equation-solving invariants was found (queries d1–d6).

Positioning sentence: *"The pairs of Theorem `thm:psl` are the classical Gassmann triple of degree 7 — the two actions of PSL(2,7) ≅ GL(3,2) on points and lines of the Fano plane, which share their permutation character without being conjugate in S₇ [perlis1977, sunada1985, kammeyer2024gassmann]. Gassmann equivalence is usually invoked to produce arithmetically equivalent fields or isospectral manifolds; here it shows that solvability of x a x² = c is not a function of any word-cycle-type invariant of the pair (a,c)."*

## 5. Theorem `thm:firstmoment` / Proposition `prop:moment` — first moment of the number of solutions in closed form

**Verdict: NEW as a statement; ingredients KNOWN.**

- Counting x with ctype(a x³) = γ combines two classical objects: (i) the number of cube roots of a permutation of given cycle type, i.e. the fibre sizes of x ↦ x³ — Chowla–Herstein–Moore, Canad. J. Math. 3 (1951), doi:10.4153/CJM-1951-038-3 (recursions for #{x : x^d = 1}); Moser, Canad. J. Math. 7 (1955), doi:10.4153/CJM-1955-021-8 (asymptotics of #{x^d = 1}); Pavlov, Math. USSR-Sb. 40 (1981) 349–362, doi:10.1070/SM1981v040n03ABEH001824 and 45 (1983) 243–255, doi:10.1070/SM1983v045n02ABEH002597 (number of solutions of x^k = a for given a, and its limit distribution); Chernoff, Discrete Math. 125 (1994) 123–127, doi:10.1016/0012-365X(94)90152-X (which permutations have p^l-th roots — the criterion on cycle multiplicities used implicitly in the closed form); Ishihara–Ochiai–Takegahara, Ann. Comb. 5 (2001) 197–210, doi:10.1007/PL00001300 (p-divisibility of #{x^p = 1}). (ii) Products of conjugacy classes / class-algebra structure constants in S_N (Frobenius' character formula): Goupil, Discrete Math. 79 (1990) 49–57, doi:10.1016/0012-365X(90)90054-L; Bédard–Goupil, Canad. Math. Bull. 35 (1992) 152–160, doi:10.4153/CMB-1992-022-9.
- No reference was found stating the identity for E[#solutions of x a x² = c] over c in a class, nor the closed form in terms of the cycle type of a (queries e1–e6). The paper should present the theorem as a routine but apparently unrecorded consequence of (i)+(ii), not as a deep result.

Positioning sentence: *"The first moment reduces to two classical counts — the fibres of the cube map on S_N [chowla1951, moser1955, pavlov1981, chernoff1994] and the structure constants of the class algebra [goupil1990, bedard1992] — and the closed form of Theorem `thm:firstmoment` follows; we have not found it stated in the literature."*

## 6. Lemma `lem:twoeq` (Conj_3 = system {x a x² = c, x a² x² = c²}), Proposition `prop:dual`, Proposition `prop:three`, Theorem `thm:transposition`

**Verdict: NEW (elementary).** No literature on positive one-variable words of x-length three over S_N (queries c1–c9 return only fixed-group equation complexity, free-group one-variable equations — Levin 1964, Chiswell–Remeslennikov 2000, Bormotov–Gilman–Myasnikov 2009 — and spherical equations). No citation needed; the general context is covered by the sentences of §1.

## 7. Journal positioning (JSC)

Recent JSC papers on permutation-group algorithms, all Crossref-verified (search restricted to container-title "Journal of Symbolic Computation", 2019–2025):

- Jefferson, Pfeiffer, Waldecker, *New refiners for permutation group search*, JSC 92 (2019) 70–92, doi:10.1016/j.jsc.2017.12.003 — backtrack search (normalisers, intersections, conjugacy of groups).
- Jefferson, Waldecker, Wilson, *Perfect refiners for permutation group backtracking algorithms*, JSC 114 (2023) 18–36, doi:10.1016/j.jsc.2022.04.007.
- Chang, Jefferson, *Disjoint direct product decompositions of permutation groups*, JSC 108 (2022) 1–16, doi:10.1016/j.jsc.2021.04.003.
- Požar, *Fast computation of the centralizer of a permutation group in the symmetric group*, JSC 123 (2024) 102287, doi:10.1016/j.jsc.2023.102287 — centralisers in S_N, the object the coset of conjugators is built from.
- Cannon, Holt, Unger, *The use of permutation representations in structural computations in large finite matrix groups*, JSC 92 (2019), doi:10.1016/j.jsc.2018.09.001 (verified in the search listing; not added to the .bib — not needed).

Positioning sentence for the cover letter / introduction: *"Algorithms for permutation groups — backtrack refiners for conjugacy and normaliser problems [jefferson2019jsc, jefferson2023jsc], centralisers in S_N [pozar2024jsc], product decompositions [chang2022jsc] — are a recurrent topic of this journal; the present paper studies the complexity of a one-variable equation over S_N whose solution set is a union of cosets of such centralisers, and shows that an order constraint on the conjugator moves the problem from polynomial to NP-complete."* The JSC scope statement itself (Elsevier page) was not fetched during this run — the journal website is outside the sandbox allowlist; the positioning above relies on the published JSC papers listed.

## 8. Query log (what was searched)

Crossref `query.bibliographic` and arXiv `search_query` for: strongly real elements / conjugate by an involution / inverted by an involution / real classes of symmetric groups (a1–a6); conjugating element of prescribed order, coset of a centraliser, conjugacy problem complexity in permutation groups (b1–b3); one-variable equations over symmetric groups, quadratic/spherical equations, equations over finite groups with group as input, permutation groups + NP-complete (c1–c9); Gassmann equivalence, almost conjugate subgroups, cycle-type-preserving isomorphism, Brauer pairs, Perlis, Sunada (d1–d6); number of solutions of x^k = a, cube roots, products of conjugacy classes, Frobenius formula, connection coefficients (e1–e6); JSC permutation-group algorithms since 2019 (f). Raw hit lists: `hits_all.json`, `hits_followup.json`; Crossref resolution records: `verified.json`.

## 9. Deviations

- **OpenAlex not queried.** The credential approval card (`host.credentials.request("openalex")` and a `credentials=['OpenAlex']` python cell) stayed pending for the whole run without a user response, so no OpenAlex request was issued. Crossref + arXiv coverage is good for the classical and the TCS references but OpenAlex's full-text/abstract search could still surface a paper on "conjugacy by an element of prescribed order" that bibliographic-title search misses; the NEW verdicts in §1–§3 and §5 are therefore **provisional**.
- Gassmann (1926) and Guralnick–Montgomery (2009) were identified but their DOIs could not be verified (no Crossref record / HTTP 429); they are mentioned in the text and excluded from the .bib.
- The JSC scope statement was not fetched (journal website not reachable from the sandbox).
