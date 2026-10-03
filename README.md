# Order-3 conjugacy and one-variable equations of length three over symmetric groups

Code, exhaustive tables and verification scripts accompanying the manuscript

> *One-variable equations of length three over symmetric groups: normal forms, the hardness of order-three conjugacy, and layered permutation walks* (G.-E. Filardo, 2026; manuscript prepared for the *Journal of Symbolic Computation*, `paper/`).

The manuscript studies the equation **x a x² = c** with unknown `x` and constants `a, c` in `S_N`, `N` part of the input (problem **P3**), and proves that its restriction to solutions of order three — *are two permutations conjugate by an element of order 3?* — is **NP-complete** (Section 4 of the paper), by reduction from Numerical 3-Dimensional Matching. The complexity of the unrestricted equation is open.

Conventions throughout: permutations are image tuples on `{0..N-1}`, products are read right to left, `(pq)(j) = p(q(j))`; the word `x a x x` applied to `j` is `x(a(x(x(j))))`. Python ≥ 3.10, standard library only unless stated (numpy in `structure/algo_track.py`; z3-solver optional in `solver/round_solver_study.py`). Seeds are derived by SHA-256 from labelled strings, so every non-timing value is reproducible; timings are local.

## Layout

| Directory | Contents | Paper |
|---|---|---|
| `hardness/` | `order3_track.py`: block lemma, affine criterion, N3DM → Conj₃ reduction, rigidity census (`order3_results.json`, `order3_notes.md`). `verify_lemmaC.py`: **independent** check of the affine criterion against exhaustive enumeration of all `L^k k!` conjugators (160/160). `verify_reduction.py`: exhaustive case analysis of the reduction for `B ≤ 40`, prime bound `L ≤ 2(12B+2)` for `B ≤ 10⁴`, 180 random N3DM instances vs brute force, exhaustive conjugator checks at `q = 1`. | §4 |
| `structure/` | `round_algebra.py`: wreath normal form (exhaustive at `m = 5`), conjugation normal form, special cases, solution statistics (`round_algebra_results.json`, notes). `algo_track.py`: exhaustive tables of P3 for `N = 3..9`, invariant tests and the `N = 6` witness, first-moment identity check, special families (`algo_results*.json`, `algo_notes.md`). | §2, §3, §5.1 |
| `solver/` | `word_solver.py`: generic chain-propagation solver for one-variable words (does **not** enforce the wreath block structure by itself). `round_solver_study.py`: instrumented solver for `x a x² = c` in `A_m wr C_2`, branching heuristics, cycle-type relaxation, conjugated form, Z3 encoding, profile (`solver_variants.json/csv`, `round_solver_profile.json`, `profile.md`). `onset_model.py`: random-coverage model of the search tree. Scaling records `xax2_scaling.json`, `cost_models.json`, `constant_order_scaling.json`. | §5.2 |
| `failed_gadgets/` | `gadget_b_track.py`: the attempted 3-partition gadget through the choice of the conjugate `b`; **fails** (unintended block splits on 77/99 instances). Kept as a record for anyone attempting the hardness of P3 itself. | §4, closing remark |
| `literature/` | `voie_A_refs.bib` (87 entries, 83 DOIs verified) and the survey `voie_A_literature.md`. | §1 |
| `paper/` | `p3_article.tex` (Elsevier `elsarticle` class, preprint layout), `p3_refs.bib`, `elsarticle.cls` and `elsarticle-num.bst` (TeX Live 2022 copies), `p3_article_compiled.pdf`, `README_overleaf.txt`. Compiles with pdflatex + bibtex, or directly on Overleaf. | — |
| `figures/` | The five manuscript figures (`fig_wreath`, `fig_statistics`, `fig_reduction`, `fig_cost_law`, `fig_variants`; PDF + PNG) and `make_figures.py`, which regenerates all of them from the recorded outputs in `figures/data/` (`algo_results.json`, `solver_variants.json`, `round_algebra_results.json`, `wreath_summary.json`, `round_solver_profile.json`). Every plotted value is read from a record; nothing is typed by hand. | §3, §4.2, §5 |

The solver and cost records for the cryptographic application (Section 5) live in the companion repository
[`cpwi-post-quantum-dynamical-inversion`](https://github.com/gautierfilardo-efrei/cpwi-post-quantum-dynamical-inversion) (release v1.3.0: `wreath_solver.py`, `round_attack.py`, `round_results/`); they are not duplicated here.

## Quick checks

```bash
python figures/make_figures.py            # rebuilds the five figures into figures/ (needs numpy, matplotlib)
```

```bash
python3 hardness/verify_lemmaC.py          # affine criterion vs exhaustive conjugators, ~5 s
python3 hardness/verify_reduction.py       # case analysis, prime bound, 180 N3DM instances, ~1 min
python3 structure/round_algebra.py --help  # wreath normal form and special cases (exhaustive at m=5)
python3 solver/round_solver_study.py --help
```
`structure/algo_track.py` regenerates the exhaustive `N ≤ 9` tables (numpy; the `N = 9` stage takes about an hour).

## What is proved, verified, open

- **Proved**: short-word classification (x-length ≤ 2 polynomial; positive x-length 3 ≡ P3); wreath normal form; bounded-support, commuting-solution and involution-pair lemmas; first-moment identity; **Conj₃ NP-complete**, hence P3 with prescribed cycle type `3^{N/3}` NP-complete and no polynomial per-type enumeration unless P = NP.
- **Verified numerically** (scripts above): every lemma on small instances, the reduction against brute force, the cost law of the propagation solver.
- **Open**: the complexity of P3 itself. On the hard family of the theorem P3 is trivially solvable (commuting solution), so hardness must come from solutions that move the constant.

## Citation and licence

Repository: https://github.com/gautierfilardo-efrei/order3-conjugacy — cite the release used (see `CITATION.cff`; a Zenodo DOI is minted for each release once the repository is linked to Zenodo). Code under the MIT licence (`LICENSE`). Generative AI tools assisted with code development and text drafting; all numerical values come from the archived programs and their recorded outputs, and the author reviewed the mathematics, the sources and the interpretation.
