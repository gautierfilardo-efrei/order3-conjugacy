# Track "Lemmes" — full proofs of lem:cent and lem:support

Deliverables: `lemmas_full.tex` (drop-in replacement for the two sketches in `p3_article.tex`, uses the paper's macros
`\Sym \ctype \Cent \supp \Pthree \ZZ` and the environments `lem`/`rem`/`proof`; it references `lem:conj`, which exists in the paper),
`lemmas_track.py` (criterion, algorithms, brute force, campaigns), `lemmas_results.json`, `lemmas_run.log`.

Conventions: image tuples, products right to left, `(pq)(j)=p(q(j))`; `x a x^2` applied to `j` is `x(a(x(x(j))))`.

## 1. lem:cent (solutions commuting with a) — PROVED, for arbitrary a

Statement (now for every `a`, both `3 | L` and `3 ∤ L`). Let `a` have `m_L` cycles of length `L`, `C(a) = ∏_L Z_L ≀ Sym(m_L)`.
Solutions commuting with `a` = cube roots of `t = a^{-1} c` in `C(a)`; need `c ∈ C(a)`. For each `L`, `t` permutes the
`L`-cycles of `a` (permutation `ρ_L`); each `ρ_L`-cycle `Δ` of length `ℓ` has a rotation `R(Δ) ∈ Z_L` defined by `t^ℓ|_C = a^{R}|_C`
(`C ∈ Δ`), equal to the sum of the rotation coordinates along `Δ` (independent of base points). Class `(L, ℓ, R)` is *free* if
`3 ∤ ℓ` and (`3 ∤ L` or `R ≡ 0 mod 3`), *bound* otherwise.

* (i) cube root exists in `C(a)` ⟺ every bound class has size divisible by 3.
* (ii) number of commuting solutions = `∏_L ∏_{(ℓ,R)} Σ_k n!/(k! 6^k (n−3k)!) (2ℓ²L²)^k σ^{n−3k}`, `n = n_L(ℓ,R)`,
  `σ = 0` (bound), `1` (free, `3 ∤ L`), `3` (free, `3 | L`).
* (iii) O(N) after cycle decompositions ⇒ polynomial.

Proof ingredients (all written out in the .tex):
1. Wreath coordinates `g = (v, π)`, `g(a^u p_i) = a^{u+v_i} p_{π(i)}`; composition `(w,σ)(v,π) = (v + w∘π, σπ)`;
   cube `(v + v∘π + v∘π², π³)`. Direct product over `L` ⇒ componentwise.
2. Cube roots of `ρ` in `Sym(m)` (Lemma `lem:cuberoots`, also written out): partition of ρ-cycles into singletons
   (`3 ∤ ℓ`, `π = ρ^k`, `3k ≡ 1`) and triples of equal length (`2ℓ²` roots each). Count formula for cube roots of a permutation.
3. The vector equation splits along π-cycles into circulant systems `u_j + u_{j+1} + u_{j+2} = s_j` over `Z_L`.
   * Triple (`λ = 3ℓ`): solvable ⟺ the three rotations are equal; then exactly `L²` solutions (telescoping
     `u_{j+3} − u_j = s_{j+1} − s_j`). **No use of gcd(L,3)=1** — the sketch's restriction was unnecessary here.
   * Singleton (`λ = ℓ`, `3 ∤ ℓ`): the system is multiplication by `Φ = 1+X+X²` in `Z_L[X]/(X^ℓ−1)`; `Res(Φ, X^ℓ−1) = 3`.
     If `3 ∤ L`: `Φ` invertible, unique solution. If `3 | L`: CRT `Z_L = Z_{3^e} × Z_{L'}`; over `Z_{3^e}`, Hensel splits
     `X^ℓ−1 = (X−1)Ψ`, `Φ ↦ 3` on the `X=1` factor and `Φ` is a unit mod `Ψ` (since `Φ ≡ (X−1)² mod 3` is coprime to `Ψ`);
     hence solvable ⟺ `R(Δ) ≡ 0 mod 3`, then exactly 3 solutions. (Necessity also by summing the equations: `R = 3Σu`.)
4. Classes are independent ⇒ criterion (i) and count (ii).

Special case of the old statement (`3 ∤ L`): bound classes = those with `3 | ℓ` — recovers the printed criterion.

VERIFIED (`lemmas_track.py`, function `verify_cent`, `verify_cent_blocks`; `lemmas_results.json`):
* All cycle types `a` for `N ≤ 8` (66 classes), all `c ∈ C(a)` (49 906 pairs, 43 206 at N=8): count from (ii) vs. enumeration of
  `C(a)` (Counter of `x³`): 0 mismatches; max count 1233 (identity `a`, N=8). Cross-check against brute force over all of
  `Sym(N)` with the commuting constraint, N ≤ 7: 1 248 pairs (including non-commuting `c`), 0 mismatches.
* Blockwise families with `3 | L`: `(L,m) ∈ {(3,1..6),(6,1..4),(9,1..3),(12,2)}`, i.e. N up to 27, `|C(a)|` up to 524 880, every
  `t ∈ C(a)` (the sketch was silent on these): 0 mismatches.

## 2. lem:support (bounded support) — PROVED

Route: `x a x² = c ⟺ c = b x³`, `b = x a x^{-1}` determined by `f = x|_S`, `S = supp a`, `|S| = s`: `b_f(f(p)) = f(a(p))`,
`supp b_f = f(S)`. Conversely for any injection `f`, `{x : x³ = t_f, x|_S = f}` with `t_f = b_f^{-1} c` is exactly the set of solutions
with `x|_S = f` (one checks `x a x^{-1} = b_f` from `x|_S = f`). So solutions = disjoint union over ≤ `N^s` injections.

Lemma `lem:prescribed` (new, PROVED): cube roots of `t` with prescribed images on `S` are counted in O(N). Constraint graph `Γ`
on the cycles of `t`, edge `D(p) → D(f(p))` with a phase; necessary conditions: equal lengths, one edge per source cycle, in-degree
≤ 1. Then components are paths/cycles inside `x̄`-orbits of size 1 or 3, hence of exactly four admissible shapes:
loop (singleton; `3 ∤ ℓ` and phase `≡ 3^{-1} mod ℓ`), directed 3-cycle (phases sum to `1 mod ℓ`), path with 2 edges (third map
forced), path with 1 edge (needs a *free* third cycle of the same length: `(n_ℓ)_{q_ℓ} ℓ^{q_ℓ}` choices). Remaining free cycles:
unconstrained cube-root count (`eq:cuberootcount`, needs `3 | (n_ℓ − q_ℓ)` for `3 | ℓ`). Formula `eq:prescribedcount`.
Total complexity O(N^{s+1}).

VERIFIED (`verify_prescribed`, `verify_support`, `verify_support_extra`):
* `lem:prescribed` directly: 396 random `(t, f)` with N ≤ 8, biased towards cycle lengths divisible by 3 and towards
  consistent prescriptions (117 with ≥ 1 root): 0 mismatches vs brute force over `Sym(N)`.
* `lem:support` exhaustive: N=5 and N=6, all `a` with `|supp a| ∈ {2,3,4}` (10+20+45, 15+40+135 permutations `a`), **all** `c`
  (146 046 instances) ; N=7, s=2, all 21 `a`, all 5040 `c` (105 840 instances): 0 mismatches of the *count* vs brute force.
* `lem:support` random: N=7,8,9 × s=2,3,4: 300/300/150 instances per (N,s) plus the first run's 40/30/12, a third planted with
  random `x`, a third random `c`, a third planted with `x` of 3-divisible cycle type (half of those perturbed by a 3-cycle):
  2 496 instances (246 + 2 250), 0 mismatches.

## 3. Deviations / caveats
* No TeX binary in the sandbox: `lemmas_full.tex` was **not** compiled here. It uses only `amsmath/amsthm` constructs and the
  paper's macros; label `lem:cuberoots`, `lem:prescribed`, equations `eq:cuberootcount`, `eq:prescribedcount`, `eq:cubeeq`,
  `eq:circ`, `eq:centcount` are new — check for clashes when pasting into v7. `\ref{lem:conj}` must resolve (it does in v6).
* Statement change vs v6: lem:cent now covers arbitrary `a`, with the extra condition `R ≡ 0 mod 3` for singletons when `3 | L`
  (missing from the sketch, which only treated `gcd(L,3)=1`); lem:support now gives the *count* and the explicit bound O(N^{s+1}).
  The sketch's parenthetical "(Verified exhaustively for s=2, N=6: 30/30)" should be replaced by the new Remark.
* Section 5 sentence "x acting on every block as the translation by (m_i−1)/3 commutes with a ... (Lemma lem:cent)" remains valid.
* Brute force at N=9 is over 9! = 362 880 permutations (numpy); random instances only (1 per ~0.1 s), not exhaustive, as planned.
