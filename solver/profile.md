# Profile of the propagation solver for x a x² = c in A_m wr C₂ — why "≈2 bits per m"

Setting. Round system in wreath normal form: unknown x = (U,V)·σ ∈ S_{2m}, constant a = (B,A),
target c = (U',V')·σ; equation x a x² = c, block structure enforced (x maps block 0 → block 1 → block 0).
Solver = the standard chain-propagation DFS (as `word_solver.solve_word` / `round_attack.solve_round`):
for every point j the chain  j → x(j) → x²(j) → a x²(j) → x(a x²(j)) = c(j)  contains three
applications of x; a chain forces its missing letter as soon as the other two are known; branching
assigns one image x(j) = v (domain = the free images in the opposite block). The instrumented,
array-based, trail-undo re-implementation is `round_solver_study.WordSolver` (≈ 70 000 nodes/s in
CPython); node counts agree with `solve_word` filtered to the block structure at m = 5 (exhaustive
solution sets identical on 4/4 instances). Numbers below: `round_solver_profile.json`
(5 planted instances per m, seeds SHA-256("CPWI-SOLV-instance|m|t"), full enumeration of the search
tree so that the planted branch is always traversed) and `results/onset_model.json`.

## 1. What the instrumentation shows (heuristic "first free point", m = 5…14)

| m | full tree nodes (median) | nodes to 1st solution | planted depth k (5 instances) | mean domain on planted path | β = nodes^{1/k} | β/m |
|---|---|---|---|---|---|---|
| 6 | 76 | 10 | 3,4,4,3,2 | 4.5 | 3.9 | 0.65 |
| 8 | 526 | 75 | 4,3,4,3,3 | 6.4 | 6.3 | 0.79 |
| 10 | 5 261 | 3 413 | 3,4,5,4,6 | 7.9 | 7.0 | 0.70 |
| 12 | 75 888 | 24 720 | 6,4,4,5,5 | 9.6 | 10.4 | 0.87 |
| 14 | 763 233 | 401 048 | 7,7,5,5,7 | 10.8 | 8.9 | 0.64 |

* **Number of free decisions.** The planted solution is reached after k(m) branching decisions,
  all the other 2m − k images being forced by propagation. Fit over m = 5…14:
  **k(m) ≈ 0.38 m + 0.6** (heuristic "first"), 0.31 m + 0.9 (heuristic "chains"). k grows linearly, not
  logarithmically, in m.
* **The first decision never propagates** (0 forced images in 42/45 planted paths at m ≥ 6, 2 forced in the 3 others), the second forces nothing in 70 % of the planted paths at m ≥ 9 with
  heuristic "first" (median 0, range 0–3; at m = 6–8 the small groups give median 2, range 0–12),
  whereas heuristic "chains" — which deliberately branches on x(v) right after x(j) = v so that a chain
  gets two known letters — always forces at least one image (m ≥ 9: median 1, range 1–4; m = 6–8:
  median 3, range 1–10). First-decision zero-propagation counts: 42/45 ("first"), 52/55 ("chains"). The cascade that completes the assignment happens at decision 3–7 and forces 10–20 images
  at once: propagation on this word is *all-or-nothing*.
* **Domains stay ≈ m.** On the planted path the domain of the d-th decision is m − (d−1) minus the few
  images forced so far: mean 4.5 (m=6) → 10.8 (m=14). The branching factor is not a constant; it is Θ(m).
* **Effective branching.** β(m) = nodes^{1/k(m)} lies between 0.6 m and 0.87 m: the tree is close to
  complete down to depth k − 1. The fill ratio (nodes at depth d ÷ m(m−1)…(m−d+1)) is
  1.0 / 0.98 / 0.93 / 0.65 / 0.30 / 0.036 at depths 1…6 for m = 12, and *increases with m at fixed
  depth* (depth 4: 0.14 at m=8, 0.38 at m=10, 0.65 at m=12, 0.82 at m=14): pruning starts later for
  larger m.
* **Where the nodes are.** 90–95 % of all nodes lie in subtrees hanging off the *root's* wrong
  children (wasted fraction at planted depth 0: 0.71 at m=5 → 0.92 at m=14), 5–9 % off the second
  planted decision, < 1 % deeper: once two or three correct images are fixed, propagation finishes the
  job; the cost is the exhaustive refutation of the wrong choices for the first 2–4 images.
* Cycle-type dependence of the constants: not observable at this sample size (5 instances/m; the
  per-m spread of node counts is a factor 3–10, e.g. m=13: 22 640…139 232), consistent with the
  algebra notes (difficulty depends on the classes of A, B only through polynomial factors).

## 2. A random-coverage model that reproduces the tree

Reveal k images of the planted x at random points and propagate (`onset_model.py`, 300 trials/point):

| m | k=2: mean forced / P(cascade complete) | k=3 | k=4 | k=5 | k=6 | k=7 |
|---|---|---|---|---|---|---|
| 8 | 1.2 / 0.01 | 4.4 / 0.18 | 9.0 / 0.65 | 10.1 / 0.85 | 9.4 / 0.89 | 8.8 / 0.95 |
| 12 | 0.6 / 0.00 | 1.9 / 0.01 | 5.2 / 0.11 | 11.9 / 0.47 | 14.1 / 0.68 | 15.5 / 0.82 |
| 16 | 0.5 / 0.00 | 1.0 / 0.00 | 2.3 / 0.00 | 6.1 / 0.06 | 13.2 / 0.36 | 17.7 / 0.59 |

and, with k *random* (wrong, injective, block-respecting) images, the probability that propagation
refutes them: P_ref(k) = 0.05 / 0.17 / 0.49 / 0.82 / 0.97 / 0.997 for k = 2…7 at m = 12.
The predicted tree size  T(m) = Σ_d m^{(d)} (1 − P_ref(d))  (m^{(d)} falling factorial) gives
T(8) ≈ 4.8·10² against the measured 5.3·10², T(10) ≈ 4.7·10³ against 5.3·10³, T(12) ≈ 5.6·10⁴ against
7.6·10⁴, T(14) ≈ 6.6·10⁵ against 7.6·10⁵ — the model is quantitatively right; no property of the solver beyond "chain propagation
with 3 letters per chain" is needed.

Why it behaves this way (in words). Each chain constraint involves three images of x, located at
j, x(j) and a(x²(j)); the last two locations themselves depend on x. A single known image never
completes two letters of any chain, so the first decision is free of consequences; k random known
images cover two letters of a given chain with probability ≈ 3(k/N)² (N = 2m), so the expected number
of forced images is ≈ 3k²/N — below one image for k < √(N/3) ≈ 2–3, and the cascade to full
determination needs the forced images to feed new chains, which (empirically) happens at
k* ≈ 0.38 m + 0.6 (random-coverage model: the median cascade-completion seed is 3,4,4,4,5,5,6,6,6,7,7 images for m = 6…16), i.e. when a constant fraction ≈ 1/5 of the 2m images is known. A *wrong* partial
assignment is refuted only through a fully known chain or an injectivity clash, which needs the
same coverage, so wrong branches survive to the same depth as the planted one. The tree is therefore
essentially the complete tree of depth k*(m) with branching m − d at depth d:

  nodes(m) ≈ fill · m^{k*(m)}  ⇒  log₂ nodes ≈ (0.38 m + 0.6)·log₂(0.7 m) + O(1).

This is an m log m law (2^{Θ(m log m)} = |A_m|^{Θ(1)}, a constant fraction ≈ 1/5 of the brute-force
exponent 2·log₂ m!), not a 2^{cm} law. Its local slope d/dm ≈ 0.38·log₂(0.7 m) + 0.55 is 1.5 (m=6),
1.7 (m=10), 1.9 (m=14), 2.3 (m=36), 2.6 (m=59) bits per unit of m: the "≈2 bits per m" measured on
m = 6…13 is the tangent of this curve in that window, and a straight-line extrapolation of it to
m = 36…67 *under*-estimates the cost (the measured half-window slopes confirm the curvature for the
full-tree counts: 1.53 bits/m over m = 5…9, 1.77 bits/m over m = 9…14, 1.68 overall; the noisier first-solution
medians on 5–9 instances do not resolve it). "Two bits per point undetermined" is thus the wrong
picture: propagation leaves not 2 bits per point but ≈ 0.4 m free *images* (each worth log₂ m bits),
because a length-3 word gives every point one constraint that only bites when two of its three
letters are already known — the constraint hypergraph has density 1 (N chains, N unknowns, 3 letters
each), exactly at the peeling threshold, where the cascade needs a linear-size seed.

## 3. Consequences for the variants (details in `solver_variants.json`)

* Better branching orders (b1, b2) reduce k(m) slightly (0.31 m instead of 0.38 m) and mainly reduce
  the fill (β/m ≈ 0.4–0.6 instead of 0.6–0.85): constant-factor gains of 10–50× at m = 13, with a
  measurably smaller slope in the window (≈1.0–1.3 vs 1.8 bits/m) but the same m log m mechanism.
* The cycle-type relaxation ctype(a x³) = ctype(c) never prunes a node (0 prunings in every run):
  c = x(a x³)x⁻¹ is checked point-wise by chain propagation, so every closed cycle of the partial a x³
  is already mapped by x onto a cycle of c of the same length, and the multiplicity constraints
  follow from injectivity. The invariant is a consequence of the local constraints, not extra
  information.
* Z3 (one-hot SAT encoding with clause learning) has the smallest measured slope (≈1.2 bits/m in
  seconds, m = 5…16) but is still exponential; the conjugated form (d) reduces to enumerating the
  S_{2m}-class of c (size |S_{2m}|/|C(c)|, 6.6·10⁶ at m = 6) and is not competitive.
