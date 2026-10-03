"""Random-coverage model: reveal k images and run propagation.
(i) correct images (from the planted x): how many further images are forced?  -> onset of the cascade
(ii) random block-respecting injective images: probability that propagation refutes them -> pruning depth
Then predicted full tree size  T(m) = sum_d m^(d)_falling * survive(d)  compared with the measured fill ratios."""
import os, sys, json, random, math, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from round_solver_study import *
rng = random.Random(12345)
out = {}
for m in range(6, 17):
    N = 2 * m
    res = {}
    for k in range(1, min(9, m) + 1):
        forced = []; refuted = 0; T = 300
        for _ in range(T):
            inst = wreath_instance(m, rng.randrange(1000))
            s = WordSolver(inst['a'], inst['c'], m)
            s.x = [NA] * N; s.xi = [NA] * N; s.trail = []; s.prop_passes = 0; s.ctype_prunes = 0
            pts = rng.sample(range(N), k)
            ok = all(s.assign(j, inst['x'][j]) for j in pts) and s.propagate()
            assert ok
            forced.append(len(s.trail) - k)
            # random wrong assignment on the same points, injective, block-respecting
            s.undo(0)
            ok = True
            used = set()
            for j in pts:
                lo, hi = (m, N) if j < m else (0, m)
                cands = [v for v in range(lo, hi) if v not in used]
                v = rng.choice(cands); used.add(v)
                ok = ok and s.assign(j, v)
            ok = ok and s.propagate()
            if not ok: refuted += 1
        res[k] = dict(mean_forced_correct=round(statistics.mean(forced), 2),
                      frac_cascade_complete=round(sum(1 for f in forced if f == N - k) / T, 3),
                      frac_forced_zero=round(sum(1 for f in forced if f == 0) / T, 3),
                      p_refuted_random=round(refuted / T, 3))
    out[m] = res
    print(m, {k: (r['mean_forced_correct'], r['frac_cascade_complete'], r['p_refuted_random']) for k, r in res.items()}, flush=True)
json.dump(out, open('results/onset_model.json', 'w'), indent=1)
