import sys, json, random, time
sys.path.insert(0, '.')
import fourlayer_track as FT
from math import factorial
rng = random.Random(2024)
out = {}
for k, nsamp in [(1, 400000), (2, 1500000)]:
    a, c = FT.affine_instance(5, tuple([2]*k)); ap, cp = FT.pad(a, c, (1,2,3,3))
    N = len(ap); target = FT.ctype(ap)
    hits = 0; t0 = time.time()
    for _ in range(nsamp):
        x = list(range(N)); rng.shuffle(x)
        # a' x^3
        y = [ap[x[x[x[j]]]] for j in range(N)]
        if FT.ctype(tuple(y)) == target: hits += 1
    p = hits / nsamp
    class_size = factorial(N) // (20**k * factorial(k))
    mean = p * factorial(N) / class_size   # = p * 20^k k!
    out[f"k{k}"] = {"N": N, "samples": nsamp, "hits": hits, "p_hat": p, "mean_solutions_over_class": mean,
                    "stderr_mean": (p*(1-p)/nsamp)**0.5 * 20**k * factorial(k), "time_s": round(time.time()-t0,1)}
    print(out[f"k{k}"], flush=True)
    json.dump(out, open("first_moment_mc.json","w"), indent=1)
