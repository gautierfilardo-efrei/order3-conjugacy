"""Gadget b track: hardness of  x a x^2 = c  (P3) with a an involution, via the choice of b = x a x^{-1}.

Gadget "B1" (block = one c-cycle of length L+1 carrying one transposition of a):
  E = index set of blocks (|E| = n), each block e has points g_e[0..L], c(g_e[i]) = g_e[i+1 mod L+1],
  a = prod_e (g_e[0]  g_e[Delta_e])  with Delta_e = L - d_e, label d_e in Z_L.
Intended solutions: t = x^3 = b c has, inside block e, one fixed point v_e = g_e[Delta_e] (= alpha'_e) and one
L-cycle B_e starting at u_e = c(v_e); z = induced permutation of x on the blocks; z has only 3-cycles
(and, when 3 does not divide L, fixed points on blocks with d_e = -3^{-1} mod L); every 3-cycle {e,f,g} of z
must satisfy  d_e + d_f + d_g = -1 (mod L).  So P3 on these instances is a "3-partition mod L".

Everything here is standard library.  The propagation solver word_solver.py (copied from the project artifact)
is used for complete enumeration when N > 9; for N <= 9 we enumerate S_N.
"""
from __future__ import annotations
import itertools, json, sys, time
from word_solver import compose, inverse, mul, evaluate, solve_word, cycle_type

def build_instance(labels, L):
    """labels: list of d_e in Z_L.  Returns (a, c, meta)."""
    n = len(labels); N = n * (L + 1)
    a = list(range(N)); c = list(range(N)); blocks = []
    for e, d in enumerate(labels):
        d %= L
        g = [e * (L + 1) + i for i in range(L + 1)]
        for i in range(L + 1): c[g[i]] = g[(i + 1) % (L + 1)]
        Delta = L - d                      # in [1, L]
        a[g[0]], a[g[Delta]] = g[Delta], g[0]
        blocks.append({"e": e, "d": d, "alpha": g[0], "alpha_prime": g[Delta], "points": g})
    return tuple(a), tuple(c), {"L": L, "n": n, "N": N, "blocks": blocks}

def perms_with_cycles_1_and_3(n, allow_fixed):
    """All z in S_n whose cycles have length 3, or length 1 on indices in allow_fixed."""
    out = []
    def rec(remaining, z):
        if not remaining:
            out.append(tuple(z)); return
        e = remaining[0]; rest = remaining[1:]
        if e in allow_fixed:
            z[e] = e; rec(rest, z)
        for f, g in itertools.permutations(rest, 2):
            z[e], z[f], z[g] = f, g, e
            rec([h for h in rest if h not in (f, g)], z)
        z[e] = -1
    rec(list(range(n)), [-1] * n)
    return out

def predicted_solutions(labels, L, a, c, meta):
    """Construct the intended solutions x from admissible block permutations z and verify them."""
    n = len(labels); N = meta["N"]
    r = None
    if L % 3 != 0:
        r = pow(3, -1, L) if L > 1 else 0
    allow_fixed = set()
    if r is not None:
        allow_fixed = {e for e, d in enumerate(labels) if (d + r) % L == 0}
    sols = []; zs = []
    for z in perms_with_cycles_1_and_3(n, allow_fixed):
        # check sums on 3-cycles
        ok = True; seen = set()
        for e in range(n):
            if e in seen: continue
            if z[e] == e: seen.add(e); continue
            f, g = z[e], z[z[e]]; seen |= {e, f, g}
            if (labels[e] + labels[f] + labels[g] + 1) % L != 0: ok = False; break
        if not ok: continue
        # build x
        x = [-1] * N
        # t = b c is not known before b, but inside block e it is: fixed point v_e, L-cycle u_e -> c(u_e) -> ...
        for e in range(n):
            blk = meta["blocks"][e]; g = blk["points"]; Delta = L - blk["d"]
            v = g[Delta]; u = g[(Delta + 1) % (L + 1)]
            Bcyc = [u]
            while len(Bcyc) < L: Bcyc.append(c[Bcyc[-1]])     # u, c u, ..., c^{L-1} u  (t = c here except last step)
            f = z[e]
            blkf = meta["blocks"][f]; gf = blkf["points"]; Deltaf = L - blkf["d"]
            vf = gf[Deltaf]; uf = gf[(Deltaf + 1) % (L + 1)]
            Bf = [uf]
            while len(Bf) < L: Bf.append(c[Bf[-1]])
            x[v] = vf
            if f == e:
                for i in range(L): x[Bcyc[i]] = Bcyc[(i + r) % L]
            else:
                # x(alpha_e) = u_f with alpha_e = t^{d_e} u_e, and x t = t x
                for i in range(L): x[Bcyc[(blk["d"] + i) % L]] = Bf[i % L]
        x = tuple(x)
        assert sorted(x) == list(range(N)), "constructed x is not a permutation"
        assert evaluate(['x', a, 'x', 'x'], x) == c, "constructed x fails x a x^2 = c"
        sols.append(x); zs.append(z)
    return sols, zs

def all_solutions(a, c, node_budget=5 * 10**6):
    N = len(c)
    if N <= 9:
        sols = [x for x in itertools.permutations(range(N)) if evaluate(['x', a, 'x', 'x'], x) == c]
        return sols, {"method": "brute", "complete": True}
    sols, st = solve_word(['x', a, 'x', 'x'], c, max_solutions=10**7, node_budget=node_budget)
    return sols, {"method": "propagation", "complete": not st["exhausted"], "nodes": st["nodes"]}

def test_instance(labels, L, node_budget=5 * 10**6):
    a, c, meta = build_instance(labels, L)
    pred, zs = predicted_solutions(labels, L, a, c, meta)
    t0 = time.time(); sols, info = all_solutions(a, c, node_budget); dt = time.time() - t0
    S = set(sols); P = set(pred)
    extra = sorted(S - P); missing = sorted(P - S)
    rec = {"L": L, "labels": list(labels), "N": meta["N"], "n_solutions": len(S), "n_predicted": len(P),
           "unintended": len(extra), "missing_predicted": len(missing), "complete": info["complete"],
           "method": info["method"], "time_s": round(dt, 3), "pass": (not extra and not missing and info["complete"]),
           "z_list": [list(z) for z in zs]}
    if extra:
        rec["unintended_examples"] = [{"x": list(x), "ctype_x": cycle_type(x),
                                      "b": list(mul(x, a, inverse(x))),
                                      "ctype_t": cycle_type(mul(x, x, x))} for x in extra[:5]]
    return rec

def run_grid(grid, node_budget=5 * 10**6, verbose=True):
    results = []
    for (L, n) in grid:
        for labels in itertools.product(range(L), repeat=n):
            rec = test_instance(labels, L, node_budget)
            results.append(rec)
            if verbose:
                print(f"L={L} n={n} labels={labels} N={rec['N']} sols={rec['n_solutions']} pred={rec['n_predicted']} "
                      f"extra={rec['unintended']} miss={rec['missing_predicted']} complete={rec['complete']} "
                      f"{rec['time_s']}s", flush=True)
    return results

from collections import Counter

def predicted_solutions2(labels, L, a, c, meta, allow_flip=True):
    """Oriented family F(labels, L): per-block orientation eps_e (which a-endpoint becomes the t-fixed point),
    effective label d_e (eps=0) or L-1-d_e (eps=1); z in S_E with 3-cycles (and fixed points on blocks whose
    effective label is -3^{-1} mod L when 3 does not divide L); 3-cycle sums = -1 mod L.  Every member is
    verified to satisfy x a x^2 = c (assertion below)."""
    n = len(labels); N = meta["N"]
    r = pow(3, -1, L) if (L % 3 != 0 and L > 1) else (0 if L == 1 else None)
    sols = []; descr = []
    for eps in itertools.product(range(2) if allow_flip else [0], repeat=n):
        eff = [(labels[e] if eps[e] == 0 else (L - 1 - labels[e])) % L for e in range(n)]
        allow_fixed = {e for e in range(n) if r is not None and (eff[e] + r) % L == 0}
        for z in perms_with_cycles_1_and_3(n, allow_fixed):
            ok = True; seen = set()
            for e in range(n):
                if e in seen: continue
                if z[e] == e: seen.add(e); continue
                f, g = z[e], z[z[e]]; seen |= {e, f, g}
                if (eff[e] + eff[f] + eff[g] + 1) % L != 0: ok = False; break
            if not ok: continue
            x = [-1] * N; data = []
            for e in range(n):
                blk = meta["blocks"][e]
                v = blk["alpha_prime"] if eps[e] == 0 else blk["alpha"]
                u = c[v]; B = [u]
                while len(B) < L: B.append(c[B[-1]])
                data.append((v, B))
            for e in range(n):
                v, B = data[e]; f = z[e]; vf, Bf = data[f]
                x[v] = vf
                if f == e:
                    for i in range(L): x[B[i]] = B[(i + r) % L]
                else:
                    for i in range(L): x[B[(eff[e] + i) % L]] = Bf[i]
            x = tuple(x)
            assert sorted(x) == list(range(N)) and evaluate(['x', a, 'x', 'x'], x) == c, (labels, L, eps, z)
            sols.append(x); descr.append((eps, z))
    return sols, descr

def cycles_of(p):
    seen = set(); out = []
    for i in range(len(p)):
        if i in seen: continue
        cyc = []; j = i
        while j not in seen: seen.add(j); cyc.append(j); j = p[j]
        out.append(cyc)
    return out

def shape(x, a, c, L):
    """(does some b-transposition cross blocks?, per-block sorted lengths of the t-cycles contained in the block)."""
    N = len(c); blk = lambda j: j // (L + 1)
    b = mul(x, a, inverse(x)); t = mul(x, x, x)
    cross = any(blk(i) != blk(b[i]) for i in range(N) if b[i] != i)
    parts = tuple(sorted(tuple(sorted(len(cy) for cy in cycles_of(t) if all(blk(j) == blk(cy[0]) == e for j in cy)))
                         for e in range(N // (L + 1))))
    return cross, parts

def test_instance3(labels, L, node_budget=5 * 10**6):
    a, c, meta = build_instance(labels, L)
    pred, descr = predicted_solutions2(labels, L, a, c, meta)
    t0 = time.time(); sols, info = all_solutions(a, c, node_budget); dt = time.time() - t0
    S = set(sols); P = set(pred); extra = sorted(S - P); missing = sorted(P - S)
    shapes = Counter(shape(x, a, c, L) for x in extra)
    fp_every_block = sum(1 for x in extra if all(1 in parts for parts in shape(x, a, c, L)[1]) and not shape(x, a, c, L)[0])
    return {"L": L, "labels": list(labels), "N": meta["N"], "n_solutions": len(S), "n_family_oriented": len(P),
            "n_family_unflipped": sum(1 for d in descr if all(v == 0 for v in d[0])),
            "unintended": len(extra), "missing_predicted": len(missing), "complete": info["complete"],
            "method": info["method"], "nodes": info.get("nodes"), "time_s": round(dt, 3),
            "unintended_with_fixed_point_in_every_block": fp_every_block,
            "unintended_shapes": [{"cross_block_b": k[0], "t_parts_per_block": [list(p) for p in k[1]], "count": v}
                                  for k, v in sorted(shapes.items(), key=str)],
            "family_pass": (not missing and info["complete"]),
            "gadget_pass": (not extra and not missing and info["complete"])}

if __name__ == "__main__":
    quick = "--quick" in sys.argv
    grid = [(2, 3), (3, 3)] + ([] if quick else [(4, 3)])
    results = []
    for (L, n) in grid:
        for labels in itertools.product(range(L), repeat=n):
            rec = test_instance3(labels, L); results.append(rec)
            print(f"L={L} labels={labels} N={rec['N']} sols={rec['n_solutions']} family={rec['n_family_oriented']} "
                  f"unintended={rec['unintended']} complete={rec['complete']} gadget_pass={rec['gadget_pass']}", flush=True)
    json.dump(results, open("gadget_b_results.json", "w"), indent=1)
