#!/usr/bin/env python3
"""lemmas_track.py -- full proofs' companion for lem:cent and lem:support of the P3 paper.

Conventions: permutations are image tuples, products right to left, (pq)(j) = p[q[j]].
x a x^2 applied to j is x(a(x(x(j)))).

Part A (lem:cent).  Solutions of x a x^2 = c commuting with a  <=>  cube roots of t = a^{-1}c in C(a).
  C(a) = prod_L  Z_L wr Sym(m_L).  Element (v, pi): x(a^k p_i) = a^{k+v_i} p_{pi(i)}  (p_i base points).
  Composition: (w,sigma)(v,pi) = (v + w o pi, sigma pi);  cube: (v + v o pi + v o pi^2, pi^3).
  Criterion (proved in lemmas_full.tex): group the rho-cycles (rho = cycle permutation induced by t) by
  (length ell, rotation sum R in Z_L).  A class may be used as singletons iff 3 ∤ ell and (3 ∤ L or R ≡ 0 mod 3);
  otherwise its size must be divisible by 3.  Count = product of per-class sums.

Part B (lem:support).  |supp a| = s.  Guess f = x|_S (<= N^s injections), b := f a f^{-1} (support f(S)),
  t := b^{-1} c; count cube roots x of t with x|_S = f  (constraint graph on the cycles of t: components are
  loops, 3-cycles, paths of 1 or 2 edges; paths of 1 edge need a free cycle of the same length).

Run:  python3 lemmas_track.py [--quick]
"""
import sys, json, time, random, itertools
from math import factorial, gcd
from collections import Counter, defaultdict
import numpy as np

# ------------------------------------------------------------------ permutation utilities
def comp(p, q):
    return tuple(p[q[j]] for j in range(len(p)))

def inv(p):
    r = [0] * len(p)
    for j, pj in enumerate(p):
        r[pj] = j
    return tuple(r)

def order(p):
    o = 1
    for c in cycles(p):
        o = o * len(c) // gcd(o, len(c))
    return o

def cycles(p):
    n = len(p); seen = [False] * n; out = []
    for j in range(n):
        if not seen[j]:
            c = []; k = j
            while not seen[k]:
                seen[k] = True; c.append(k); k = p[k]
            out.append(c)
    return out

def ctype(p):
    return tuple(sorted((len(c) for c in cycles(p)), reverse=True))

def identity(n):
    return tuple(range(n))

def rand_perm(n, rng):
    l = list(range(n)); rng.shuffle(l); return tuple(l)

def support(p):
    return [j for j in range(len(p)) if p[j] != j]

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

def perm_of_type(lam, n):
    p = list(range(n)); pos = 0
    for L in lam:
        for j in range(L):
            p[pos + j] = pos + (j + 1) % L
        pos += L
    return tuple(p)

# ------------------------------------------------------------------ Part A: centraliser criterion
def cube_count_class(n, sigma, tau):
    """Number of ways to split n labelled objects into singletons (weight sigma) and unordered triples
    (weight tau), summed: sum_k n!/(k! 6^k (n-3k)!) tau^k sigma^(n-3k).  sigma=0 forbids singletons."""
    tot = 0
    for k in range(n // 3 + 1):
        rest = n - 3 * k
        if sigma == 0 and rest > 0:
            continue
        tot += factorial(n) // (factorial(k) * 6 ** k * factorial(rest)) * tau ** k * (sigma ** rest if rest else 1)
    return tot

def cent_data(a, c):
    """Return None if c does not commute with a; else dict L -> list of (ell, R) over rho-cycles."""
    n = len(a)
    if comp(a, c) != comp(c, a):
        return None
    t = comp(inv(a), c)
    cyc = cycles(a)
    where = {}            # point -> (cycle index, exponent k with point = a^k(base))
    for i, cy in enumerate(cyc):
        for k, pt in enumerate(cy):
            where[pt] = (i, k)
    m = len(cyc)
    rho = [None] * m; r = [None] * m
    for i, cy in enumerate(cyc):
        j, k = where[t[cy[0]]]
        assert len(cyc[j]) == len(cy)
        rho[i] = j; r[i] = k
    out = defaultdict(list)
    for rc in cycles(tuple(rho)):
        L = len(cyc[rc[0]])
        R = sum(r[i] for i in rc) % L
        out[L].append((len(rc), R))
    return out

def cent_criterion(a, c):
    """(exists, count) of x in C(a) with x a x^2 = c, i.e. x^3 = a^{-1}c in C(a).  Proved criterion."""
    data = cent_data(a, c)
    if data is None:
        return False, 0
    total = 1
    for L, lst in data.items():
        classes = Counter(lst)
        for (ell, R), nn in classes.items():
            single_ok = (ell % 3 != 0) and (L % 3 != 0 or R % 3 == 0)
            sigma = (3 if L % 3 == 0 else 1) if single_ok else 0
            tau = 2 * ell * ell * L * L
            cnt = cube_count_class(nn, sigma, tau)
            if cnt == 0:
                return False, 0
            total *= cnt
    return True, total

def enumerate_centraliser(a):
    """Yield every x in C(a) as an image tuple (prod_L Z_L wr Sym(m_L))."""
    n = len(a)
    byL = defaultdict(list)
    for cy in cycles(a):
        byL[len(cy)].append(cy)
    blocks = []
    for L, cys in byL.items():
        m = len(cys)
        opts = []
        for pi in itertools.permutations(range(m)):
            for v in itertools.product(range(L), repeat=m):
                maps = []
                for i in range(m):
                    src = cys[i]; dst = cys[pi[i]]
                    for k in range(L):
                        maps.append((src[k], dst[(k + v[i]) % L]))
                opts.append(maps)
        blocks.append(opts)
    for choice in itertools.product(*blocks):
        x = [0] * n
        for maps in choice:
            for s_, d_ in maps:
                x[s_] = d_
        yield tuple(x)

def cube_counter_centraliser(a):
    """Counter: t -> #{x in C(a): x^3 = t}."""
    cnt = Counter()
    for x in enumerate_centraliser(a):
        cnt[comp(x, comp(x, x))] += 1
    return cnt

# ------------------------------------------------------------------ Part B: prescribed cube roots
def prescribed_cube_roots(t, f):
    """Number of x in Sym(N) with x^3 = t and x(p) = f[p] for all p in f (f a partial injection).
    Polynomial-time: constraint graph on the cycles of t."""
    n = len(t)
    cyc = cycles(t)
    cid = [0] * n; pos = [0] * n
    for i, cy in enumerate(cyc):
        for k, pt in enumerate(cy):
            cid[pt] = i; pos[pt] = k
    m = len(cyc)
    # edges: cycle i -> (cycle j, phase d) meaning x(t^k p_i) = t^{k+d} p_j  (p_i = cyc[i][0])
    out_edge = {}
    for p, q in f.items():
        i, j = cid[p], cid[q]
        if len(cyc[i]) != len(cyc[j]):
            return 0
        d = (pos[q] - pos[p]) % len(cyc[i])
        if i in out_edge and out_edge[i] != (j, d):
            return 0
        out_edge[i] = (j, d)
    in_edge = {}
    for i, (j, d) in out_edge.items():
        if j in in_edge:
            return 0          # two cycles would map onto the same cycle
        in_edge[j] = i
    touched = set(out_edge) | set(in_edge)
    # loops
    comp_free = Counter(len(cy) for i, cy in enumerate(cyc) if i not in touched)
    need_third = Counter()       # ell -> number of 1-edge paths needing a free third cycle
    visited = set()
    for i in touched:
        if i in visited:
            continue
        # walk back to the start of the component
        start = i
        seen_back = {start}
        while start in in_edge and in_edge[start] not in seen_back:
            start = in_edge[start]; seen_back.add(start)
        # collect the forward chain from start
        chain = [start]; phases = []
        cur = start
        is_cycle = False
        while cur in out_edge:
            j, d = out_edge[cur]
            if j == start:
                is_cycle = True; phases.append(d); break
            if j in chain:
                return 0      # cannot happen (in-degree <= 1) but be safe
            chain.append(j); phases.append(d); cur = j
        visited.update(chain)
        ell = len(cyc[start])
        if is_cycle:
            if len(chain) == 1:                        # loop: singleton cycle, x = t^k, 3k = 1 mod ell
                if ell % 3 == 0:
                    return 0
                k = pow(3, -1, ell) if ell > 1 else 0
                if phases[0] % ell != k % ell:
                    return 0
            elif len(chain) == 3:                      # full triple: composition of phases must be t
                if sum(phases) % ell != 1 % ell:
                    return 0
            else:
                return 0
        else:
            if len(chain) == 2:                         # one edge: need a free third cycle of length ell
                need_third[ell] += 1
            elif len(chain) == 3:                       # two edges: third map determined, always consistent
                pass
            else:
                return 0
    total = 1
    for ell in set(comp_free) | set(need_third):
        nfree = comp_free[ell]; q = need_third[ell]
        if nfree < q:
            return 0
        ways = 1
        for j in range(q):
            ways *= (nfree - j) * ell                   # choose the free third cycle, then its phase
        rest = nfree - q
        sigma = 0 if ell % 3 == 0 else 1
        ways *= cube_count_class(rest, sigma, 2 * ell * ell)
        if ways == 0:
            return 0
        total *= ways
    return total

def support_count(a, c):
    """Number of solutions of x a x^2 = c via the bounded-support algorithm (sum over injections f on supp a)."""
    n = len(a); S = support(a)
    total = 0
    for img in itertools.permutations(range(n), len(S)):
        f = dict(zip(S, img))
        b = list(range(n))
        for p in S:
            b[f[p]] = f[a[p]]
        b = tuple(b)
        t = comp(inv(b), c)
        total += prescribed_cube_roots(t, f)
    return total

# ------------------------------------------------------------------ brute force (numpy over Sym(N))
_PERMS = {}
def all_perms(n):
    if n not in _PERMS:
        _PERMS[n] = np.array(list(itertools.permutations(range(n))), dtype=np.int8 if n < 120 else np.int16)
    return _PERMS[n]

def brute_count(a, c, commuting=False):
    P = all_perms(len(a)); M = P.shape[0]
    rows = np.arange(M)[:, None]
    A = np.array(a); C = np.array(c)
    XX = P[rows, P]                 # x o x
    AXX = A[XX]                     # a o x o x
    W = P[rows, AXX]                # x o a o x o x
    ok = np.all(W == C, axis=1)
    if commuting:
        XA = P[rows, A[None, :].repeat(M, 0)]   # x o a
        AX = A[P]                                # a o x
        ok &= np.all(XA == AX, axis=1)
    return int(ok.sum())

def brute_cube_prescribed(t, f):
    P = all_perms(len(t)); M = P.shape[0]
    rows = np.arange(M)[:, None]
    X3 = P[rows, P[rows, P]]
    ok = np.all(X3 == np.array(t), axis=1)
    for p, q in f.items():
        ok &= (P[:, p] == q)
    return int(ok.sum())

# ------------------------------------------------------------------ verification campaigns
def verify_cent(Nmax, rng, results, quick=False):
    """For every N <= Nmax and every cycle type a, every c in C(a): compare criterion with the enumeration of C(a).
    Additionally, for N <= 7, compare with brute force over Sym(N) for a sample of c (commuting solutions)."""
    rep = {"instances": 0, "classes": 0, "mismatch": 0, "solvable": 0, "max_count": 0,
           "sym_crosscheck": 0, "sym_mismatch": 0, "per_N": {}}
    for N in range(1, Nmax + 1):
        t0 = time.time(); inst = 0; cls = 0
        for lam in partitions(N):
            a = perm_of_type(lam, N)
            cnt = cube_counter_centraliser(a)
            cls += 1
            elems = list(enumerate_centraliser(a))
            for t in elems:
                c = comp(a, t)                     # c = a t, so a^{-1}c = t in C(a)
                ex, k = cent_criterion(a, c)
                kb = cnt.get(t, 0)
                inst += 1
                if k != kb or ex != (kb > 0):
                    rep["mismatch"] += 1
                    rep.setdefault("examples", []).append({"a": a, "c": c, "crit": k, "brute": kb})
                if kb > 0: rep["solvable"] += 1
                rep["max_count"] = max(rep["max_count"], kb)
            # c not commuting with a: criterion must say no solution commuting with a
            if N <= 7:
                sample = elems if len(elems) <= 60 else rng.sample(elems, 60)
                for t in sample:
                    c = comp(a, t)
                    kb = brute_count(a, c, commuting=True)
                    rep["sym_crosscheck"] += 1
                    if kb != cent_criterion(a, c)[1]:
                        rep["sym_mismatch"] += 1
                for _ in range(10):
                    c = rand_perm(N, rng)
                    kb = brute_count(a, c, commuting=True)
                    rep["sym_crosscheck"] += 1
                    if kb != cent_criterion(a, c)[1]:
                        rep["sym_mismatch"] += 1
        rep["instances"] += inst; rep["classes"] += cls
        rep["per_N"][N] = {"classes": cls, "instances": inst, "seconds": round(time.time() - t0, 2)}
        print(f"[cent] N={N} classes={cls} instances={inst} mismatch={rep['mismatch']} t={time.time()-t0:.1f}s", flush=True)
    results["cent_exhaustive"] = rep
    return rep

def verify_cent_blocks(rng, results):
    """Blockwise families a = m cycles of length L with 3 | L, beyond N = 8 (C(a) enumerated directly)."""
    fams = [(3, 1), (3, 2), (3, 3), (3, 4), (6, 1), (6, 2), (6, 3), (9, 1), (9, 2), (9, 3), (3, 5), (12, 2), (6, 4), (3, 6)]
    rep = []
    for L, m in fams:
        N = L * m
        a = perm_of_type((L,) * m, N)
        t0 = time.time()
        cnt = cube_counter_centraliser(a)
        size = 0; mism = 0; solv = 0
        for t in enumerate_centraliser(a):
            size += 1
            c = comp(a, t)
            ex, k = cent_criterion(a, c)
            kb = cnt.get(t, 0)
            if k != kb: mism += 1
            if kb: solv += 1
        rep.append({"L": L, "m": m, "N": N, "centraliser_size": size, "cubes_in_centraliser": solv,
                    "mismatch": mism, "seconds": round(time.time() - t0, 2)})
        print(f"[cent-blocks] L={L} m={m} |C|={size} cubes={solv} mismatch={mism} t={time.time()-t0:.1f}s", flush=True)
    results["cent_blocks_3_divides_L"] = rep
    return rep

def verify_prescribed(rng, results, Nmax=8, trials=400):
    """Prescribed cube roots: random t (biased towards cycles of length divisible by 3), random partial injection."""
    rep = {"trials": 0, "mismatch": 0, "nonzero": 0, "examples": []}
    for N in range(3, Nmax + 1):
        for _ in range(trials // (Nmax - 2)):
            # t: either random, or built from a random cube root (so solvable), or with many 3-divisible cycles
            mode = rng.randrange(3)
            if mode == 0:
                t = rand_perm(N, rng)
            elif mode == 1:
                x = rand_perm(N, rng); t = comp(x, comp(x, x))
            else:
                lam = rng.choice([p for p in partitions(N) if any(L % 3 == 0 for L in p)] or list(partitions(N)))
                t = perm_of_type(lam, N)
                g = rand_perm(N, rng); t = comp(g, comp(t, inv(g)))
            s = rng.randrange(0, 5)
            S = rng.sample(range(N), min(s, N))
            if rng.random() < 0.6 and mode == 1:
                f = {p: x[p] for p in S}       # consistent with a root
            else:
                img = rng.sample(range(N), len(S)); f = dict(zip(S, img))
            ka = prescribed_cube_roots(t, f); kb = brute_cube_prescribed(t, f)
            rep["trials"] += 1
            if ka != kb:
                rep["mismatch"] += 1
                if len(rep["examples"]) < 5: rep["examples"].append({"t": t, "f": f, "alg": ka, "brute": kb})
            if kb: rep["nonzero"] += 1
    print(f"[prescribed] trials={rep['trials']} nonzero={rep['nonzero']} mismatch={rep['mismatch']}", flush=True)
    results["prescribed_cube_roots"] = rep
    return rep

def verify_support(rng, results, quick=False):
    rep = {"exhaustive": [], "random": [], "mismatch": 0, "instances": 0}
    # exhaustive: N = 5, 6 ; all a with |supp a| = s ; all c
    for N in ([5] if quick else [5, 6]):
        for s in (2, 3, 4):
            t0 = time.time(); inst = 0; mism = 0; solv = 0
            P = all_perms(N)
            As = [tuple(int(v) for v in p) for p in P if sum(1 for j in range(N) if p[j] != j) == s]
            for a in As:
                for c in P:
                    c = tuple(int(v) for v in c)
                    ka = support_count(a, c); kb = brute_count(a, c)
                    inst += 1
                    if ka != kb: mism += 1
                    if kb: solv += 1
            rep["exhaustive"].append({"N": N, "s": s, "num_a": len(As), "instances": inst, "solvable": solv,
                                      "mismatch": mism, "seconds": round(time.time() - t0, 1)})
            rep["mismatch"] += mism; rep["instances"] += inst
            print(f"[support-exh] N={N} s={s} #a={len(As)} inst={inst} solvable={solv} mismatch={mism} t={time.time()-t0:.1f}s", flush=True)
    # random: N = 7, 8, 9 ; s = 2, 3, 4 ; half planted
    plan = {7: 40, 8: 30, 9: 12} if not quick else {7: 10, 8: 4, 9: 2}
    for N, ntr in plan.items():
        for s in (2, 3, 4):
            t0 = time.time(); mism = 0; solv = 0
            for i in range(ntr):
                S = rng.sample(range(N), s)
                # random a with support exactly S: derangement of S
                while True:
                    img = S[:]; rng.shuffle(img)
                    if all(img[k] != S[k] for k in range(s)): break
                a = list(range(N))
                for k in range(s): a[S[k]] = img[k]
                a = tuple(a)
                if i % 2 == 0:
                    x = rand_perm(N, rng); c = comp(x, comp(a, comp(x, x)))
                else:
                    c = rand_perm(N, rng)
                ka = support_count(a, c); kb = brute_count(a, c)
                if ka != kb: mism += 1
                if kb: solv += 1
            rep["random"].append({"N": N, "s": s, "instances": ntr, "solvable": solv, "mismatch": mism,
                                  "seconds": round(time.time() - t0, 1)})
            rep["mismatch"] += mism; rep["instances"] += ntr
            print(f"[support-rand] N={N} s={s} inst={ntr} solvable={solv} mismatch={mism} t={time.time()-t0:.1f}s", flush=True)
    results["support"] = rep
    return rep

def save(results, path="lemmas_results.json"):
    def conv(o):
        if isinstance(o, dict): return {str(k): conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [conv(v) for v in o]
        if isinstance(o, (np.integer,)): return int(o)
        return o
    with open(path, "w") as fh:
        json.dump(conv(results), fh, indent=1)

if __name__ == "__main__" and "--extra" not in sys.argv:
    quick = "--quick" in sys.argv
    rng = random.Random(20261007)
    results = {"seed": 20261007, "quick": quick, "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    T0 = time.time()
    verify_cent(7 if quick else 8, rng, results, quick); save(results)
    verify_cent_blocks(rng, results); save(results)
    verify_prescribed(rng, results, Nmax=7 if quick else 8, trials=150 if quick else 400); save(results)
    verify_support(rng, results, quick); save(results)
    results["total_seconds"] = round(time.time() - T0, 1)
    save(results)
    print("done", results["total_seconds"], "s")

def verify_support_extra(rng, results):
    """Extra campaign: exhaustive N=7, s=2 (all a of support 2, all c); larger random grids at N=7,8,9;
    and a 3-divisible stress test (c chosen so that b^{-1}c has many cycles of length 3 or 6)."""
    rep = {"exhaustive": [], "random": [], "mismatch": 0, "instances": 0}
    N = 7; s = 2; t0 = time.time(); P = all_perms(N); inst = 0; mism = 0; solv = 0
    As = [tuple(int(v) for v in p) for p in P if sum(1 for j in range(N) if p[j] != j) == s]
    for a in As:
        for c in P:
            c = tuple(int(v) for v in c)
            ka = support_count(a, c); kb = brute_count(a, c); inst += 1
            if ka != kb: mism += 1
            if kb: solv += 1
    rep["exhaustive"].append({"N": N, "s": s, "num_a": len(As), "instances": inst, "solvable": solv, "mismatch": mism,
                              "seconds": round(time.time() - t0, 1)})
    rep["mismatch"] += mism; rep["instances"] += inst
    print(f"[support-exh2] N={N} s={s} #a={len(As)} inst={inst} solvable={solv} mismatch={mism} t={time.time()-t0:.1f}s", flush=True)
    plan = {7: 300, 8: 300, 9: 150}
    for N, ntr in plan.items():
        for s in (2, 3, 4):
            t0 = time.time(); mism = 0; solv = 0
            for i in range(ntr):
                S = rng.sample(range(N), s)
                while True:
                    img = S[:]; rng.shuffle(img)
                    if all(img[k] != S[k] for k in range(s)): break
                a = list(range(N))
                for k in range(s): a[S[k]] = img[k]
                a = tuple(a)
                mode = i % 3
                if mode == 0:
                    x = rand_perm(N, rng); c = comp(x, comp(a, comp(x, x)))
                elif mode == 1:
                    c = rand_perm(N, rng)
                else:   # stress: x of cycle type with lengths divisible by 3 where possible
                    lam = rng.choice([p for p in partitions(N) if any(L % 3 == 0 for L in p)])
                    g = rand_perm(N, rng); x = perm_of_type(lam, N); x = comp(g, comp(x, inv(g)))
                    c = comp(x, comp(a, comp(x, x)))
                    if rng.random() < 0.5:   # perturb c by a random 3-cycle so that it is usually unsolvable
                        u = rng.sample(range(N), 3); z = list(range(N)); z[u[0]] = u[1]; z[u[1]] = u[2]; z[u[2]] = u[0]
                        c = comp(tuple(z), c)
                ka = support_count(a, c); kb = brute_count(a, c)
                if ka != kb: mism += 1
                if kb: solv += 1
            rep["random"].append({"N": N, "s": s, "instances": ntr, "solvable": solv, "mismatch": mism,
                                  "seconds": round(time.time() - t0, 1)})
            rep["mismatch"] += mism; rep["instances"] += ntr
            print(f"[support-rand2] N={N} s={s} inst={ntr} solvable={solv} mismatch={mism} t={time.time()-t0:.1f}s", flush=True)
    results["support_extra"] = rep
    return rep

if __name__ == "__main__" and "--extra" in sys.argv:
    rng = random.Random(777)
    results = json.load(open("lemmas_results.json"))
    verify_support_extra(rng, results); save(results)
    print("extra done")
