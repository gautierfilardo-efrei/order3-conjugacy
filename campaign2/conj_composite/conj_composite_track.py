#!/usr/bin/env python3
"""
conj_composite_track.py -- Conj_n for composite n (n = 4, 6, 8, 9, 12, ...).

Conj_n(a, c): is there x in Sym(N) with x^n = 1 and x a x^{-1} = c ?
Conventions: image tuples, (pq)(j) = p(q(j)).

Stages (each dumps into conj_composite_results.json as soon as it finishes):
  V1  exhaustive conjugator-order profiles for all pairs (a, c) in Sym(N), N <= NMAX_EXH
  V1b disjoint-union lemma checked by brute force over Sym(N1+N2) on random pairs
  V2  affine criterion for x^n = 1 (Lemma A) vs enumeration of all conjugators on the
      blockwise family, (i) via Sym(N) for N <= 9, (ii) via the full list of L^k k! conjugators
  V3  encoding lemma (Lemma B): all subsets S with |S| | n of the exponent multiset
  V4  end-to-end reduction N3DM -> Conj_n (n = 3, 4, 6, 9) on random small instances,
      with the certificate x built as a permutation and checked (n = 3, 4, 6)

Standard library + numpy only.  Run:  python conj_composite_track.py [stages]   (default: all)
"""
import itertools, json, math, random, sys, time
from collections import Counter
import numpy as np

NMAX_EXH = 9          # exhaustive Sym(N) tables for N <= NMAX_EXH
NMAX_ORDER = 12       # Conj_n tabulated for n = 1..NMAX_ORDER
OUT = "conj_composite_results.json"
RESULTS = {"meta": {"script": "conj_composite_track.py", "NMAX_EXH": NMAX_EXH, "NMAX_ORDER": NMAX_ORDER}}

def dump():
    with open(OUT, "w") as f:
        json.dump(RESULTS, f, indent=1, default=str)

# ---------------------------------------------------------------- basic permutation helpers
def comp(p, q):
    return tuple(p[q[j]] for j in range(len(q)))

def inv(p):
    r = [0] * len(p)
    for i, v in enumerate(p):
        r[v] = i
    return tuple(r)

def cycles(p):
    seen = [False] * len(p); out = []
    for i in range(len(p)):
        if not seen[i]:
            cyc = []; j = i
            while not seen[j]:
                seen[j] = True; cyc.append(j); j = p[j]
            out.append(cyc)
    return out

def ctype(p):
    return tuple(sorted((len(c) for c in cycles(p)), reverse=True))

def order(p):
    o = 1
    for c in cycles(p):
        o = o * len(c) // math.gcd(o, len(c))
    return o

def cyc_str(p):
    s = "".join("(" + " ".join(str(j + 1) for j in c) + ")" for c in cycles(p) if len(c) > 1)
    return s or "()"

def perm_of_ctype(lam):
    p = []; base = 0
    for l in lam:
        p += [base + (i + 1) % l for i in range(l)]; base += l
    return tuple(p)

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

# ---------------------------------------------------------------- numpy machinery over all of Sym(N)
_PERMS = {}
def all_perms(N):
    if N not in _PERMS:
        P = np.array(list(itertools.permutations(range(N))), dtype=np.int8)
        _PERMS[N] = (P, perm_orders(P))
    return _PERMS[N]

def perm_orders(P):
    R, N = P.shape
    rows = np.arange(R)
    ordv = np.ones(R, dtype=np.int64)
    for j in range(N):
        cur = P[:, j].astype(np.int64); cnt = np.ones(R, dtype=np.int64)
        for _ in range(N - 1):
            act = cur != j
            if not act.any(): break
            cnt += act
            cur = np.where(act, P[rows, cur], cur)
        ordv = np.lcm(ordv, cnt)
    return ordv

def conj_all(P, a):
    """rows x -> x a x^{-1}  (as image tuples):  c[x[j]] = x[a[j]]."""
    a = np.asarray(a); R, N = P.shape
    C = np.empty_like(P)
    C[np.arange(R)[:, None], P.astype(np.int64)] = P[:, a]
    return C

def keys_of(C):
    N = C.shape[1]
    return (C.astype(np.int64) * (N ** np.arange(N, dtype=np.int64))).sum(axis=1)

def yes_set_from_orders(orders_present, nmax=NMAX_ORDER):
    """{n <= nmax : some conjugator has order dividing n}."""
    return sorted(n for n in range(1, nmax + 1) if any(n % o == 0 for o in orders_present))

def conj_orders_bruteforce(a, c):
    """Set of orders of all x in Sym(N) with x a x^{-1} = c (N <= 9)."""
    N = len(a); P, ordv = all_perms(N)
    mask = (conj_all(P, a) == np.asarray(c, dtype=np.int8)).all(axis=1)
    return sorted(set(ordv[mask].tolist()))

# ================================================================ V1: exhaustive profiles N <= 9
def stage_V1():
    t0 = time.time(); res = {}
    for N in range(1, NMAX_EXH + 1):
        P, ordv = all_perms(N); R = len(P); rows = np.arange(R)
        maxo = int(ordv.max())
        per_n_pairs = Counter(); strict = {}; yesset_dist = Counter(); total_pairs = 0
        for lam in partitions(N):
            a = perm_of_ctype(lam)
            C = conj_all(P, a); k = keys_of(C)
            uniq, invidx = np.unique(k, return_inverse=True)
            G = len(uniq); cls = G                       # |class(a)| = number of distinct conjugates
            total_pairs += cls * cls
            M = np.zeros((G, maxo + 1), dtype=bool)
            M[invidx, ordv] = True
            yes = {}
            for n in range(1, NMAX_ORDER + 1):
                yes[n] = M[:, divisors(n)].any(axis=1) if max(divisors(n)) <= maxo else M[:, [d for d in divisors(n) if d <= maxo]].any(axis=1)
                per_n_pairs[n] += int(yes[n].sum()) * cls
            # distribution of yes-sets (as tuple) weighted by cls
            packed = np.zeros(G, dtype=np.int64)
            for n in range(1, NMAX_ORDER + 1):
                packed |= yes[n].astype(np.int64) << n
            for v, cnt in zip(*np.unique(packed, return_counts=True)):
                ys = tuple(n for n in range(1, NMAX_ORDER + 1) if (int(v) >> n) & 1)
                yesset_dist[ys] += int(cnt) * cls
            # strict witnesses: Conj_n yes but Conj_d no for every proper divisor d
            for n in range(2, NMAX_ORDER + 1):
                w = yes[n].copy()
                for d in divisors(n)[:-1]:
                    w &= ~yes[d]
                if w.any():
                    g = int(np.flatnonzero(w)[0]); row = int(np.flatnonzero(invidx == g)[0])
                    c = tuple(int(v) for v in C[row])
                    ords = sorted(int(o) for o in np.flatnonzero(M[g]))
                    strict.setdefault(n, {"count_pairs": 0, "example": None})
                    strict[n]["count_pairs"] += int(w.sum()) * cls
                    if strict[n]["example"] is None:
                        strict[n]["example"] = {"a": cyc_str(a), "c": cyc_str(c), "conjugator_orders": ords}
        res[N] = {
            "same_class_pairs": total_pairs,
            "conj_n_yes_pairs": {n: per_n_pairs[n] for n in range(1, NMAX_ORDER + 1)},
            "conj_n_yes_fraction": {n: per_n_pairs[n] / total_pairs for n in range(1, NMAX_ORDER + 1)},
            "strict_witnesses": strict,
            "yes_set_distribution": {",".join(map(str, k)) if k else "none": v for k, v in sorted(yesset_dist.items())},
        }
        print(f"V1 N={N}: pairs={total_pairs}  yes4={per_n_pairs[4]} yes2={per_n_pairs[2]}  strict4={strict.get(4,{}).get('count_pairs',0)} strict6={strict.get(6,{}).get('count_pairs',0)}  {time.time()-t0:.1f}s", flush=True)
    RESULTS["V1_exhaustive_profiles"] = res; dump()

# ================================================================ V1b: disjoint-union lemma
def random_pair(N, rng):
    a = list(range(N)); rng.shuffle(a); a = tuple(a)
    x = list(range(N)); rng.shuffle(x); x = tuple(x)
    c = comp(comp(x, a), inv(x))
    return a, c

def stage_V1b(trials=150, seed=1):
    rng = random.Random(seed); t0 = time.time(); checked = 0; fails = []
    while checked < trials:
        N1 = rng.randint(1, 7); N2 = rng.randint(1, 8 - N1 + 1)
        if N1 + N2 > 9: continue
        a1, c1 = random_pair(N1, rng); a2, c2 = random_pair(N2, rng)
        if set(ctype(a1)) & set(ctype(a2)): continue        # cycle lengths (incl. 1) must be disjoint
        a = a1 + tuple(N1 + v for v in a2); c = c1 + tuple(N1 + v for v in c2)
        Y = yes_set_from_orders(conj_orders_bruteforce(a, c))
        Y1 = yes_set_from_orders(conj_orders_bruteforce(a1, c1))
        Y2 = yes_set_from_orders(conj_orders_bruteforce(a2, c2))
        pred = sorted(set(Y1) & set(Y2))
        checked += 1
        if Y != pred:
            fails.append({"a1": cyc_str(a1), "c1": cyc_str(c1), "a2": cyc_str(a2), "c2": cyc_str(c2), "Y": Y, "pred": pred})
    # the counterexamples to the naive "order-lowering" idea: Conj_4 yes, Conj_2 no (5-cycle, square)
    a = perm_of_ctype((5,)); c = comp(a, a)
    o5 = conj_orders_bruteforce(a, c)
    a7 = perm_of_ctype((7,)); c7 = tuple(a7[a7[a7[j]]] for j in range(7))      # a^3, slope 3 of order 6 mod 7
    o7 = conj_orders_bruteforce(a7, c7)
    RESULTS["V1b_disjoint_union_lemma"] = {"trials": checked, "failures": fails,
        "five_cycle_square_conjugator_orders": o5, "seven_cycle_cube_conjugator_orders": o7,
        "time_s": round(time.time() - t0, 1)}
    print(f"V1b: {checked} unions checked, failures={len(fails)}; orders((12345)->a^2)={o5}; orders((1..7)->a^3)={o7}", flush=True)
    dump()

# ================================================================ Lemma A: affine criterion for x^n = 1
def blockwise(L, k, m):
    a = [0] * (L * k); c = [0] * (L * k)
    for i in range(k):
        for z in range(L):
            a[i * L + z] = i * L + (z + 1) % L
            c[i * L + z] = i * L + (z + m[i]) % L
    return tuple(a), tuple(c)

def criterion_A(L, k, m, n):
    """exists sigma in Sym(k), sigma^n = 1, with (prod_{i in C} m_i)^{n/|C|} = 1 mod L on every cycle C."""
    for sig in itertools.permutations(range(k)):
        ok = True
        for C in cycles(sig):
            l = len(C)
            if n % l: ok = False; break
            R = 1
            for i in C: R = R * m[i] % L
            if pow(R, n // l, L) != 1 % L: ok = False; break
        if ok: return True
    return False

def all_conjugators_blockwise(L, k, m):
    """All L^k k! maps x: B_i -> B_{sigma(i)}, z -> m_{sigma(i)} z + t_i, as an int array (rows, N)."""
    N = L * k; sigs = list(itertools.permutations(range(k))); ts = list(itertools.product(range(L), repeat=k))
    X = np.empty((len(sigs) * len(ts), N), dtype=np.int64)
    z = np.arange(L)
    r = 0
    for sig in sigs:
        lin = np.empty(N, dtype=np.int64); toff = np.empty(N, dtype=np.int64)
        for i in range(k):
            lin[i * L:(i + 1) * L] = (m[sig[i]] * z) % L
            toff[i * L:(i + 1) * L] = sig[i] * L
        T = np.array(ts, dtype=np.int64)                    # (L^k, k)
        blk = np.repeat(T, L, axis=1)                        # (L^k, N): t_i on block i
        X[r:r + len(ts)] = toff[None, :] + (lin[None, :] + blk) % L
        r += len(ts)
    return X

def stage_V2(seed=2):
    rng = random.Random(seed); t0 = time.time(); out = {"sym_N_bruteforce": [], "full_conjugator_list": []}
    # (i) genuinely independent: brute force over Sym(N), N <= 9, all unit vectors m
    for (L, k) in [(3, 1), (4, 1), (5, 1), (6, 1), (7, 1), (8, 1), (9, 1), (3, 2), (4, 2), (3, 3)]:
        units = [u for u in range(1, L) if math.gcd(u, L) == 1]
        cnt = 0; bad = 0
        for m in itertools.product(units, repeat=k):
            a, c = blockwise(L, k, m)
            Y = yes_set_from_orders(conj_orders_bruteforce(a, c))
            Yc = [n for n in range(1, NMAX_ORDER + 1) if criterion_A(L, k, m, n)]
            cnt += 1; bad += (Y != Yc)
        out["sym_N_bruteforce"].append({"L": L, "k": k, "instances": cnt, "mismatches": bad})
        print(f"V2(i) L={L} k={k}: {cnt} instances, mismatches={bad}", flush=True)
    # (ii) full list of conjugators (block lemma, proved in the paper), actual powers of permutations
    grid = [(5, 3), (7, 3), (11, 3), (13, 3), (5, 4), (7, 4), (11, 4), (5, 5), (7, 5)]
    for (L, k) in grid:
        units = [u for u in range(1, L) if math.gcd(u, L) == 1]
        ninst = 12 if L ** k * math.factorial(k) < 3e5 else 6
        bad = 0
        for _ in range(ninst):
            m = [rng.choice(units) for _ in range(k)]
            if rng.random() < 0.3: m[rng.randrange(k)] = 1
            a, c = blockwise(L, k, m); N = L * k
            X = all_conjugators_blockwise(L, k, m)
            assert (conj_all(X, a) == np.asarray(c)).all(), "not a conjugator"
            rows = np.arange(len(X))[:, None]; pw = X.copy(); ident = np.arange(N)
            Y = []
            for n in range(1, NMAX_ORDER + 1):
                if n > 1: pw = X[rows, pw]
                if (pw == ident).all(axis=1).any(): Y.append(n)
            Yc = [n for n in range(1, NMAX_ORDER + 1) if criterion_A(L, k, m, n)]
            bad += (Y != Yc)
        out["full_conjugator_list"].append({"L": L, "k": k, "instances": ninst, "conjugators_each": L ** k * math.factorial(k), "mismatches": bad})
        print(f"V2(ii) L={L} k={k}: {ninst} instances x {L**k*math.factorial(k)} conjugators, mismatches={bad}  {time.time()-t0:.0f}s", flush=True)
    out["time_s"] = round(time.time() - t0, 1)
    RESULTS["V2_affine_criterion_lemmaA"] = out; dump()

# ================================================================ Lemma B: the exponent encoding
def encode(n, classes, Bp, T=None):
    """classes: list of n lists of sizes s(u) in [1, Bp].  Returns (T, exponents list of (class j, s, e))."""
    M = n * Bp + 1
    w = [(n + 1) ** j for j in range(n - 1)]          # w_1..w_{n-1}
    W = sum(w)
    T0 = M * (n + 1) ** (n - 1)
    if T is None: T = T0
    E = []
    for j, cl in enumerate(classes):
        for s in cl:
            e = s + M * w[j] if j < n - 1 else s - M * W - Bp
            E.append((j, s, e))
    return T, T0, E

def lemma_B_check(n, classes, Bp, T=None):
    T, T0, E = encode(n, classes, Bp, T)
    k = len(E); bad = 0; checked = 0
    for l in divisors(n):
        for S in itertools.combinations(range(k), l):
            se = sum(E[i][2] for i in S)
            cond = ((n // l) * se) % T == 0
            good = (l == n and len({E[i][0] for i in S}) == n and sum(E[i][1] for i in S) == Bp)
            checked += 1; bad += (cond != good)
    return checked, bad, T

def random_n3dm(n, q, Bp, rng, planted):
    """n classes of q sizes in [1,Bp]; classes 4..n are padding (size 1) as in the paper."""
    if planted:
        W, X, Y = [], [], []
        for _ in range(q):
            while True:
                w = rng.randint(1, Bp - 2 - (n - 3)); x = rng.randint(1, Bp - 1 - w - (n - 3)); y = Bp - (n - 3) - w - x
                if 1 <= y <= Bp: break
            W.append(w); X.append(x); Y.append(y)
    else:
        W = [rng.randint(1, Bp) for _ in range(q)]; X = [rng.randint(1, Bp) for _ in range(q)]; Y = [rng.randint(1, Bp) for _ in range(q)]
    classes = [W, X, Y] + [[1] * q for _ in range(n - 3)]
    return classes

def n3dm_bruteforce(classes, Bp):
    q = len(classes[0]); n = len(classes)
    for perms in itertools.product(*[itertools.permutations(range(q)) for _ in range(n - 1)]):
        if all(classes[0][i] + sum(classes[j + 1][perms[j][i]] for j in range(n - 1)) == Bp for i in range(q)):
            return True
    return False

def stage_V3(seed=3):
    rng = random.Random(seed); t0 = time.time(); out = []
    for n, q, Bp_max, ninst in [(3, 3, 9, 40), (4, 3, 8, 30), (5, 2, 8, 20), (6, 2, 8, 20), (8, 2, 7, 8), (9, 2, 7, 6)]:
        tot = 0; bad = 0
        for i in range(ninst):
            Bp = rng.randint(max(3, n), Bp_max + n); classes = random_n3dm(n, q, Bp, rng, planted=(i % 2 == 0))
            # T = T0 and also a larger T (prime-minus-one style), the lemma needs only T >= T0
            for T in (None, 2 * encode(n, classes, Bp)[1] + 1):
                ch, b, T_used = lemma_B_check(n, classes, Bp, T); tot += ch; bad += b
        out.append({"n": n, "q": q, "instances": ninst, "subsets_checked": tot, "mismatches": bad})
        print(f"V3 n={n} q={q}: {tot} subsets checked over {ninst} instances, mismatches={bad}  {time.time()-t0:.0f}s", flush=True)
    RESULTS["V3_encoding_lemmaB"] = {"grid": out, "time_s": round(time.time() - t0, 1)}; dump()

# ================================================================ V4: end-to-end reduction
def least_prime_ge(x):
    def isprime(v):
        if v < 2: return False
        if v % 2 == 0: return v == 2
        f = 3
        while f * f <= v:
            if v % f == 0: return False
            f += 2
        return True
    while not isprime(x): x += 1
    return x

def primitive_root(L):
    phi = L - 1; fs = set()
    v = phi; f = 2
    while f * f <= v:
        while v % f == 0: fs.add(f); v //= f
        f += 1
    if v > 1: fs.add(v)
    for g in range(2, L):
        if all(pow(g, phi // p, L) != 1 for p in fs): return g
    raise ValueError

def exact_cover_blocks(n, m, L):
    """Lemma A as an exact-cover search: partition blocks into sets S, |S| | n, (prod m)^{n/|S|} = 1 mod L."""
    k = len(m)
    adm = {}
    for l in divisors(n):
        for S in itertools.combinations(range(k), l):
            R = 1
            for i in S: R = R * m[i] % L
            if pow(R, n // l, L) == 1 % L: adm.setdefault(S[0], []).append(S)
    # backtracking: cover the smallest free block
    free = set(range(k))
    def rec():
        if not free: return []
        i = min(free)
        for S in adm.get(i, []):
            if all(j in free for j in S):
                for j in S: free.discard(j)
                r = rec()
                for j in S: free.add(j)
                if r is not None: return [S] + r
        return None
    return rec()

def stage_V4(seed=4):
    rng = random.Random(seed); t0 = time.time(); out = []
    for n, q, Bp_max, ninst, build in [(3, 3, 7, 24, True), (4, 2, 6, 24, True), (4, 3, 5, 10, False), (6, 2, 5, 16, True), (9, 2, 4, 6, False)]:
        bad = 0; yes_cnt = 0; Nmax = 0
        for i in range(ninst):
            Bp = rng.randint(n, Bp_max + n); classes = random_n3dm(n, q, Bp, rng, planted=(i % 2 == 0))
            truth = n3dm_bruteforce(classes, Bp)
            T, T0, E = encode(n, classes, Bp)
            L = least_prime_ge(T0 + 1); Tp = L - 1          # prime choice: T = L - 1 >= T0
            g = primitive_root(L)
            m = [pow(g, e % Tp, L) for (_, _, e) in E]
            assert all(mi != 1 for mi in m)
            k = len(m); N = k * L; Nmax = max(Nmax, N)
            sol = exact_cover_blocks(n, m, L)
            got = sol is not None
            if got:
                yes_cnt += 1
                assert all(len(S) == n for S in sol)
                if build:
                    # certificate as a permutation: sigma cycles the blocks of each n-set, t = 0
                    a = np.empty(N, dtype=np.int64); c = np.empty(N, dtype=np.int64); x = np.empty(N, dtype=np.int64)
                    z = np.arange(L)
                    for b in range(k):
                        a[b * L:(b + 1) * L] = b * L + (z + 1) % L
                        c[b * L:(b + 1) * L] = b * L + (z + m[b]) % L
                    for S in sol:
                        for idx, b in enumerate(S):
                            tb = S[(idx + 1) % n]
                            x[b * L:(b + 1) * L] = tb * L + (m[tb] * z) % L
                    assert sorted(x.tolist()) == list(range(N))
                    xinv = np.empty(N, dtype=np.int64); xinv[x] = np.arange(N)
                    assert (x[a[xinv]] == c).all(), "x a x^-1 != c"
                    pw = x.copy()
                    for _ in range(n - 1): pw = x[pw]
                    assert (pw == np.arange(N)).all(), "x^n != 1"
                    # x^d != 1 for proper divisors d (order exactly n)
                    for d in divisors(n)[:-1]:
                        pd = x.copy()
                        for _ in range(d - 1): pd = x[pd]
                        assert not (pd == np.arange(N)).all()
            bad += (got != truth)
        out.append({"n": n, "q": q, "instances": ninst, "yes_instances": yes_cnt, "mismatches": bad, "max_N": Nmax, "certificate_built_and_checked": build})
        print(f"V4 n={n} q={q}: {ninst} instances ({yes_cnt} yes), mismatches={bad}, max N={Nmax}  {time.time()-t0:.0f}s", flush=True)
    RESULTS["V4_end_to_end_reduction"] = {"grid": out, "time_s": round(time.time() - t0, 1)}; dump()

if __name__ == "__main__":
    stages = sys.argv[1:] or ["V1", "V1b", "V2", "V3", "V4"]
    try:
        with open(OUT) as f: RESULTS.update(json.load(f))
    except Exception: pass
    T0 = time.time()
    for s in stages:
        globals()["stage_" + s]()
    RESULTS["meta"]["total_time_s"] = round(time.time() - T0, 1); dump()
    print("done", round(time.time() - T0, 1), "s")
