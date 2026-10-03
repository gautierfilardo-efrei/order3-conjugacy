"""algo_track.py -- "Algorithme" track for P3:  x a x^2 = c  in S_N.

Conventions: permutations are image tuples / rows, product right to left,
(pq)(j) = p(q(j)).  The word  x a x x  applied to j is  x(a(x(x(j)))).

Part A  exhaustive tabulation (numpy, all of S_N for N <= NMAX):
        for every cycle-type representative a and every c, the number of
        solutions x and the set of cycle types lambda = ctype(x) realised.
Part B  invariant tests: does a polynomial-time invariant of the pair (a,c)
        (cycle types of short words in a, c; orbit structure of <a,c>) decide
        solvability?  Same test per fixed lambda (the p(N)-enumeration route).
Part C  special families of a (N-cycle, distinct cycle lengths, coprime
        cycle lengths): same questions restricted to those a.
Part D  propagation-solver cost, a = N-cycle versus random a (sanity check).

Standard library + numpy.  Results are written to algo_results.json.
"""
from __future__ import annotations
import itertools, json, math, sys, time, hashlib, random
from collections import Counter, defaultdict
import numpy as np

# ----------------------------------------------------------------- basics
def compose(p, q):
    return tuple(p[q[j]] for j in range(len(p)))

def inverse(p):
    o = [0] * len(p)
    for i, v in enumerate(p): o[v] = i
    return tuple(o)

def mul(*ps):
    r = ps[0]
    for p in ps[1:]: r = compose(r, p)
    return r

def cycle_type(p):
    seen = [False] * len(p); c = []
    for i in range(len(p)):
        if not seen[i]:
            j = i; L = 0
            while not seen[j]: seen[j] = True; j = p[j]; L += 1
            c.append(L)
    return tuple(sorted(c, reverse=True))

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

def perm_of_type(lam):
    """canonical representative of cycle type lam (cycles on consecutive points)"""
    N = sum(lam); p = list(range(N)); s = 0
    for L in lam:
        for i in range(L): p[s + i] = s + (i + 1) % L
        s += L
    return tuple(p)

def z_lambda(lam):
    """centraliser order of cycle type lam"""
    z = 1
    for L, m in Counter(lam).items(): z *= (L ** m) * math.factorial(m)
    return z

def all_perms(N):
    """all permutations of range(N) as an int8 array of shape (N!, N)"""
    P = np.zeros((1, 1), dtype=np.int8)
    for n in range(2, N + 1):
        P = np.concatenate([np.insert(P, k, n - 1, axis=1) for k in range(n)], axis=0)
    return P

def compose_rows(P, Q):
    """row-wise composition (P o Q)[i] = P[i] o Q[i]  :  j -> P[i][Q[i][j]]"""
    return np.take_along_axis(P, Q, axis=1)

def inverse_rows(P):
    return np.argsort(P, axis=1).astype(np.int8)

def ctype_codes(P, N):
    """per-row cycle-type code: sum_d m_d (N+1)^d  (m_d = number of d-cycles)"""
    M = P.shape[0]
    idx = np.arange(N, dtype=np.int8)
    f = np.zeros((N + 1, M), dtype=np.int64)  # f[k] = fixed points of P^k
    Pk = P
    for k in range(1, N + 1):
        f[k] = (Pk == idx).sum(axis=1)
        if k < N: Pk = np.take_along_axis(P, Pk, axis=1)
    m = np.zeros((N + 1, M), dtype=np.int64)
    for k in range(1, N + 1):
        s = np.zeros(M, dtype=np.int64)
        for d in range(1, k):
            if k % d == 0: s += d * m[d]
        m[k] = (f[k] - s) // k
    code = np.zeros(M, dtype=np.int64)
    for d in range(1, N + 1): code += m[d] * (N + 1) ** d
    return code

def code_of_type(lam, N):
    c = 0
    for L, mL in Counter(lam).items(): c += mL * (N + 1) ** L
    return c

def type_of_code(code, N):
    lam = []
    for d in range(N, 0, -1):
        md = code // (N + 1) ** d; code -= md * (N + 1) ** d
        lam += [d] * int(md)
    return tuple(lam)

def rank_rows(P, N):
    w = (N ** np.arange(N, dtype=np.int64))
    return P.astype(np.int64) @ w

def orbit_codes(P, a, N):
    """per row c=P[i]: multiset of orbit sizes of <a, c>, encoded as a partition code"""
    M = P.shape[0]
    a = np.asarray(a, dtype=np.int8); ainv = np.asarray(inverse(tuple(int(v) for v in a)), dtype=np.int8)
    Pinv = inverse_rows(P)
    lab = np.tile(np.arange(N, dtype=np.int8), (M, 1))
    for _ in range(N):
        lab = np.minimum.reduce([lab, lab[:, a], lab[:, ainv],
                                 np.take_along_axis(lab, P, 1), np.take_along_axis(lab, Pinv, 1)])
    sizes = np.zeros((M, N), dtype=np.int64)
    for l in range(N):
        sizes[:, l] = (lab == l).sum(axis=1)
    sizes.sort(axis=1)
    code = np.zeros(M, dtype=np.int64)
    for l in range(N):  # multiset code: count of each size s in base (N+1)
        pass
    # encode multiset of sizes: sum over sizes s>0 of (N+1)^s * multiplicity
    for s in range(1, N + 1):
        code += (sizes == s).sum(axis=1) * (N + 1) ** s
    return code

# ----------------------------------------------------------------- Part A/B core
WORDS_SHORT = {  # words in a (A = a^-1), c (C = c^-1), read as product left-to-right for the label
    'ac': ('a', 'c'), 'Ac': ('A', 'c'),
}
def word_rows(word, P, Pinv, a, ainv):
    """rows of the word product; letters applied right to left: (w)(j)"""
    N = P.shape[1]; M = P.shape[0]
    cur = np.tile(np.arange(N, dtype=np.int8), (M, 1))
    for letter in reversed(word):
        if letter == 'a': cur = np.asarray(a, dtype=np.int8)[cur]
        elif letter == 'A': cur = np.asarray(ainv, dtype=np.int8)[cur]
        elif letter == 'c': cur = np.take_along_axis(P, cur, 1)
        elif letter == 'C': cur = np.take_along_axis(Pinv, cur, 1)
    return cur

def word_list(maxlen):
    """words of length <= maxlen in a,A,c,C containing at least one c/C, modulo obvious redundancy
    (we keep all -- cheap at N<=8)"""
    ws = []
    for L in range(1, maxlen + 1):
        for w in itertools.product('aAcC', repeat=L):
            if any(ch in 'cC' for ch in w): ws.append(''.join(w))
    return ws

def combine_codes(cols):
    """combine several int64 label columns into one group id (via np.unique on rows)"""
    A = np.stack(cols, axis=1)
    _, inv = np.unique(A, axis=0, return_inverse=True)
    return inv.reshape(-1)

def decides(group_id, member):
    """True iff `member` (bool per row) is constant on each group; else returns a witness pair
    of rows (i, j) with the same group but different membership, chosen with minimal group id."""
    key = group_id.astype(np.int64) * 2 + member.astype(np.int64)
    uk = np.unique(key)
    g = uk // 2
    dup = g[1:] == g[:-1]
    if not dup.any(): return True, None
    bad_g = g[:-1][dup][0]
    i = int(np.flatnonzero((group_id == bad_g) & ~member)[0])
    j = int(np.flatnonzero((group_id == bad_g) & member)[0])
    return False, (i, j)

def tabulate_N(N, maxword=3, do_orbits=True, per_lambda=True, verbose=True):
    t0 = time.time()
    P = all_perms(N); M = P.shape[0]
    Pinv = inverse_rows(P)
    Pcode = rank_rows(P, N); order = np.argsort(Pcode); sorted_codes = Pcode[order]
    ct = ctype_codes(P, N)                       # ctype code of every row (as x or as c)
    parts = list(partitions(N)); part_code = {lam: code_of_type(lam, N) for lam in parts}
    code_part = {v: k for k, v in part_code.items()}
    x2 = compose_rows(P, P)
    words = [w for w in word_list(maxword) if w not in ('c', 'C')]
    out = {'N': N, 'n_perms': int(M), 'per_a': {}, 'invariant_tests': {}, 'lambda_tests': {}}
    inv_fail = {}      # invariant name -> smallest witness (over a)
    lam_fail = {}
    total_solvable_pairs = 0; total_pairs = 0
    nsol_hist_total = Counter()
    for lam_a in parts:
        a = perm_of_type(lam_a); ainv = inverse(a)
        xa = P[:, np.asarray(a)]                 # (x o a)[i] : j -> x[a[j]]
        c = np.take_along_axis(xa, x2, 1)        # (x a x x)
        ccode = rank_rows(c, N)
        crow = order[np.searchsorted(sorted_codes, ccode)]   # row index of c in P
        nsol = np.bincount(crow, minlength=M)
        solvable = nsol > 0
        cls_size = math.factorial(N) // z_lambda(lam_a)
        total_solvable_pairs += int(solvable.sum()) * cls_size; total_pairs += M * cls_size
        hist = Counter(nsol.tolist())
        for k, v in hist.items(): nsol_hist_total[k] += v * cls_size
        # solvability by ctype(c)
        by_ct = defaultdict(lambda: [0, 0, Counter()])
        for code in np.unique(ct):
            mask = ct == code
            lam_c = code_part[int(code)]
            by_ct[lam_c] = [int(mask.sum()), int((solvable & mask).sum()), Counter(nsol[mask].tolist())]
        rec = {'ctype_a': lam_a, 'n_c': int(M), 'n_solvable_c': int(solvable.sum()),
               'nsol_hist': {str(k): int(v) for k, v in sorted(hist.items())},
               'by_ctype_c': {str(k): {'n': v[0], 'solvable': v[1], 'nsol_hist': {str(kk): int(vv) for kk, vv in sorted(v[2].items())}}
                              for k, v in by_ct.items()}}
        # ---- invariants of the pair (a,c): cycle types of words
        wct = {}
        for w in words:
            wct[w] = ctype_codes(word_rows(w, P, Pinv, a, ainv), N)
        invariants = {
            'I1:ctype(c)': [ct],
            'I2:+ctype(ac)': [ct, wct['ac']],
            'I3:+ctype(Ac)': [ct, wct['ac'], wct['Ac']],
            'I4:all words len<=%d' % maxword: [ct] + [wct[w] for w in words],
        }
        if do_orbits:
            orb = orbit_codes(P, a, N)
            invariants['I5:I4+orbits<a,c>'] = invariants['I4:all words len<=%d' % maxword] + [orb]
            invariants['I0:orbits<a,c> only'] = [orb]
        rec['invariants'] = {}
        gids = {name: combine_codes(cols) for name, cols in invariants.items()}
        for name, gid in gids.items():
            ok, wit = decides(gid, solvable)
            entry = {'decides': bool(ok), 'n_groups': int(gid.max() + 1)}
            if not ok:
                i, j = wit
                entry['witness'] = {'a': list(a), 'c_unsolvable': [int(v) for v in P[i]], 'c_solvable': [int(v) for v in P[j]],
                                    'nsol': [int(nsol[i]), int(nsol[j])],
                                    'ctype_c': list(code_part[int(ct[i])])}
                if name not in inv_fail: inv_fail[name] = dict(entry['witness'], N=N)
            rec['invariants'][name] = entry
        # ---- per lambda = ctype(x)
        if per_lambda:
            lam_rows = {}
            for code in np.unique(ct):
                xs = ct == code                      # rows x of type lambda
                sol_c = np.zeros(M, dtype=bool); sol_c[crow[xs]] = True
                lam = code_part[int(code)]
                cnt_c = np.bincount(crow[xs], minlength=M)
                e = {'n_x': int(xs.sum()), 'n_solvable_c': int(sol_c.sum()),
                     'nsol_hist': {str(k): int(v) for k, v in sorted(Counter(cnt_c.tolist()).items())}, 'invariants': {}}
                for name in ('I1:ctype(c)', 'I3:+ctype(Ac)', 'I4:all words len<=%d' % maxword) + (('I5:I4+orbits<a,c>',) if do_orbits else ()):
                    gid = gids[name]
                    ok, wit = decides(gid, sol_c)
                    e['invariants'][name] = bool(ok)
                    if not ok:
                        key = (name, str(lam))
                        if key not in lam_fail:
                            i, j = wit
                            lam_fail[key] = {'N': N, 'lambda': list(lam), 'a': list(a), 'c_unsolvable': [int(v) for v in P[i]],
                                             'c_solvable': [int(v) for v in P[j]], 'ctype_c': list(code_part[int(ct[i])])}
                lam_rows[str(lam)] = e
            rec['per_lambda'] = lam_rows
            # number of distinct lambda realised per solvable c
            pairs = np.unique(np.stack([crow, ct], 1), axis=0)
            nlam = np.bincount(pairs[:, 0], minlength=M)
            rec['n_lambda_per_c_hist'] = {str(k): int(v) for k, v in sorted(Counter(nlam[solvable].tolist()).items())}
        out['per_a'][str(lam_a)] = rec
        if verbose: print(f'N={N} a={lam_a} solvable {int(solvable.sum())}/{M}  t={time.time()-t0:.1f}s', flush=True)
    out['fraction_solvable_pairs'] = total_solvable_pairs / total_pairs
    out['nsol_hist_all_pairs'] = {str(k): int(v) for k, v in sorted(nsol_hist_total.items())}
    out['first_failure_per_invariant'] = inv_fail
    out['first_failure_per_lambda'] = {f'{k[0]} | lambda={k[1]}': v for k, v in lam_fail.items()}
    out['time_s'] = time.time() - t0
    return out

if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    res = {}
    for N in range(3, nmax + 1):
        res[str(N)] = tabulate_N(N, maxword=3 if N >= 9 else 4, do_orbits=True)
    json.dump(res, open('algo_results_partA.json', 'w'), indent=1, default=str)
