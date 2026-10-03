"""Constraint-propagation solver for one-variable word equations  w(x) = c  in S_N.

The word is a list of letters: 'x' (the unknown), 'X' (its inverse) or a constant
permutation (tuple).  The product is read right to left, (pq)(j) = p(q(j)), so
the word [x, a, x, x] denotes  x a x x  and  (x a x x)(j) = x(a(x(x(j)))).

For every point j the equation gives a chain  j -> ... -> c(j) through the
letters read right to left.  Propagation: when all letters on a chain are
determined except one application of x (or X), that application is forced.
Branching: on an undetermined image x(j), smallest domain first.
Standard library only.
"""
from __future__ import annotations
import itertools, random, time

NA = -1

def compose(p, q): return tuple(p[q[j]] for j in range(len(p)))
def inverse(p):
    o = [0] * len(p)
    for i, v in enumerate(p): o[v] = i
    return tuple(o)
def mul(*ps):
    r = ps[0]
    for p in ps[1:]: r = compose(r, p)
    return r
def parity(p):
    s = 0; seen = [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j = i; L = 0
            while not seen[j]: seen[j] = True; j = p[j]; L += 1
            s += L - 1
    return s % 2
def cycle_type(p):
    seen = [False] * len(p); c = []
    for i in range(len(p)):
        if not seen[i]:
            j = i; L = 0
            while not seen[j]: seen[j] = True; j = p[j]; L += 1
            c.append(L)
    return tuple(sorted(c, reverse=True))
def random_perm(N, rng):
    p = list(range(N)); rng.shuffle(p); return tuple(p)
def random_even(N, rng):
    while True:
        p = random_perm(N, rng)
        if parity(p) == 0: return p
def evaluate(word, x):
    xi = inverse(x)
    r = tuple(range(len(x)))
    for L in reversed(word):
        p = x if L == 'x' else xi if L == 'X' else L
        r = compose(p, r)
    return r


def solve_word(word, c, max_solutions=1, node_budget=10**6, order=None):
    """Enumerate x in S_N with w(x) = c.  Returns (solutions, stats)."""
    N = len(c)
    letters = list(reversed(word))          # application order
    sols = []
    st = {"nodes": 0, "exhausted": False, "decisions": 0, "forced": 0}

    def assign(x, xi, j, v):
        if x[j] != NA: return x[j] == v
        if xi[v] != NA: return False
        x[j] = v; xi[v] = j
        return True

    def propagate(x, xi):
        changed = True
        while changed:
            changed = False
            for j in range(N):
                # forward walk from j: track known prefix; then backward walk from c(j)
                # positions: value after applying letters[0..k-1]
                vals = [j]; ok = True
                for L in letters:
                    v = vals[-1]
                    if L == 'x': nv = x[v]
                    elif L == 'X': nv = xi[v]
                    else: nv = L[v]
                    if nv == NA: ok = False; break
                    vals.append(nv)
                if ok:
                    if vals[-1] != c[j]: return False
                    continue
                k = len(vals) - 1               # first unknown letter index
                # backward from c(j) through letters[k+1..]
                back = [c[j]]; okb = True
                for L in reversed(letters[k + 1:]):
                    v = back[-1]
                    if L == 'x': nv = xi[v]
                    elif L == 'X': nv = x[v]
                    else: nv = inverse_cache[id(L)][v]
                    if nv == NA: okb = False; break
                    back.append(nv)
                if not okb: continue
                # letters[k] applied to vals[k] must give back[-1]
                L = letters[k]; src = vals[k]; dst = back[-1]
                if L == 'x':
                    if not assign(x, xi, src, dst): return False
                else:  # 'X': xi[src] = dst  <=>  x[dst] = src
                    if not assign(x, xi, dst, src): return False
                st["forced"] += 1
                changed = True
        return True

    inverse_cache = {id(L): inverse(L) for L in word if not isinstance(L, str)}

    def dfs(x, xi):
        st["nodes"] += 1
        if st["nodes"] > node_budget: st["exhausted"] = True; return
        free = [j for j in range(N) if x[j] == NA]
        if not free:
            xt = tuple(x)
            if evaluate(word, xt) == c: sols.append(xt)
            return
        j = free[0] if order is None else min(free, key=order)
        st["decisions"] += 1
        for v in range(N):
            if xi[v] != NA: continue
            x2, xi2 = x[:], xi[:]
            if assign(x2, xi2, j, v) and propagate(x2, xi2):
                dfs(x2, xi2)
                if len(sols) >= max_solutions or st["exhausted"]: return

    x0 = [NA] * N; xi0 = [NA] * N
    if propagate(x0, xi0): dfs(x0, xi0)
    return sols, st


def wreath_instance(m, rng, even=True):
    """A CPWI round as x a x^2 = c in A_m wr C_2 <= S_{2m}; returns (a, c, x_planted)."""
    R = random_even if even else random_perm
    U, V, A, B = [R(m, rng) for _ in range(4)]
    def emb(P, Q):
        return tuple(list(P) + [q + m for q in Q])
    sig = tuple((j + m) % (2 * m) for j in range(2 * m))
    x = compose(emb(U, V), sig); a = emb(B, A)
    c = mul(x, a, x, x)
    return a, c, x


def generic_instance(N, rng):
    a = random_perm(N, rng); x = random_perm(N, rng)
    return a, mul(x, a, x, x), x
