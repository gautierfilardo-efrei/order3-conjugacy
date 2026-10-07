"""fourlayer_track.py -- the 'quatre couches' (M-layer) padding  Conj_3 -> P3  and an exhaustive
bitmask-domain solver for  x a x^2 = c  in Sym(N).

Conventions (as in the paper): permutations are image tuples, products right to left,
(pq)(j) = p(q(j)); the word x a x x applied to j is x(a(x(x(j)))).

Instance.  Omega = k blocks B_i = {iL,...,iL+L-1}, a = translation by 1 on every block (k disjoint
L-cycles), c = a^{m_i} on block i.  Padded instance on Omega x Z_M, point (j,r) encoded as r*kL + j:
    a'(j,r) = (a^{e_r} j, r+1),   c'(j,r) = (c^{e_r} j, r+1)      (rho = the M-cycle r -> r+1).

Solver.  Domains X[p] (bitmasks) for x(p); propagation = all-different (naked/hidden singles)
+ generalised arc consistency on the relaxed chain  p -> q=x(p) -> s=x(q) -> w=a'(s) -> x(w)=c'(p)
(three positions).  Branching on the smallest domain.  `complete` is True iff the whole search tree
was explored (no node/time budget hit), so the returned solution list is then exhaustive.
Standard library only.
"""
from __future__ import annotations
import itertools, time, json, sys

# ---------------------------------------------------------------- permutations
def comp(p, q): return tuple(p[q[j]] for j in range(len(p)))
def inv(p):
    r = [0]*len(p)
    for j, v in enumerate(p): r[v] = j
    return tuple(r)
def power(p, k):
    n = len(p); r = tuple(range(n))
    if k < 0: p = inv(p); k = -k
    for _ in range(k): r = comp(p, r)
    return r
def ctype(p):
    seen = [False]*len(p); t = []
    for j in range(len(p)):
        if not seen[j]:
            l = 0; v = j
            while not seen[v]: seen[v] = True; v = p[v]; l += 1
            t.append(l)
    return tuple(sorted(t, reverse=True))
def cycles(p):
    seen = [False]*len(p); out = []
    for j in range(len(p)):
        if not seen[j]:
            cyc = []; v = j
            while not seen[v]: seen[v] = True; cyc.append(v); v = p[v]
            out.append(cyc)
    return out
def is_p3_solution(a, c, x):
    """x a x x == c ?"""
    xx = comp(x, x)
    return comp(x, comp(a, xx)) == c
def order(p):
    n = len(p); r = p; o = 1
    idp = tuple(range(n))
    while r != idp: r = comp(p, r); o += 1
    return o

# ---------------------------------------------------------------- instances
def affine_instance(L, m):
    """a = k disjoint L-cycles (translation by 1 on each block), c = a^{m_i} on block i."""
    k = len(m); N = k*L
    a = [0]*N; c = [0]*N
    for i in range(k):
        for z in range(L):
            a[i*L+z] = i*L + (z+1) % L
            c[i*L+z] = i*L + (z+m[i]) % L
    return tuple(a), tuple(c)

def pad(a, c, e, rho=None):
    """Layered padding on Omega x Z_M (M = len(e)); point (j,r) -> r*N0 + j."""
    N0 = len(a); M = len(e)
    if rho is None: rho = [(r+1) % M for r in range(M)]
    apow = {E: power(a, E) for E in set(e)}
    cpow = {E: power(c, E) for E in set(e)}
    ap = [0]*(N0*M); cp = [0]*(N0*M)
    for r in range(M):
        for j in range(N0):
            ap[r*N0+j] = rho[r]*N0 + apow[e[r]][j]
            cp[r*N0+j] = rho[r]*N0 + cpow[e[r]][j]
    return tuple(ap), tuple(cp)

def lift(y, M):
    """y (x) id on Omega x Z_M."""
    N0 = len(y)
    return tuple(r*N0 + y[j] for r in range(M) for j in range(N0))

def conj3_affine_witnesses(L, m):
    """All order-3 conjugators y a y^-1 = c on the affine family (block lemma): y maps B_i -> B_sigma(i)
    by z -> m_sigma(i) z + t_i.  Enumerates sigma with sigma^3 = 1 and all translations (L^k); fine for small k."""
    k = len(m); a, c = affine_instance(L, m); out = []
    for sigma in itertools.permutations(range(k)):
        s3 = [sigma[sigma[sigma[i]]] for i in range(k)]
        if s3 != list(range(k)): continue
        for t in itertools.product(range(L), repeat=k):
            y = [0]*(k*L)
            for i in range(k):
                for z in range(L):
                    y[i*L+z] = sigma[i]*L + (m[sigma[i]]*z + t[i]) % L
            y = tuple(y)
            if comp(comp(y, a), inv(y)) == c and power(y, 3) == tuple(range(k*L)):
                out.append(y)
    return out

# ---------------------------------------------------------------- bitmask utilities
def bits(mask):
    while mask:
        b = mask & -mask
        yield b.bit_length()-1
        mask ^= b

def map_mask(mask, perm):
    """image of a set (bitmask) under a permutation."""
    r = 0
    while mask:
        b = mask & -mask
        r |= 1 << perm[b.bit_length()-1]
        mask ^= b
    return r

# ---------------------------------------------------------------- the solver
class P3Solver:
    def __init__(self, a, c):
        assert len(a) == len(c)
        self.a, self.c = tuple(a), tuple(c)
        self.N = len(a)
        self.ainv = inv(self.a); self.cinv = inv(self.c)
        self.full = (1 << self.N) - 1
        self.stats = {"nodes": 0, "complete": True, "time": 0.0}

    def propagate(self, X):
        N, a, c, ainv = self.N, self.a, self.c, self.ainv
        full = self.full
        while True:
            changed = False
            # --- all-different: naked singles
            for p in range(N):
                d = X[p]
                if d == 0: return None
                if d & (d-1) == 0:
                    for q in range(N):
                        if q != p and X[q] & d:
                            X[q] &= ~d; changed = True
                            if X[q] == 0: return None
            # inverse domains
            Xi = [0]*N
            for p in range(N):
                for v in bits(X[p]): Xi[v] |= 1 << p
            for v in range(N):
                if Xi[v] == 0: return None
                if Xi[v] & (Xi[v]-1) == 0:          # hidden single: only p can map to v
                    p = Xi[v].bit_length()-1
                    if X[p] != 1 << v:
                        X[p] = 1 << v; changed = True
            if changed: continue
            # --- chain  p -> q -> s -> w = a'(s) -> c'(p)
            # W[p] = feasible s for position 3 : a'(s) in Xi[c'(p)]
            W = [map_mask(Xi[c[p]], ainv) for p in range(N)]
            # S2[p] = union of X[q], q in X[p]  (possible x^2(p))
            S2 = [0]*N
            for p in range(N):
                u = 0
                for q in bits(X[p]): u |= X[q]
                S2[p] = u & W[p]
                if S2[p] == 0: return None
            # position 1: q in X[p] needs X[q] & W[p] != 0
            for p in range(N):
                d = X[p]; nd = 0
                for q in bits(d):
                    if X[q] & W[p]: nd |= 1 << q
                if nd != d:
                    X[p] = nd; changed = True
                    if nd == 0: return None
            # position 2: s in X[q] needs some p in Xi[q] with s in W[p]
            for q in range(N):
                u = 0
                for p in bits(Xi[q]): u |= W[p]
                nd = X[q] & u
                if nd != X[q]:
                    X[q] = nd; changed = True
                    if nd == 0: return None
            # position 3: x(w) = c'(p) needs a'^{-1}(w) in S2[p]
            T = [0]*N                             # T[s] = {p : s in S2[p]}
            for p in range(N):
                for s in bits(S2[p]): T[s] |= 1 << p
            for w in range(N):
                allowed = map_mask(T[ainv[w]], c)
                nd = X[w] & allowed
                if nd != X[w]:
                    X[w] = nd; changed = True
                    if nd == 0: return None
            if not changed: return X

    def solve(self, max_solutions=None, node_budget=None, time_cap=None, X0=None, verbose=False):
        N = self.N
        sols = []
        st = self.stats = {"nodes": 0, "complete": True, "time": 0.0}
        t0 = time.time()
        a, c = self.a, self.c

        def dfs(X):
            st["nodes"] += 1
            if node_budget is not None and st["nodes"] > node_budget: st["complete"] = False; return
            if time_cap is not None and time.time() - t0 > time_cap: st["complete"] = False; return
            if verbose and st["nodes"] % 2000 == 0:
                print(f"  nodes={st['nodes']} sols={len(sols)} t={time.time()-t0:.0f}s", flush=True)
            # choose variable
            best = -1; bsz = N+1
            for p in range(N):
                d = X[p]; sz = bin(d).count("1")
                if 1 < sz < bsz: best, bsz = p, sz
            if best < 0:
                x = tuple(X[p].bit_length()-1 for p in range(N))
                if is_p3_solution(a, c, x): sols.append(x)
                return
            for v in bits(X[best]):
                Y = X[:]; Y[best] = 1 << v
                Y = self.propagate(Y)
                if Y is not None:
                    dfs(Y)
                    if max_solutions is not None and len(sols) >= max_solutions: return
                    if not st["complete"]: return

        X = [self.full]*N if X0 is None else list(X0)
        X = self.propagate(X)
        if X is not None: dfs(X)
        st["time"] = time.time() - t0
        return sols, dict(st)

# ---------------------------------------------------------------- analysis of solutions
def layer_structure(x, N0, M):
    """Return (sigma_layers, per-layer maps) if x maps every layer Omega x {r} onto a layer, else None."""
    lay = []
    for r in range(M):
        imgs = {x[r*N0+j] // N0 for j in range(N0)}
        if len(imgs) != 1: return None
        lay.append(imgs.pop())
    return lay

def is_layer_affine(x, L, k, M):
    """x(j,r) = (affine block map, r + const)?  Returns dict or None."""
    N0 = k*L
    lay = layer_structure(x, N0, M)
    if lay is None: return None
    shifts = {(lay[r]-r) % M for r in range(M)}
    if len(shifts) != 1: return {"layers": lay, "shift": None, "affine": False}
    info = {"layers": lay, "shift": shifts.pop(), "affine": True, "maps": []}
    for r in range(M):
        for i in range(k):
            imgs = [x[r*N0 + i*L + z] % N0 for z in range(L)]
            blocks = {v // L for v in imgs}
            if len(blocks) != 1: info["affine"] = False; return info
            zs = [v % L for v in imgs]
            s = (zs[1]-zs[0]) % L
            if any((zs[z]-zs[0]) % L != s*z % L for z in range(L)): info["affine"] = False; return info
            info["maps"].append((r, i, blocks.pop(), s, zs[0]))
    return info

def ord_mod(g, L):
    g %= L; r = g; o = 1
    while r != 1 % L:
        r = r*g % L; o += 1
        if o > L: return None
    return o
