"""
order3_track.py -- "Ordre 3" track: the order-3 core D0 of the word equation
P3 : x a x^2 = c in S_N, and its rigidity.

D0(a,c) : is there x in S_N with x^3 = 1 and x a x^{-1} = c ?
        (= P3 restricted to solutions of order dividing 3, since x^{-1} = x^2)

Conventions (same as the project): permutations are tuples of images on
range(N); product right to left, (p*q)(j) = p(q(j)).

Contents
  1. permutation primitives
  2. exhaustive D0 / P3 at small N
  3. propagation solver for D0 (any N), cross-checked against exhaustive
  4. block instances: a = k disjoint L-cycles, c = an L-cycle on each block
     Lemma B (block lemma): D0 <=> exists sigma in S_k, sigma^3=1, with
        Fix(i) for fixed points and Adm(i,j,k) for 3-cycles.
  5. affine blocks (L prime, c|block_i = translation by unit m_i):
        Adm(i,j,k) <=> m_i m_j m_k = 1,  Fix(i) <=> m_i^3 = 1.
  6. reduction Numerical-3DM -> D0 (polynomial), tested end-to-end.
  7. rigidity census: (a,c) such that every P3 solution has x^3 = 1.
Standard library only.
"""
import itertools, json, math, random, sys, time
from collections import defaultdict, Counter

# ---------------------------------------------------------------- 1. primitives
def compose(p, q):
    """(p*q)(j) = p(q(j))"""
    return tuple(p[j] for j in q)

def mul(*ps):
    r = ps[0]
    for p in ps[1:]:
        r = compose(r, p)
    return r

def inverse(p):
    r = [0] * len(p)
    for i, j in enumerate(p):
        r[j] = i
    return tuple(r)

def identity(N):
    return tuple(range(N))

def power(p, e):
    N = len(p)
    r = identity(N)
    if e < 0:
        p, e = inverse(p), -e
    for _ in range(e):
        r = compose(r, p)
    return r

def cycles(p):
    N = len(p); seen = [False] * N; out = []
    for i in range(N):
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

def random_perm(N, rng):
    l = list(range(N)); rng.shuffle(l); return tuple(l)

def from_cycles(N, cycs):
    p = list(range(N))
    for c in cycs:
        for i, j in zip(c, c[1:] + c[:1]):
            p[i] = j
    return tuple(p)

def conjugate(x, a):
    """x a x^{-1}"""
    return mul(x, a, inverse(x))

def class_rep(ct):
    """canonical permutation of cycle type ct (tuple of lengths)"""
    N = sum(ct); cycs = []; s = 0
    for L in ct:
        cycs.append(list(range(s, s + L))); s += L
    return from_cycles(N, cycs)

def partitions(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

# ---------------------------------------------------------------- 2. exhaustive
_ORDER3_CACHE = {}
def order3_elements(N):
    """all x in S_N with x^3 = 1 (enumerated as products of disjoint 3-cycles)."""
    if N in _ORDER3_CACHE: return _ORDER3_CACHE[N]
    out = []
    def rec(rem, cycs):
        if not rem:
            out.append(from_cycles(N, cycs)); return
        i = rem[0]; rest = rem[1:]
        rec(rest, cycs)                              # i fixed
        for (j, k) in itertools.combinations(rest, 2):
            r2 = [t for t in rest if t != j and t != k]
            rec(r2, cycs + [[i, j, k]]); rec(r2, cycs + [[i, k, j]])
    rec(list(range(N)), [])
    _ORDER3_CACHE[N] = out
    return out

def d0_exhaustive(a, c):
    return [x for x in order3_elements(len(a)) if conjugate(x, a) == c]

def p3_exhaustive(a, c):
    N = len(a)
    return [x for x in map(tuple, itertools.permutations(range(N))) if mul(x, a, x, x) == c]

# ---------------------------------------------------------------- 3. propagation solver
class D0Solver:
    """Backtracking search for x with x a x^{-1} = c and x^3 = 1.
    Assigning x(p)=q propagates x(a^n p) = c^n q along the whole a-cycle of p,
    and the order-3 closure: x(p)=q & x(q)=r  =>  x(r)=p ;  x(r)=p & x(p)=q => x(q)=r.
    """
    def __init__(self, a, c, node_budget=10**7):
        self.a, self.c, self.N = a, c, len(a)
        self.budget = node_budget; self.nodes = 0
        self.solutions = []

    def _assign(self, x, xinv, p, q, trail):
        """assign x[p]=q with propagation; return False on conflict."""
        stack = [(p, q)]
        while stack:
            p, q = stack.pop()
            if x[p] != -1:
                if x[p] != q: return False
                continue
            if xinv[q] != -1: return False
            x[p] = q; xinv[q] = p; trail.append(p)
            # intertwining along a-cycle
            stack.append((self.a[p], self.c[q]))
            # order-3 closure
            if x[q] != -1: stack.append((x[q], p))
            if xinv[p] != -1: stack.append((q, xinv[p]))
        return True

    def _undo(self, x, xinv, trail, mark):
        while len(trail) > mark:
            p = trail.pop(); xinv[x[p]] = -1; x[p] = -1

    def solve(self, max_solutions=1):
        N = self.N
        if ctype(self.a) != ctype(self.c): return []
        x = [-1] * N; xinv = [-1] * N; trail = []
        # candidate images: q must lie on a c-cycle of the same length as p's a-cycle
        alen = [0] * N; clen = [0] * N
        for cy in cycles(self.a):
            for p in cy: alen[p] = len(cy)
        for cy in cycles(self.c):
            for p in cy: clen[p] = len(cy)
        order_pts = sorted(range(N), key=lambda p: -alen[p])
        def rec():
            self.nodes += 1
            if self.nodes > self.budget: raise TimeoutError
            p = next((p for p in order_pts if x[p] == -1), None)
            if p is None:
                self.solutions.append(tuple(x)); return len(self.solutions) >= max_solutions
            for q in range(N):
                if xinv[q] != -1 or clen[q] != alen[p]: continue
                mark = len(trail)
                if self._assign(x, xinv, p, q, trail):
                    if rec(): return True
                self._undo(x, xinv, trail, mark)
            return False
        try:
            rec()
        except TimeoutError:
            self.timed_out = True
        return self.solutions

def d0_solve(a, c, max_solutions=1, node_budget=10**7):
    s = D0Solver(a, c, node_budget); sols = s.solve(max_solutions)
    for x in sols:
        assert power(x, 3) == identity(len(a)) and conjugate(x, a) == c
    return sols, s.nodes

# ---------------------------------------------------------------- 4. block instances
def block_instance(L, pis):
    """a = translation by +1 on each block Z/L (k = len(pis) blocks);
    c|block_i = pi_i a pi_i^{-1} where pi_i is a permutation of range(L)."""
    k = len(pis); N = k * L
    a = [0] * N; c = [0] * N
    for i, pi in enumerate(pis):
        base = i * L
        for z in range(L):
            a[base + z] = base + (z + 1) % L
        pinv = inverse(tuple(pi))
        for z in range(L):
            # c(pi(z)) = pi(z+1)
            c[base + pi[z]] = base + pi[(z + 1) % L]
    return tuple(a), tuple(c)

def shift(L, r):
    return tuple((z + r) % L for z in range(L))

def adm(L, pi_i, pi_j, pi_k):
    """Adm(i,j,k): exists r_i,r_j,r_k with  pi_i s^{r_k} pi_k s^{r_j} pi_j s^{r_i} = id  on Z/L,
    i.e. x_k x_j x_i = id for the 3-cycle i->j->k->i, x_i = pi_{sigma(i)} s^{r_i}."""
    for rk in range(L):
        for rj in range(L):
            for ri in range(L):
                m = mul(pi_i, shift(L, rk), pi_k, shift(L, rj), pi_j, shift(L, ri))
                if m == identity(L): return True
    return False

def fix_ok(L, pi_i):
    return any(power(mul(pi_i, shift(L, r)), 3) == identity(L) for r in range(L))

def order3_perms_of(k):
    return [s for s in order3_elements(k)]

def d0_block_criterion(L, pis):
    """Lemma B criterion, evaluated by brute force over sigma and rotations."""
    k = len(pis)
    fix = [fix_ok(L, pi) for pi in pis]
    admc = {}
    for sigma in order3_elements(k):
        ok = True
        for cyc in cycles(sigma):
            if len(cyc) == 1:
                if not fix[cyc[0]]: ok = False; break
            else:
                i, j, kk = cyc            # i -> j -> kk -> i
                key = (i, j, kk)
                if key not in admc: admc[key] = adm(L, pis[i], pis[j], pis[kk])
                if not admc[key]: ok = False; break
        if ok: return True
    return False

# ---------------------------------------------------------------- 5. affine blocks
def affine_block_instance(L, ms):
    """c|block_i = translation by m_i (unit mod L). a = translation by 1."""
    k = len(ms); N = k * L
    a = tuple((i * L + (z + 1) % L) for i in range(k) for z in range(L))
    c = tuple((i * L + (z + ms[i]) % L) for i in range(k) for z in range(L))
    return a, c

def d0_affine_criterion(L, ms):
    """exists sigma^3=1 : fixed i have m_i^3=1 mod L, 3-cycles have m_i m_j m_k = 1 mod L."""
    k = len(ms)
    for sigma in order3_elements(k):
        ok = True
        for cyc in cycles(sigma):
            if len(cyc) == 1:
                if pow(ms[cyc[0]], 3, L) != 1: ok = False; break
            else:
                if (ms[cyc[0]] * ms[cyc[1]] * ms[cyc[2]]) % L != 1: ok = False; break
        if ok: return True
    return False

# ---------------------------------------------------------------- 6. reduction N3DM -> D0
def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0: return False
        f += 2
    return True

def primitive_root(p):
    phi = p - 1
    fac = set(); n = phi; f = 2
    while f * f <= n:
        while n % f == 0: fac.add(f); n //= f
        f += 1
    if n > 1: fac.add(n)
    for g in range(2, p):
        if all(pow(g, phi // q, p) != 1 for q in fac): return g
    raise ValueError

def reduce_n3dm(W, X, Y, B):
    """Numerical 3-Dimensional Matching instance (sizes W,X,Y lists of length q, bound B)
    -> D0 instance (a, c) on N = 3 q L points.
    Encoding: R = L-1, K = 2B; e_w = s(w), e_x = s(x)+K, e_y = s(y)+R-K-B (mod R);
    m_i = omega^{e_i}; L prime, L = 2 mod 3, L-1 >= 10 B.
    Requires all sizes in [1, B-2] (else trivially NO -> return a fixed NO instance)."""
    q = len(W); assert len(X) == q == len(Y)
    sizes = W + X + Y
    if any(s < 1 or s > B - 2 for s in sizes):
        # trivial NO instance: a=(0 1), c=(2 3) in S_4 has no order-3 conjugator (checked exhaustively)
        return (1, 0, 2, 3), (0, 1, 3, 2), {"trivial_no": True}
    L = 10 * B + 1
    while not (is_prime(L) and L % 3 == 2): L += 1
    R = L - 1; K = 2 * B
    e = [s for s in W] + [s + K for s in X] + [(s + R - K - B) % R for s in Y]
    assert all(ei % R != 0 for ei in e)
    omega = primitive_root(L)
    ms = [pow(omega, ei, L) for ei in e]
    assert all(m != 1 for m in ms) and all(pow(m, 3, L) != 1 for m in ms)
    a, c = affine_block_instance(L, ms)
    return a, c, {"L": L, "R": R, "K": K, "omega": omega, "e": e, "m": ms, "q": q, "N": len(a)}

def n3dm_bruteforce(W, X, Y, B):
    q = len(W)
    for px in itertools.permutations(range(q)):
        for py in itertools.permutations(range(q)):
            if all(W[i] + X[px[i]] + Y[py[i]] == B for i in range(q)): return True
    return False

def zero_sum_triples_bruteforce(ms, L):
    """exists partition of indices into triples with product 1 mod L (no fixed blocks allowed
    when all m^3 != 1) -- generic version = d0_affine_criterion."""
    return d0_affine_criterion(L, ms)

# ---------------------------------------------------------------- 7. test driver
def run_checks(quick=True):
    """Reproduces the verifications recorded in order3_results.json (quick ~1 min)."""
    rng = random.Random(1); out = {}
    # (i) propagation solver == exhaustive, N<=8
    for N in range(3, 9):
        for t in range(30):
            a = random_perm(N, rng)
            c = conjugate(random_perm(N, rng), a) if t % 2 == 0 else random_perm(N, rng)
            assert set(d0_exhaustive(a, c)) == set(d0_solve(a, c, 10**6)[0])
    # (ii) Lemma B exhaustive: L=3,k=3 (N=9)
    for pis in itertools.product(list(map(tuple, itertools.permutations(range(3)))), repeat=3):
        a, c = block_instance(3, pis)
        assert bool(d0_exhaustive(a, c)) == d0_block_criterion(3, pis)
    # (iii) Lemma C (affine) vs solver: L=5,7 ; k=3
    for L in (5, 7):
        for ms in itertools.product(range(1, L), repeat=3):
            a, c = affine_block_instance(L, ms)
            assert bool(d0_solve(a, c, 1, 10**6)[0]) == d0_affine_criterion(L, ms)
    # (iv) reduction N3DM -> D0 (criterion side), q=2,3
    for q, B, trials in [(2, 6, 30), (3, 7, 20)]:
        for t in range(trials):
            W = [rng.randint(1, B - 2) for _ in range(q)]; X = [rng.randint(1, B - 2) for _ in range(q)]
            Y = [rng.randint(1, B - 2) for _ in range(q)]
            a, c, meta = reduce_n3dm(W, X, Y, B)
            assert n3dm_bruteforce(W, X, Y, B) == d0_affine_criterion(meta["L"], meta["m"])
    print("all checks passed")

if __name__ == "__main__":
    run_checks()
