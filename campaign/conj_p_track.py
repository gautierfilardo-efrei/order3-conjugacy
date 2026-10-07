#!/usr/bin/env python3
"""conj_p_track.py -- Conj_p track: order-p conjugacy in the implicit conjugator coset.

Conj_p(a,c): given a,c in Sym(N), is there x with x^p = 1 and x a x^{-1} = c ?
Conventions: permutations are image tuples on {0..N-1}, (pq)(j) = p(q(j)).
Standard library only.  Run:  python conj_p_track.py [--quick]
Writes conj_p_results.json.
"""
import itertools, json, math, random, sys, time
from collections import Counter

QUICK = '--quick' in sys.argv
T0 = time.time()
def log(*a):
    print('[%6.1fs]' % (time.time()-T0), *a); sys.stdout.flush()

# ---------------- permutation primitives ----------------
def compose(p, q):            # (pq)(j) = p(q(j))
    return tuple(p[q[j]] for j in range(len(q)))
def inverse(p):
    inv = [0]*len(p)
    for i, j in enumerate(p): inv[j] = i
    return tuple(inv)
def identity(n): return tuple(range(n))
def ppow(p, e):
    r = identity(len(p))
    for _ in range(e): r = compose(p, r)
    return r
def conj(x, a):               # x a x^{-1}
    return compose(compose(x, a), inverse(x))
def perm_cycles(p):
    n = len(p); seen = [False]*n; out = []
    for s in range(n):
        if not seen[s]:
            cyc = []; u = s
            while not seen[u]:
                seen[u] = True; cyc.append(u); u = p[u]
            out.append(cyc)
    return out
def ctype(p):
    return tuple(sorted((len(c) for c in perm_cycles(p)), reverse=True))
def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    f = 3
    while f*f <= n:
        if n % f == 0: return False
        f += 2
    return True
def mult_order(g, L):
    o, x = 1, g % L
    while x != 1:
        x = x*g % L; o += 1
    return o

# ---------------- blockwise family ----------------
def family(L, k, m):
    """a = k disjoint L-cycles (translation z->z+1 on block B_i), c|B_i = a^{m_i}."""
    N = L*k; a = [0]*N; c = [0]*N
    for i in range(k):
        for z in range(L):
            a[i*L+z] = i*L + (z+1) % L
            c[i*L+z] = i*L + (z+m[i]) % L
    return tuple(a), tuple(c)
def affine_perm(L, k, m, sigma, t):
    x = [0]*(L*k)
    for i in range(k):
        for z in range(L):
            x[i*L+z] = sigma[i]*L + (m[sigma[i]]*z + t[i]) % L
    return tuple(x)
def as_affine(x, L, k, m):
    """Return (sigma,t) if x is blockwise affine with slopes m_{sigma(i)}, else None."""
    sigma = []; t = []
    for i in range(k):
        j, r = divmod(x[i*L], L)
        for z in range(L):
            jj, rr = divmod(x[i*L+z], L)
            if jj != j or rr != (m[j]*z + r) % L: return None
        sigma.append(j); t.append(r)
    return tuple(sigma), tuple(t)
def affine_criterion_sigma(p, L, m, sigma):
    """sigma^p = 1; on each p-cycle the slope product is 1 mod L; on each fixed block m_i^p = 1 mod L."""
    for cyc in perm_cycles(sigma):
        if len(cyc) == 1:
            if pow(m[cyc[0]], p, L) != 1: return False
        elif len(cyc) == p:
            pr = 1
            for i in cyc: pr = pr*m[i] % L
            if pr != 1: return False
        else:
            return False
    return True
def criterion_decides(p, L, k, m):
    return any(affine_criterion_sigma(p, L, m, s) for s in itertools.permutations(range(k)))
def units(L): return [u for u in range(1, L) if math.gcd(u, L) == 1]

# ---------------- Verification A: brute force over Sym(N), N <= 9 ----------------
def verifyA(p, L, k, m):
    a, c = family(L, k, m); N = L*k; ident = identity(N)
    conjugators = []
    for x in itertools.permutations(range(N)):
        ok = True
        for j in range(N):
            if x[a[j]] != c[x[j]]: ok = False; break
        if ok: conjugators.append(x)
    non_affine = sum(1 for x in conjugators if as_affine(x, L, k, m) is None)
    direct_sols = [x for x in conjugators if ppow(x, p) == ident]
    sig_direct = set(as_affine(x, L, k, m)[0] for x in direct_sols)
    sig_crit = set(s for s in itertools.permutations(range(k)) if affine_criterion_sigma(p, L, m, s))
    return dict(p=p, L=L, k=k, m=list(m), N=N, n_conjugators=len(conjugators),
                expected_conjugators=L**k*math.factorial(k), non_affine=non_affine,
                n_solutions=len(direct_sols), sigmas_agree=(sig_direct == sig_crit),
                exists_direct=bool(direct_sols), exists_criterion=bool(sig_crit),
                all_fpf_order_p=all(ctype(x) == tuple([p]*(N//p)) for x in direct_sols) if direct_sols else None)

# ---------------- Verification B: exhaustive (sigma, t) enumeration ----------------
def verifyB(p, L, k, m):
    N = L*k; ident = identity(N); mism = 0; n_sig_ok = 0; n_sol = 0
    for sigma in itertools.permutations(range(k)):
        crit = affine_criterion_sigma(p, L, m, sigma)
        found = False
        for t in itertools.product(range(L), repeat=k):
            x = affine_perm(L, k, m, sigma, t)
            if ppow(x, p) == ident:
                found = True; n_sol += 1
                if not crit: break
        if found != crit: mism += 1
        n_sig_ok += crit
    return dict(p=p, L=L, k=k, m=list(m), N=N, enumerated=L**k*math.factorial(k),
                mismatching_sigmas=mism, admissible_sigmas=n_sig_ok, n_solutions=n_sol)

# ---------------- The reduction N3DM -> Conj_p ----------------
def choose_modulus(p, T0, kind):
    """kind='prime': smallest prime L with L-1 >= T0 and p not dividing L-1 (g = primitive root, T = L-1).
       kind='pow2' : L = 2^r with T = 2^(r-2) >= T0, g = 5."""
    if kind == 'prime':
        L = T0 + 1
        while not (is_prime(L) and (L-1) % p != 0): L += 1
        T = L - 1
        g = next(u for u in range(2, L) if mult_order(u, L) == T)
        return L, T, g
    r = 3
    while 2**(r-2) < T0: r += 1
    return 2**r, 2**(r-2), 5

def reduction(p, sizes, B, kind='prime'):
    """sizes: list of p lists (classes 1..p), each of q integers in [1,B].  Returns instance data.
       (Padding from N3DM: classes 4..p all of size 1 and B' = B + p - 3 -- done by the caller.)"""
    q = len(sizes[0]); k = p*q
    M = p*B + 1
    w = [(p+1)**j for j in range(p-1)]          # w_1..w_{p-1} = 1, p+1, (p+1)^2, ...
    W = sum(w)
    T0 = M*(p+1)**(p-1)
    L, T, g = choose_modulus(p, T0, kind)
    e = []; cls = []; s_all = []
    for j in range(p):
        for u in range(q):
            s = sizes[j][u]
            if j < p-1: ee = (s + M*w[j]) % T
            else:       ee = (s - M*W - B) % T
            e.append(ee); cls.append(j); s_all.append(s)
    m = [pow(g, ee, L) for ee in e]
    return dict(p=p, q=q, k=k, B=B, M=M, T0=T0, L=L, T=T, g=g, e=e, m=m, cls=cls, s=s_all, kind=kind)

def admissible_sets_check(inst):
    """all p-subsets of blocks: (sum e == 0 mod T)  <=>  (one per class and s-sum == B)."""
    p, k, T = inst['p'], inst['k'], inst['T']
    bad = 0; n_adm = 0
    for S in itertools.combinations(range(k), p):
        adm = sum(inst['e'][i] for i in S) % T == 0
        match = (sorted(inst['cls'][i] for i in S) == list(range(p))) and sum(inst['s'][i] for i in S) == inst['B']
        if adm != match: bad += 1
        n_adm += adm
    return bad, n_adm, (min(inst['e']) > 0)

def matching_exists(sizes, B):
    """brute-force Numerical p-DM: partition into p-tuples (one per class) of size-sum B."""
    p = len(sizes); q = len(sizes[0])
    def rec(used, u0):
        if u0 == q: return True
        # element u0 of class 0 must be matched
        for combo in itertools.product(*[[v for v in range(q) if not used[j][v]] for j in range(1, p)]):
            if sizes[0][u0] + sum(sizes[j][combo[j-1]] for j in range(1, p)) == B:
                for j in range(1, p): used[j][combo[j-1]] = True
                if rec(used, u0+1): return True
                for j in range(1, p): used[j][combo[j-1]] = False
        return False
    return rec([[False]*q for _ in range(p)], 0)

def admissible_partition(inst):
    """brute force: partition blocks into p-sets with exponent sum 0 mod T; returns partition or None."""
    p, k, T, e = inst['p'], inst['k'], inst['T'], inst['e']
    blocks = list(range(k))
    def rec(rem):
        if not rem: return []
        first = rem[0]
        for rest in itertools.combinations(rem[1:], p-1):
            S = (first,) + rest
            if sum(e[i] for i in S) % T == 0:
                sub = rec([b for b in rem if b not in S])
                if sub is not None: return [S] + sub
        return None
    return rec(blocks)

def certificate_from_partition(inst, part):
    """Build x of cycle type p^{N/p}: sigma cycles the blocks of each admissible p-set; translations chosen
       so that x^p = 1 (solve the single linear condition per cycle)."""
    p, k, L, m = inst['p'], inst['k'], inst['L'], inst['m']
    sigma = [0]*k; t = [0]*k
    for S in part:
        for idx in range(p): sigma[S[idx]] = S[(idx+1) % p]
    # translation: on a p-cycle i1->i2->...->ip->i1, x^p|B_{i1}(0) = sum_{r} (prod of later slopes) t_{i_r}; set all t=0 except last chosen
    # simplest: with all t=0, x^p|B_{i1} is z -> (prod m) z = z.  So t = 0 works.
    x = affine_perm(L, k, m, tuple(sigma), tuple(t))
    return x

def random_n3dm(q, B, planted, rng):
    if planted:
        W, X, Y = [], [], []
        for _ in range(q):
            w = rng.randint(1, B-2); xx = rng.randint(1, B-1-w); y = B - w - xx
            W.append(w); X.append(xx); Y.append(y)
        rng.shuffle(X); rng.shuffle(Y)
        return [W, X, Y]
    return [[rng.randint(1, B) for _ in range(q)] for _ in range(3)]

def pad(sizes3, B, p):
    q = len(sizes3[0])
    return sizes3 + [[1]*q for _ in range(p-3)], B + (p-3)

# ---------------- Conj_2: polynomial algorithm ----------------
def g_orbits(a, c):
    N = len(a); idx = [-1]*N; orbs = []
    for s0 in range(N):
        if idx[s0] < 0:
            k = len(orbs); orb = [s0]; idx[s0] = k; st = [s0]
            while st:
                u = st.pop()
                for g in (a, c):
                    w = g[u]
                    if idx[w] < 0: idx[w] = k; orb.append(w); st.append(w)
            orbs.append(orb)
    return orbs, idx

def equivariant_map(a, c, j0, v, twisted):
    """Partial bijection x on the <a,c>-orbit of j0 with x(j0)=v and x g = g' x for (g,g') in pairs; None if inconsistent."""
    pairs = ((a, c), (c, a)) if twisted else ((a, a), (c, c))
    x = {j0: v}; st = [j0]
    while st:
        u = st.pop()
        for g, gp in pairs:
            u2 = g[u]; v2 = gp[x[u]]
            if u2 in x:
                if x[u2] != v2: return None
            else:
                x[u2] = v2; st.append(u2)
    if len(set(x.values())) != len(x): return None
    return x

def conj2_solve(a, c):
    """Return an involution x with x a x^{-1} = c, or None.  Polynomial (O(N^3))."""
    N = len(a)
    orbs, idx = g_orbits(a, c); r = len(orbs)
    # plain isomorphism classes of orbits (as <a,c>-sets)
    cls = [-1]*r; reps = []
    for i in range(r):
        for ci, rep in enumerate(reps):
            if len(orbs[rep]) == len(orbs[i]) and any(equivariant_map(a, c, orbs[rep][0], v, False) for v in orbs[i]):
                cls[i] = ci; break
        else:
            cls[i] = len(reps); reps.append(i)
    members = [[i for i in range(r) if cls[i] == ci] for ci in range(len(reps))]
    # twisted partner class
    partner = []
    for ci, rep in enumerate(reps):
        pj = None
        for cj, rep2 in enumerate(reps):
            if len(orbs[rep]) == len(orbs[rep2]) and any(equivariant_map(a, c, orbs[rep][0], v, True) for v in orbs[rep2]):
                pj = cj; break
        if pj is None: return None
        partner.append(pj)
    x = [None]*N
    def install(O, xm):            # xm: dict on O; extend by inverse on image
        for u, v in xm.items():
            x[u] = v; x[v] = u
    done = set()
    for ci in range(len(reps)):
        if ci in done: continue
        cj = partner[ci]
        if cj != ci:
            if len(members[ci]) != len(members[cj]): return None
            for i, j in zip(members[ci], members[cj]):
                O = orbs[i]
                xm = next((mm for v in orbs[j] for mm in [equivariant_map(a, c, O[0], v, True)] if mm), None)
                install(O, xm)
            done.update({ci, cj})
        else:
            mem = members[ci]
            for i, j in zip(mem[0::2], mem[1::2]):
                O = orbs[i]
                xm = next((mm for v in orbs[j] for mm in [equivariant_map(a, c, O[0], v, True)] if mm), None)
                install(O, xm)
            if len(mem) % 2 == 1:
                O = orbs[mem[-1]]; j0 = O[0]; found = None
                for v in O:
                    mm = equivariant_map(a, c, j0, v, True)
                    if mm is not None and mm[mm[j0]] == j0: found = mm; break
                if found is None: return None
                install(O, found)
            done.add(ci)
    return tuple(x)

def conj2_decide(a, c): return conj2_solve(a, c) is not None

def involutions(N):
    """all x in Sym(N) with x^2 = 1 (including identity)."""
    out = []
    def rec(x, free):
        if not free: out.append(tuple(x)); return
        u = free[0]
        x[u] = u; rec(x, free[1:])
        for v in free[1:]:
            x[u] = v; x[v] = u
            rec(x, [w for w in free[1:] if w != v])
            x[v] = v
    rec(list(range(N)), list(range(N)))
    return out

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-first, first): yield (first,) + rest
def class_rep(lam):
    a = []; s = 0
    for ln in lam:
        a += [s + (i+1) % ln for i in range(ln)]; s += ln
    return tuple(a)

def verify_conj2_exhaustive(N):
    invs = involutions(N); allperms = list(itertools.permutations(range(N)))
    total = 0; mism = 0; n_yes = 0; cert_fail = 0
    for lam in partitions(N):
        a = class_rep(lam)
        yes = set(conj(x, a) for x in invs)
        klass = set(conj(h, a) for h in allperms)
        for c in klass:
            x = conj2_solve(a, c)
            dec = x is not None
            if dec:
                if ppow(x, 2) != identity(N) or conj(x, a) != c: cert_fail += 1
            if dec != (c in yes): mism += 1
            total += 1; n_yes += (c in yes)
    return dict(N=N, pairs_tested=total, yes_instances=n_yes, mismatches=mism, certificate_failures=cert_fail,
                note='one representative a per conjugacy class, all c in the class (WLOG by simultaneous conjugation)')

def verify_conj2_random(N, trials, rng):
    invs = involutions(N); mism = 0; n_yes = 0; cert_fail = 0
    for t in range(trials):
        a = tuple(rng.sample(range(N), N))
        if t % 2 == 0:
            c = conj(rng.choice(invs), a)
        else:
            h = tuple(rng.sample(range(N), N)); c = conj(h, a)
        truth = any(conj(x, a) == c for x in invs)
        x = conj2_solve(a, c); dec = x is not None
        if dec and (ppow(x, 2) != identity(N) or conj(x, a) != c): cert_fail += 1
        mism += (dec != truth); n_yes += truth
    return dict(N=N, trials=trials, yes_instances=n_yes, mismatches=mism, certificate_failures=cert_fail)

# ---------------- main ----------------
def main():
    rng = random.Random(20261007)
    R = {'quick': QUICK}
    # ---- A: brute force over Sym(N) ----
    log('A: brute force over Sym(N)')
    A = []
    casesA = [(3, 5, 1), (3, 7, 1), (3, 4, 2), (3, 8, 1), (5, 3, 2), (5, 3, 3), (5, 7, 1), (5, 4, 2), (2, 3, 2), (2, 4, 2), (2, 5, 1), (7, 3, 3), (7, 8, 1), (3, 3, 3)]
    if QUICK: casesA = [cs for cs in casesA if cs[1]**cs[2]*1 <= 8 or cs[1]*cs[2] <= 8]
    for (p, L, k) in casesA:
        U = units(L)
        mlist = [tuple(rng.choice(U) for _ in range(k)) for _ in range(2 if L*k >= 9 else 3)]
        if k == 1: mlist = [(u,) for u in U]
        for m in mlist:
            A.append(verifyA(p, L, k, m))
        log('  A done', p, L, k)
    R['A_bruteforce_SymN'] = A
    R['A_summary'] = dict(instances=len(A),
                          block_lemma_ok=all(r['non_affine'] == 0 and r['n_conjugators'] == r['expected_conjugators'] for r in A),
                          criterion_ok=all(r['sigmas_agree'] and r['exists_direct'] == r['exists_criterion'] for r in A))
    json.dump(R, open('conj_p_results.json', 'w'), indent=1)
    # ---- B: exhaustive affine enumeration ----
    log('B: exhaustive (sigma,t) enumeration')
    Bres = []
    casesB = [(3, 5, 4, 3), (3, 11, 3, 3), (3, 7, 4, 2), (3, 5, 5, 2), (5, 3, 5, 3), (5, 4, 5, 2), (5, 7, 4, 2), (5, 11, 3, 2), (2, 7, 4, 2), (7, 3, 4, 1)]
    if QUICK: casesB = [(3, 5, 4, 2), (5, 3, 5, 2), (5, 4, 5, 1), (2, 7, 4, 1)]
    for (p, L, k, ntr) in casesB:
        U = units(L)
        for tr in range(ntr):
            # bias towards instances where admissible sigmas exist
            m = [rng.choice(U) for _ in range(k)]
            if tr == 0 and k >= p:
                pr = 1
                for i in range(p-1): pr = pr*m[i] % L
                m[p-1] = pow(pr, -1, L)
            Bres.append(verifyB(p, L, k, tuple(m)))
        log('  B done', p, L, k)
    R['B_affine_enumeration'] = Bres
    R['B_summary'] = dict(instances=len(Bres), total_enumerated=sum(r['enumerated'] for r in Bres),
                          criterion_ok=all(r['mismatching_sigmas'] == 0 for r in Bres))
    json.dump(R, open('conj_p_results.json', 'w'), indent=1)
    # ---- C: the reduction ----
    log('C: reduction N3DM -> Conj_p')
    C = []
    ntr = 6 if QUICK else 24
    for p in (3, 5):
        ncert = {'prime': 0, 'pow2': 0}
        for tr in range(ntr):
            q = rng.choice([2, 3]) if p == 3 else 2
            B = rng.choice([3, 4, 5, 6])
            sizes3 = random_n3dm(q, B, planted=(tr % 2 == 0), rng=rng)
            truth = matching_exists(sizes3, B)
            sizes, Bp = pad(sizes3, B, p)
            assert matching_exists(sizes, Bp) == truth
            for kind in ('prime', 'pow2'):
                inst = reduction(p, sizes, Bp, kind)
                bad, n_adm, nonzero = admissible_sets_check(inst)
                part = admissible_partition(inst)
                dec = part is not None
                assert bad == 0 and nonzero, (p, sizes, kind)
                assert dec == truth
                rec = dict(p=p, q=q, B=B, Bp=Bp, kind=kind, L=inst['L'], T=inst['T'], N=inst['k']*inst['L'],
                           truth=truth, reduction_decides=dec, bad_subsets=bad, admissible_subsets=n_adm, all_m_nontrivial=nonzero)
                if dec and ncert[kind] < (1 if QUICK else 3):
                    ncert[kind] += 1
                    a, c = family(inst['L'], inst['k'], inst['m'])
                    x = certificate_from_partition(inst, part)
                    rec['certificate_checked'] = (ppow(x, p) == identity(len(a))) and (conj(x, a) == c) and ctype(x) == tuple([p]*(len(a)//p))
                    rec['x_cycle_type'] = str(Counter(ctype(x)))
                C.append(rec)
        log('  C done p =', p)
    R['C_reduction'] = C
    R['C_summary'] = dict(instances=len(C), yes=sum(r['truth'] for r in C), all_agree=all(r['truth'] == r['reduction_decides'] for r in C),
                          all_subset_checks_ok=all(r['bad_subsets'] == 0 and r['all_m_nontrivial'] for r in C),
                          certificates_ok=all(r.get('certificate_checked', True) for r in C),
                          certificates_built=sum('certificate_checked' in r for r in C))
    json.dump(R, open('conj_p_results.json', 'w'), indent=1)
    # ---- D: Conj_2 ----
    log('D: Conj_2 exhaustive')
    D = []
    for N in range(1, 8 if QUICK else 10):
        D.append(verify_conj2_exhaustive(N)); log('  D done N =', N, D[-1]['pairs_tested'], D[-1]['mismatches'])
        R['D_conj2_exhaustive'] = D
        json.dump(R, open('conj_p_results.json', 'w'), indent=1)
    Drand = [verify_conj2_random(N, 60 if QUICK else 300, rng) for N in (9, 10, 11)]
    R['D_conj2_random'] = Drand
    R['D_summary'] = dict(exhaustive_pairs=sum(d['pairs_tested'] for d in D), exhaustive_mismatches=sum(d['mismatches'] for d in D),
                          random_trials=sum(d['trials'] for d in Drand), random_mismatches=sum(d['mismatches'] for d in Drand),
                          certificate_failures=sum(d['certificate_failures'] for d in D+Drand))
    # Conj_2 on the blockwise family: perfect-matching criterion vs algorithm (small sanity)
    fam = []
    for (L, k) in [(5, 3), (7, 3), (5, 4)]:
        for tr in range(3):
            m = tuple(rng.choice(units(L)) for _ in range(k))
            a, c = family(L, k, m)
            dec = conj2_decide(a, c)
            crit = criterion_decides(2, L, k, m)
            fam.append(dict(L=L, k=k, m=list(m), conj2_algorithm=dec, affine_criterion_p2=crit))
    R['D_conj2_family'] = fam
    R['D_family_agree'] = all(f['conj2_algorithm'] == f['affine_criterion_p2'] for f in fam)
    R['wall_time_s'] = round(time.time()-T0, 1)
    json.dump(R, open('conj_p_results.json', 'w'), indent=1)
    log('done; summaries:', R['A_summary'], R['B_summary'], R['C_summary'], R['D_summary'], R['D_family_agree'])

if __name__ == '__main__':
    main()
