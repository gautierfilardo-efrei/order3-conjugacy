"""Padding track: polynomial reduction Conj_3 -> P3 ?   (P3: x a x^2 = c in Sym(N))

Conventions: permutations are image tuples on {0..N-1}, (pq)(j) = p(q(j)); the word x a x x applied to j is
x(a(x(x(j)))).  Standard library only (word_solver.py is the project's propagation solver).

Contents
  * exact tools: all P3 solutions (brute force N <= 9, propagation solver beyond, completeness flag checked),
    order-3 conjugators (Conj_3) by brute force or by the block lemma on the affine family;
  * "layered" constructions  Omega x Z_M  with a' = (a^{e_r} on layer r, then layer permutation rho),
    c' = (c^{f_r}, rho') -- the non-disjoint paddings compatible with a planted lift of an order-3 conjugator
    (Lemma P1 in padding_notes.md shows disjoint paddings are inert);
  * the affine family (a = k L-cycles, c = a^{m_i} blockwise) and the exact description of its P3 solutions;
  * the repaired gadget B (second colour class absorbing the t-fixed points) of the gadget_b track.
Run:  python padding_track.py  [quick]   -> writes padding_results.json
"""
from __future__ import annotations
import itertools, json, sys, time
from word_solver import compose, inverse, mul, evaluate, solve_word, cycle_type

def identity(n): return tuple(range(n))
def power(p, e):
    n = len(p); r = identity(n)
    if e < 0: p = inverse(p); e = -e
    for _ in range(e): r = compose(p, r)
    return r
def order(p):
    from math import lcm
    o = 1
    for L in cycle_type(p): o = lcm(o, L)
    return o

# ----------------------------------------------------------------------------------------------------------------
# exact solvers
# ----------------------------------------------------------------------------------------------------------------
def p3_solutions(a, c, node_budget=5 * 10**6, max_solutions=10**6):
    """All x with x a x^2 = c.  Returns (sols, info); info['complete'] must be True for any claim."""
    N = len(c)
    if N <= 9:
        sols = [x for x in itertools.permutations(range(N)) if evaluate(['x', a, 'x', 'x'], x) == c]
        return sols, {"method": "brute", "complete": True, "N": N}
    sols, st = solve_word(['x', a, 'x', 'x'], c, max_solutions=max_solutions, node_budget=node_budget)
    return sols, {"method": "propagation", "complete": (not st["exhausted"]) and len(sols) < max_solutions,
                  "nodes": st["nodes"], "N": N}

def conj3_solutions(a, c):
    """All x with x^3 = 1 and x a x^{-1} = c (brute force; n <= 8)."""
    n = len(a)
    if cycle_type(a) != cycle_type(c): return []
    out = []
    for x in itertools.permutations(range(n)):
        if mul(x, x, x) == identity(n) and mul(x, a, inverse(x)) == c: out.append(x)
    return out

# ----------------------------------------------------------------------------------------------------------------
# affine family  a = k L-cycles,  c = a^{m_i} on block i
# ----------------------------------------------------------------------------------------------------------------
def affine_instance(L, ms):
    k = len(ms); N = k * L
    a = [0] * N; c = [0] * N
    for i, m in enumerate(ms):
        for z in range(L):
            a[i * L + z] = i * L + (z + 1) % L
            c[i * L + z] = i * L + (z + m) % L
    return tuple(a), tuple(c)

def conj3_affine(L, ms):
    """Block lemma (Theorem B): order-3 conjugators <-> sigma with sigma^3 = 1, m_i m_j m_l = 1 on 3-cycles,
    m_i^3 = 1 on fixed blocks.  Returns the admissible sigma's (each gives L^k conjugators)."""
    k = len(ms); out = []
    for sigma in itertools.permutations(range(k)):
        ok = all(sigma[sigma[sigma[i]]] == i for i in range(k))
        for i in range(k):
            if not ok: break
            if sigma[i] == i: ok = pow(ms[i], 3, L) == 1
            else:
                j, l = sigma[i], sigma[sigma[i]]
                ok = (ms[i] * ms[j] * ms[l]) % L == 1
        if ok: out.append(sigma)
    return out

def affine_decompose(x, L, k):
    """If x maps every block onto a block by an affine map z -> m z + t (mod L), return (sigma, slopes, shifts); else None."""
    sigma = [None] * k; slopes = [None] * k; shifts = [None] * k
    for i in range(k):
        imgs = [x[i * L + z] for z in range(L)]
        blocks = {v // L for v in imgs}
        if len(blocks) != 1: return None
        j = blocks.pop(); vals = [v - j * L for v in imgs]
        t = vals[0]; m = (vals[1] - vals[0]) % L
        if any(vals[z] != (m * z + t) % L for z in range(L)): return None
        sigma[i] = j; slopes[i] = m; shifts[i] = t
    return tuple(sigma), tuple(slopes), tuple(shifts)

def block_shape(x, L, k):
    """Does x map blocks onto blocks (as sets)?  Returns the block permutation or None."""
    sigma = []
    for i in range(k):
        bl = {x[i * L + z] // L for z in range(L)}
        if len(bl) != 1: return None
        sigma.append(bl.pop())
    return tuple(sigma)

# ----------------------------------------------------------------------------------------------------------------
# layered constructions on Omega x Z_M   (point (j, r) -> j + r n)
# ----------------------------------------------------------------------------------------------------------------
def layered(a, c, exps_a, exps_c, rho_a, rho_c):
    """a'(j, r) = (a^{exps_a[r]} j, rho_a[r]),  c'(j, r) = (c^{exps_c[r]} j, rho_c[r])."""
    n = len(a); M = len(exps_a)
    pa = [power(a, e) for e in exps_a]; pc = [power(c, e) for e in exps_c]
    A = [0] * (n * M); C = [0] * (n * M)
    for r in range(M):
        for j in range(n):
            A[j + r * n] = pa[r][j] + rho_a[r] * n
            C[j + r * n] = pc[r][j] + rho_c[r] * n
    return tuple(A), tuple(C)

def lift(y, n, M, pi):
    """x'(j, r) = (y j, pi[r])."""
    X = [0] * (n * M)
    for r in range(M):
        for j in range(n): X[j + r * n] = y[j] + pi[r] * n
    return tuple(X)

def layer_shape(x, n, M):
    """Return the layer map r -> layer(x(j, r)) if x maps layers onto layers, else None."""
    pi = []
    for r in range(M):
        imgs = {x[j + r * n] // n for j in range(n)}
        if len(imgs) != 1: return None
        pi.append(imgs.pop())
    return tuple(pi)

CONSTRUCTIONS = {
    # name: (exps_a, exps_c, rho_a, rho_c, intended layer map pi of the lift x' = (y, pi))
    "M2_swap_e12":   ((1, 2), (1, 2), (1, 0), (1, 0), (0, 1)),      # diagonal system y a y^2=c, y a^2 y^2=c^2 (Lemma P2)
    "M2_swap_e1m1":  ((1, -1), (1, -1), (1, 0), (1, 0), (0, 1)),
    "M2_swap_e11":   ((1, 1), (1, 1), (1, 0), (1, 0), (0, 1)),      # control: a (x) swap, every P3 solution lifts
    "M2_swap_e10":   ((1, 0), (1, 0), (1, 0), (1, 0), (0, 1)),      # one content layer + one transfer layer
    "M2_id_e12_shift": ((1, 2), (2, 1), (0, 1), (0, 1), (1, 0)),    # layer-preserving a', c'; intended x' swaps layers
    "M3_shift_e120": ((1, 2, 0), (0, 1, 2), (1, 2, 0), (1, 2, 0), (1, 2, 0)),
    "M3_shift_e121": ((1, 2, 1), (1, 1, 2), (1, 2, 0), (1, 2, 0), (1, 2, 0)),
    "M3_id_e120":    ((1, 2, 0), (1, 2, 0), (0, 1, 2), (0, 1, 2), (1, 2, 0)),
}

def check_construction(name, a, c, node_budget=5 * 10**6, conj=None):
    exps_a, exps_c, rho_a, rho_c, pi = CONSTRUCTIONS[name]
    n = len(a); M = len(exps_a)
    A, C = layered(a, c, exps_a, exps_c, rho_a, rho_c)
    if conj is None and n <= 7: conj = conj3_solutions(a, c)
    sols, info = p3_solutions(A, C, node_budget=node_budget)
    lifted = {lift(y, n, M, pi) for y in (conj or [])}
    S = set(sols)
    rec = {"construction": name, "n": n, "N": n * M, "ctype_a": cycle_type(a), "ctype_c": cycle_type(c),
           "conj3": len(conj) if conj is not None else None, "p3_padded": len(S), "complete": info["complete"],
           "lift_ok": lifted <= S, "n_lifted": len(lifted), "spurious": len(S - lifted)}
    if conj is not None: rec["equivalent"] = (len(conj) > 0) == (len(S) > 0)
    shapes = {}
    for x in S - lifted:
        sh = layer_shape(x, n, M)
        key = ("layered " + str(sh)) if sh is not None else "non-layered"
        key += " ord3" if mul(x, x, x) == identity(n * M) else " ord>3"
        shapes[key] = shapes.get(key, 0) + 1
    rec["spurious_shapes"] = shapes
    return rec, sols, A, C

def all_pairs(n, conjugate_only=True):
    """a = one representative per cycle type, c over the whole class of a (or all of S_n)."""
    reps = {}
    for p in itertools.permutations(range(n)):
        t = cycle_type(p)
        if t not in reps: reps[t] = p
    for t, a in sorted(reps.items()):
        for c in itertools.permutations(range(n)):
            if conjugate_only and cycle_type(c) != t: continue
            yield a, c

def grid_test(name, n, conjugate_only=True):
    t0 = time.time(); recs = []; bad = 0
    for a, c in all_pairs(n, conjugate_only):
        rec, sols, A, C = check_construction(name, a, c)
        rec["a"] = list(a); rec["c"] = list(c)
        recs.append(rec)
        if not rec["equivalent"] or not rec["lift_ok"]: bad += 1
    agg = {"name": name, "n": n, "pairs": len(recs), "non_equivalent": bad,
           "sum_spurious": sum(r["spurious"] for r in recs), "all_complete": all(r["complete"] for r in recs),
           "lift_ok_all": all(r["lift_ok"] for r in recs), "time_s": round(time.time() - t0, 1)}
    return agg, recs

# ----------------------------------------------------------------------------------------------------------------
# gadget B of the gadget_b track + "absorber" repair (second colour class of blocks)
# ----------------------------------------------------------------------------------------------------------------
def gadget_b(labels, L, absorbers=0, Lp=3):
    """Gadget B1 blocks (c-cycle of length L+1, chord (g0, g_{L-d})) plus `absorbers` blocks of a second colour:
    c-cycle of length Lp+1 with the chord on c-adjacent points (label 0)."""
    n = len(labels); N = n * (L + 1) + absorbers * (Lp + 1); a = list(range(N)); c = list(range(N)); p = 0
    for d in labels:
        g = list(range(p, p + L + 1)); p += L + 1
        for i in range(L + 1): c[g[i]] = g[(i + 1) % (L + 1)]
        D = L - (d % L); a[g[0]], a[g[D]] = g[D], g[0]
    for _ in range(absorbers):
        g = list(range(p, p + Lp + 1)); p += Lp + 1
        for i in range(Lp + 1): c[g[i]] = g[(i + 1) % (Lp + 1)]
        a[g[0]], a[g[Lp]] = g[Lp], g[0]
    return tuple(a), tuple(c)

if __name__ == "__main__":
    # minimal reproduction of the key numbers (see padding_notes.md); full runs were done interactively
    CONSTRUCTIONS["M4_cyc_e1233"] = ((1, 2, 3, 3), (1, 2, 3, 3), (1, 2, 3, 0), (1, 2, 3, 0), (0, 1, 2, 3))
    for m in [2, 3]:
        a, c = affine_instance(5, [m])
        rec, sols, A, C = check_construction("M4_cyc_e1233", a, c, node_budget=2 * 10**6)
        print("M4_cyc_e1233 L=5 m=%d: conj3=%s p3_padded=%d complete=%s" % (m, rec["conj3"], rec["p3_padded"], rec["complete"]))
    for (L, ms) in [(5, [2, 3]), (5, [2, 2, 4])]:
        a, c = affine_instance(L, ms); sols, info = p3_solutions(a, c)
        print(L, ms, len(sols), "all affine:", all(affine_decompose(x, L, len(ms)) is not None for x in sols), info["complete"])
