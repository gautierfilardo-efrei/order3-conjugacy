#!/usr/bin/env python3
"""round_algebra.py -- exact structural analysis of the CPWI single-round system

    R(U',V'; A,B):   U' = U A V U,   V' = V B U V,   (U,V) in A_m x A_m.

Convention: permutations of {0..m-1} stored as image tuples, product read right
to left: (p q)(j) = p(q(j)).  Own toolkit (standard library only); the
reference solver round_attack.solve_round is called optionally (part B/C).

Parts
  K  wreath-product normal form x a x^2 = c in A_m wr C_2 (primary normal form)
  A  normal forms / invariants (identities verified exhaustively at m=5 over
     all 3600 (U,V) for many (A,B), randomly at m=6,7; untwisted cube-root
     reduction; obstruction; cycle-type invariant pruning)
  B  special cases (A=B, A=B^-1, A=id, B=id, [A,B]=1, ord 2, U'=V', [U',V']=1,
     A~B): reductions, exhaustive checks at m=5, frequencies
  C  solution statistics at m=5 (all 3600 (A,B)) and m=6 (random (A,B))

Usage:  python round_algebra.py [--m7-systems K] [--out round_algebra_results.json]
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import random
import sys
import time
from collections import Counter, defaultdict

# --------------------------------------------------------------------------
# permutation toolkit
# --------------------------------------------------------------------------

def compose(p, q):
    """p o q  (apply q first)."""
    return tuple(p[j] for j in q)


def mul(*ps):
    """Left-to-right product p1 p2 ... pk with (pq)(j) = p(q(j))."""
    r = ps[0]
    for p in ps[1:]:
        r = compose(r, p)
    return r


def inverse(p):
    out = [0] * len(p)
    for i, img in enumerate(p):
        out[img] = i
    return tuple(out)


def identity(m):
    return tuple(range(m))


def parity(p):
    """0 even, 1 odd (via cycle decomposition)."""
    seen = [False] * len(p)
    par = 0
    for i in range(len(p)):
        if not seen[i]:
            j, ln = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                ln += 1
            par ^= (ln - 1) & 1
    return par


def cycle_type(p):
    seen = [False] * len(p)
    ct = []
    for i in range(len(p)):
        if not seen[i]:
            j, ln = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                ln += 1
            ct.append(ln)
    return tuple(sorted(ct, reverse=True))


def order(p):
    o = 1
    for ln in set(cycle_type(p)):
        o = o * ln // math.gcd(o, ln)
    return o


def power(p, k):
    r = identity(len(p))
    for _ in range(k):
        r = compose(r, p)
    return r


def alternating_group(m):
    return [p for p in itertools.permutations(range(m)) if parity(p) == 0]


def random_even(m, rng):
    p = list(range(m))
    rng.shuffle(p)
    if parity(p):
        p[0], p[1] = p[1], p[0]
    return tuple(p)


def conj(g, x):
    """g x g^-1"""
    return mul(g, x, inverse(g))


def commute(p, q):
    return compose(p, q) == compose(q, p)


class GroupTable:
    """A_m with multiplication table on indices (m <= 6 is instant; m = 7 ok)."""

    def __init__(self, m):
        self.m = m
        self.elems = alternating_group(m)
        self.idx = {p: i for i, p in enumerate(self.elems)}
        n = self.n = len(self.elems)
        idx = self.idx
        self.mul = [[idx[compose(p, q)] for q in self.elems] for p in self.elems]
        self.inv = [idx[inverse(p)] for p in self.elems]
        self.e = idx[identity(m)]
        self.ctype = [cycle_type(p) for p in self.elems]
        self.order = [order(p) for p in self.elems]
        M = self.mul
        self.cube = [M[M[i][i]][i] for i in range(n)]
        self.square = [M[i][i] for i in range(n)]
        self.cube_roots = defaultdict(list)
        for i, c in enumerate(self.cube):
            self.cube_roots[c].append(i)
        # S_m-conjugacy classes intersected with A_m == same cycle type
        self.class_of_ctype = defaultdict(list)
        for i, ct in enumerate(self.ctype):
            self.class_of_ctype[ct].append(i)
        # A_m conjugacy classes (orbits under A_m conjugation)
        self.am_class = [-1] * n
        k = 0
        for i in range(n):
            if self.am_class[i] < 0:
                for g in range(n):
                    self.am_class[M[M[g][i]][self.inv[g]]] = k
                k += 1
        self.n_am_classes = k

    def prod(self, *ids):
        r = ids[0]
        M = self.mul
        for i in ids[1:]:
            r = M[r][i]
        return r

    def round(self, u, v, a, b):
        M = self.mul
        return M[M[M[u][a]][v]][u], M[M[M[v][b]][u]][v]

    def preimage_counts(self, a, b):
        """counts[up*n+vp] = number of (u,v) with round(u,v)=(up,vp)."""
        n, M = self.n, self.mul
        counts = [0] * (n * n)
        for u in range(n):
            Mu = M[u]
            ua = M[u][a]
            Mua = M[ua]
            for v in range(n):
                up = M[Mua[v]][u]
                vp = M[M[M[v][b]][u]][v]
                counts[up * n + vp] += 1
        return counts

    def solutions(self, up, vp, a, b):
        """Exhaustive solution set (list of (u,v) indices)."""
        n, M = self.n, self.mul
        # use elimination: v = a^-1 u^-1 up u^-1, then check second equation
        ai = self.inv[a]
        sols = []
        for u in range(n):
            ui = self.inv[u]
            v = M[M[M[ai][ui]][up]][ui]
            if M[M[M[v][b]][u]][v] == vp and M[M[M[u][a]][v]][u] == up:
                sols.append((u, v))
        return sols


def pmf_from_counts(counts):
    c = Counter(counts)
    tot = sum(c.values())
    return {int(k): v / tot for k, v in sorted(c.items())}


def planted_pmf_from_counts(counts):
    """Size-biased: a planted (U,V) lands on target t with prob 1/N^2; the
    number of solutions of its target is counts[t]; pmf over k weighted by k."""
    c = Counter(counts)
    tot = sum(k * v for k, v in c.items())
    return {int(k): k * v / tot for k, v in sorted(c.items()) if k > 0}


def poisson_pmf(lam, kmax):
    return {k: math.exp(-lam) * lam ** k / math.factorial(k) for k in range(kmax + 1)}


def size_biased_poisson_pmf(lam, kmax):
    # P(k) = k e^-lam lam^k / (k! lam) = e^-lam lam^(k-1)/(k-1)!, k>=1
    return {k: math.exp(-lam) * lam ** (k - 1) / math.factorial(k - 1) for k in range(1, kmax + 1)}


def tv_distance(p, q):
    keys = set(p) | set(q)
    return 0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in keys)


def moments(pmf):
    mean = sum(k * v for k, v in pmf.items())
    var = sum((k - mean) ** 2 * v for k, v in pmf.items())
    return mean, var


# --------------------------------------------------------------------------
# Part A: identities on tuples (work for any m)
# --------------------------------------------------------------------------

def check_identities(U, V, A, B, g=None, h=None):
    Ai, Bi, Ui, Vi = inverse(A), inverse(B), inverse(U), inverse(V)
    Up = mul(U, A, V, U)
    Vp = mul(V, B, U, V)
    Upi = inverse(Up)
    ok = {}
    # (i) elimination lemma and single-unknown equation (U: 3 times, U^-1: 2 times)
    ok["elimination_V"] = V == mul(Ai, Ui, Up, Ui)
    ok["elimination_U"] = U == mul(Bi, Vi, Vp, Vi)
    ok["single_unknown_eq"] = mul(B, U, Ai, Ui, Up, Ui) == mul(U, Upi, U, A, Vp)
    # (ii) symmetric form Y=AV, Z=BU
    Y, Z = mul(A, V), mul(B, U)
    ok["sym_Y"] = mul(Y, Z, Ai, Y) == mul(A, Vp)
    ok["sym_Z"] = mul(Z, Y, Bi, Z) == mul(B, Up)
    ok["UpA"] = mul(Up, A) == mul(mul(U, A), V, mul(U, A))
    ok["VpB"] = mul(Vp, B) == mul(mul(V, B), U, mul(V, B))
    # (iv) cycle-type invariant
    P, Q = mul(Up, Ui), mul(Vp, Vi)
    ok["PB_conj"] = mul(P, B) == mul(U, mul(A, Q), Ui)
    ok["QA_conj"] = mul(Q, A) == mul(V, mul(B, P), Vi)
    ok["ctype_inv"] = cycle_type(mul(Up, Ui, B)) == cycle_type(mul(A, Vp, Vi))
    # conjugation normal form: (A,B;U',V') ~ (g^-1 A g, B; U' g, g^-1 V') via (U g, g^-1 V)
    if g is not None:
        gi = inverse(g)
        A2, Uh, Vh = mul(gi, A, g), mul(U, g), mul(gi, V)
        ok["normal_form_A"] = (mul(Uh, A2, Vh, Uh) == mul(Up, g)
                               and mul(Vh, B, Uh, Vh) == mul(gi, Vp))
    if h is not None:
        hi = inverse(h)
        B2, Uh, Vh = mul(hi, B, h), mul(hi, U), mul(V, h)
        ok["normal_form_B"] = (mul(Uh, A, Vh, Uh) == mul(hi, Up)
                               and mul(Vh, B2, Uh, Vh) == mul(Vp, h))
    # swap symmetry
    ok["swap"] = (mul(V, B, U, V), mul(U, A, V, U)) == (Vp, Up)
    # twist defect: the untwisted reduction goes through iff  A Y B^-1 = Y  iff V^-1 A V = B
    ok["defect_equiv"] = (mul(A, Y, Bi) == Y) == (mul(Vi, A, V) == B)
    return ok


def check_untwisted(Y, Z):
    c1, c2 = mul(Y, Z, Y), mul(Z, Y, Z)
    W = mul(Y, Z)
    ok = {}
    ok["YZ_cubed"] = power(W, 3) == mul(c1, c2)
    ok["ZY_cubed"] = power(mul(Z, Y), 3) == mul(c2, c1)
    ok["Y_from_W"] = Y == mul(inverse(W), c1)
    ok["Z_from_W"] = Z == mul(inverse(c1), W, W)
    return ok


def check_A_id_reduction(U, V, B):
    """A = id: with X = U'^-1 U, C = V'U', D = U'^-1 B U':  X^2 C X = D."""
    m = len(U)
    I = identity(m)
    Up, Vp = mul(U, I, V, U), mul(V, B, U, V)
    X = mul(inverse(Up), U)
    C = mul(Vp, Up)
    D = mul(inverse(Up), B, Up)
    return mul(X, X, C, X) == D


def part_A(tables, rng, results, exhaustive_AB=60, random_trials=3000):
    t0 = time.time()
    A = {}
    # ---- (i),(ii),(iv): exhaustive over all (U,V) at m=5 for many (A,B)
    T5 = tables[5]
    fails = Counter()
    n_checked = 0
    ab_pairs = [(rng.randrange(T5.n), rng.randrange(T5.n)) for _ in range(exhaustive_AB)]
    ab_pairs += [(T5.e, T5.e), (T5.e, rng.randrange(T5.n)), (rng.randrange(T5.n), T5.e)]
    for a, b in ab_pairs:
        Ap, Bp = T5.elems[a], T5.elems[b]
        g, h = T5.elems[rng.randrange(T5.n)], T5.elems[rng.randrange(T5.n)]
        for Up_ in T5.elems:
            for Vp_ in T5.elems:
                ok = check_identities(Up_, Vp_, Ap, Bp, g, h)
                n_checked += 1
                for k, v in ok.items():
                    if not v:
                        fails[k] += 1
    A["identities_exhaustive_m5"] = {
        "n_AB_pairs": len(ab_pairs), "n_UV_pairs_each": T5.n * T5.n,
        "total_checks": n_checked, "failures": dict(fails),
        "all_hold": not fails}
    # random at m=6,7 (tuples)
    rand = {}
    for m in (6, 7):
        f = Counter()
        for _ in range(random_trials):
            U, V, Ap, Bp = (random_even(m, rng) for _ in range(4))
            g, h = random_even(m, rng), random_even(m, rng)
            for k, v in check_identities(U, V, Ap, Bp, g, h).items():
                if not v:
                    f[k] += 1
        rand[str(m)] = {"trials": random_trials, "failures": dict(f), "all_hold": not f}
    A["identities_random_m6_m7"] = rand
    A["single_unknown_equation"] = {
        "equation": "B U A^-1 U^-1 U' U^-1 = U U'^-1 U A V'",
        "occurrences_U": 3, "occurrences_U_inverse": 2, "total_occurrences": 5,
        "derivation": "substitute V = A^-1 U^-1 U' U^-1 into V' = V B U V, left-multiply by U U'^-1 U A"}

    # ---- (iii) untwisted system Y Z Y = c1, Z Y Z = c2 at m=5, exhaustive
    n, M = T5.n, T5.mul
    img = Counter()
    bad_cube = 0
    for y in range(n):
        for z in range(n):
            c1 = M[M[y][z]][y]
            c2 = M[M[z][y]][z]
            w = M[y][z]
            if M[M[w][w]][w] != M[c1][c2]:
                bad_cube += 1
            img[(c1, c2)] += 1
    # count of solutions == number of cube roots of c1 c2 in A_5, for ALL 3600 targets
    mismatch = 0
    n_targets_with_solution = 0
    for c1 in range(n):
        for c2 in range(n):
            pred = len(T5.cube_roots.get(M[c1][c2], []))
            if pred != img.get((c1, c2), 0):
                mismatch += 1
            if pred:
                n_targets_with_solution += 1
    # explicit solution set check for 200 random targets
    set_mismatch = 0
    for _ in range(200):
        c1, c2 = rng.randrange(n), rng.randrange(n)
        c = M[c1][c2]
        pred = set()
        for w in T5.cube_roots.get(c, []):
            y = M[T5.inv[w]][c1]
            z = M[M[T5.inv[c1]][w]][w]
            pred.add((y, z))
        actual = {(y, z) for y in range(n) for z in range(n)
                  if M[M[y][z]][y] == c1 and M[M[z][y]][z] == c2}
        if pred != actual:
            set_mismatch += 1
    # random tuple checks m=6,7
    unt_rand = {}
    for m in (6, 7):
        f = Counter()
        for _ in range(random_trials):
            Y, Z = random_even(m, rng), random_even(m, rng)
            for k, v in check_untwisted(Y, Z).items():
                if not v:
                    f[k] += 1
        unt_rand[str(m)] = {"trials": random_trials, "failures": dict(f)}
    cube_root_hist = Counter(len(T5.cube_roots.get(c, [])) for c in range(n))
    A["untwisted_reduction_m5"] = {
        "YZ_cubed_equals_c1c2_failures": bad_cube,
        "solution_count_equals_cube_root_count_mismatches_over_3600_targets": mismatch,
        "explicit_solution_set_mismatches_200_targets": set_mismatch,
        "targets_with_solution": n_targets_with_solution, "targets_total": n * n,
        "cube_root_count_histogram_A5": {str(k): v for k, v in sorted(cube_root_hist.items())},
        "solution_set": "{(Y,Z) = (W^-1 c1, c1^-1 W^2) : W in A_m, W^3 = c1 c2}",
        "consistency_condition": "c1 c2 is a cube in A_m (then #solutions = #cube roots of c1 c2; each W gives exactly one solution)",
        "random_m6_m7": unt_rand}

    # ---- (iii) obstruction for the twisted (CPWI) shape
    # F(target) := cube-root family built from c1 = A V', c2' = B U' A^-1 (pretending the twist is absent)
    # Claim:  F ∩ Sol = Sol ∩ {V^-1 A V = B} = F ∩ {V^-1 A V = B}
    obs = {"instances": 0, "claim_F_cap_Sol_eq_Sol_cap_locus": 0, "claim_F_cap_Sol_eq_F_cap_locus": 0,
           "total_solutions": 0, "solutions_on_locus": 0, "F_elements_total": 0,
           "F_elements_that_are_solutions": 0, "solutions_on_locus_U": 0, "solutions_on_both_loci": 0}
    inv = T5.inv
    ab_list = [(rng.randrange(n), rng.randrange(n)) for _ in range(20)]
    # include A=B and A~B instances (locus non-empty) and A=B=id
    for _ in range(10):
        a = rng.randrange(n)
        g = rng.randrange(n)
        ab_list.append((a, a))
        ab_list.append((a, M[M[inv[g]][a]][g]))
    ab_list.append((T5.e, T5.e))
    for a, b in ab_list:
        counts = T5.preimage_counts(a, b)
        ai, bi = inv[a], inv[b]
        for t in range(n * n):
            if counts[t] == 0 and rng.random() > 0.05:
                continue  # sample a few empty targets, all non-empty ones
            up, vp = divmod(t, n)
            sol = set(T5.solutions(up, vp, a, b))
            assert len(sol) == counts[t]
            c1 = M[a][vp]
            c2p = M[M[b][up]][ai]
            F = set()
            for w in T5.cube_roots.get(M[c1][c2p], []):
                y = M[inv[w]][c1]
                zt = M[M[inv[c1]][w]][w]
                v = M[ai][y]
                u = M[M[bi][zt]][a]
                F.add((u, v))
            locus = lambda uv: M[M[inv[uv[1]]][a]][uv[1]] == b
            sol_locus = {s for s in sol if locus(s)}
            F_locus = {s for s in F if locus(s)}
            FS = F & sol
            obs["instances"] += 1
            obs["claim_F_cap_Sol_eq_Sol_cap_locus"] += (FS == sol_locus)
            obs["claim_F_cap_Sol_eq_F_cap_locus"] += (FS == F_locus)
            obs["total_solutions"] += len(sol)
            obs["solutions_on_locus"] += len(sol_locus)
            obs["F_elements_total"] += len(F)
            obs["F_elements_that_are_solutions"] += len(FS)
            locU = {s for s in sol if M[M[inv[s[0]]][b]][s[0]] == a}
            obs["solutions_on_locus_U"] += len(locU)
            obs["solutions_on_both_loci"] += len(locU & sol_locus)
    obs["all_claims_hold"] = (obs["claim_F_cap_Sol_eq_Sol_cap_locus"] == obs["instances"]
                              == obs["claim_F_cap_Sol_eq_F_cap_locus"])
    # non-invariance: #solutions is not a function of conjugacy data of (c1, c2') / of the target
    a, b = rng.randrange(n), rng.randrange(n)
    counts = T5.preimage_counts(a, b)
    by_ct = defaultdict(set)
    by_pair_ct = defaultdict(set)
    for t in range(n * n):
        up, vp = divmod(t, n)
        c1 = M[a][vp]
        c2p = M[M[b][up]][T5.inv[a]]
        by_ct[T5.ctype[M[c1][c2p]]].add(counts[t])
        by_pair_ct[(T5.ctype[up], T5.ctype[vp], T5.ctype[M[up][vp]])].add(counts[t])
    A["obstruction"] = {
        "statement": ("With Zt = Z A^-1 the system reads Y Zt Y = c1 := A V' and Zt (A Y B^-1) Zt = c2' := B U' A^-1. "
                      "The inner letter of the second equation is A Y B^-1, not Y; the identity (Y Zt)^3 = (Y Zt Y)(Zt Y Zt) "
                      "therefore produces c1 * Zt Y Zt, which is NOT determined by the data. The defect A Y B^-1 Y^-1 vanishes "
                      "iff V^-1 A V = B. Symmetrically (Yt = Y B^-1) the defect vanishes iff U^-1 B U = A."),
        "numerical_m5": obs,
        "number_of_solutions_not_determined_by_ctype_of_c1c2p": {
            "instance_AB": [T5.elems[a], T5.elems[b]],
            "ctype_classes_with_several_solution_counts":
                sum(1 for s in by_ct.values() if len(s) > 1), "ctype_classes_total": len(by_ct),
            "example": {str(k): sorted(v) for k, v in list(by_ct.items())[:6]}},
        "number_of_solutions_not_determined_by_ctypes_of_(Up,Vp,UpVp)": {
            "classes_with_several_counts": sum(1 for s in by_pair_ct.values() if len(s) > 1),
            "classes_total": len(by_pair_ct)},
    }

    # ---- (iv) cycle-type invariant: pruning power at m=5,6
    prune = {}
    for m in (5, 6):
        T = tables[m]
        n, M, inv = T.n, T.mul, T.inv
        # exact value Σ_λ p(λ)^2 over even cycle types
        p = Counter(T.ctype)
        sum_p2 = sum((v / n) ** 2 for v in p.values())
        frac_U_pass = []
        frac_U_pass_nonsol = []
        pair_survive = []
        trials = 60 if m == 5 else 30
        for _ in range(trials):
            a, b, u, v = (rng.randrange(n) for _ in range(4))
            up, vp = T.round(u, v, a, b)
            ai = inv[a]
            lamU = [T.ctype[M[M[up][inv[x]]][b]] for x in range(n)]
            muV = [T.ctype[M[M[a][vp]][inv[x]]] for x in range(n)]
            cU, cV = Counter(lamU), Counter(muV)
            pair_survive.append(sum(cU[k] * cV.get(k, 0) for k in cU) / (n * n))
            # per-U test with V eliminated
            passing = 0
            sols = 0
            for x in range(n):
                xi = inv[x]
                vx = M[M[M[ai][xi]][up]][xi]
                if lamU[x] == T.ctype[M[M[a][vp]][inv[vx]]]:
                    passing += 1
                    if M[M[M[vx][b]][x]][vx] == vp:
                        sols += 1
            frac_U_pass.append(passing / n)
            frac_U_pass_nonsol.append((passing - sols) / n)
        prune[str(m)] = {
            "trials": trials,
            "sum_p_lambda_squared_exact": sum_p2,
            "mean_fraction_UV_pairs_surviving_ctype_match": sum(pair_survive) / trials,
            "mean_fraction_U_passing_ctype_test_after_elimination": sum(frac_U_pass) / trials,
            "mean_fraction_U_passing_but_not_solution": sum(frac_U_pass_nonsol) / trials,
            "n_even_cycle_types": len(p)}
    # asymptotic remark values for m=7 (from table) and m=8,9 (from partition enumeration)
    for m in (7, 8, 9, 10, 12, 16, 20):
        tot = 0.0
        cnt = 0
        for lam in partitions(m):
            if sum(l - 1 for l in lam) % 2 == 0:
                z = 1
                for l, mult in Counter(lam).items():
                    z *= l ** mult * math.factorial(mult)
                pl = 2.0 / z  # |class| / |A_m| = (m!/z) / (m!/2)
                tot += pl * pl
                cnt += 1
        prune[f"sum_p2_m{m}"] = tot
    A["ctype_invariant_pruning"] = prune
    A["elapsed_s"] = time.time() - t0
    results["A"] = A
    return A


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


# --------------------------------------------------------------------------
# Part B: special cases
# --------------------------------------------------------------------------

def count_stats(counts, n):
    c = Counter(counts)
    nz = sum(v for k, v in c.items() if k > 0)
    return {"targets": n * n, "targets_with_solution": nz,
            "fraction_targets_with_solution": nz / (n * n),
            "max_solutions": max(c), "planted_mean_solutions": sum(k * k * v for k, v in c.items()) / (n * n),
            "pmf_random_target": {str(k): v / (n * n) for k, v in sorted(c.items())}}


def load_solver():
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in [os.path.join(here, "repo", "cpwi-post-quantum-dynamical-inversion"),
                 os.path.join(here, "cpwi-post-quantum-dynamical-inversion"), here]:
        if os.path.exists(os.path.join(cand, "round_attack.py")):
            sys.path.insert(0, cand)
            try:
                import round_attack  # noqa
                return round_attack.solve_round
            except Exception as exc:  # pragma: no cover
                print("solver import failed:", exc, file=sys.stderr)
                return None
    return None


def solver_nodes(solve_round, T, instances, rng, budget=300_000):
    """Average nodes to enumerate ALL solutions (max_solutions large) on planted instances."""
    if solve_round is None:
        return None
    nodes = []
    for a, b in instances:
        u, v = rng.randrange(T.n), rng.randrange(T.n)
        up, vp = T.round(u, v, a, b)
        sols, st = solve_round(T.elems[up], T.elems[vp], T.elems[a], T.elems[b],
                               max_solutions=10 ** 6, node_budget=budget)
        nodes.append(st["nodes"])
    return {"instances": len(nodes), "mean_nodes": sum(nodes) / len(nodes),
            "median_nodes": sorted(nodes)[len(nodes) // 2], "max_nodes": max(nodes)}


def part_B(tables, rng, results, solve_round):
    t0 = time.time()
    T = tables[5]
    n, M, inv, e = T.n, T.mul, T.inv, T.e
    B = {}

    def ab_case_lists():
        cases = {}
        allpairs = [(a, b) for a in range(n) for b in range(n)]
        cases["generic_sample"] = [(rng.randrange(n), rng.randrange(n)) for _ in range(40)]
        cases["A_eq_B"] = [(a, a) for a in range(n)]
        cases["A_eq_Binv"] = [(a, inv[a]) for a in range(n)]
        cases["A_eq_id"] = [(e, b) for b in range(n)]
        cases["B_eq_id"] = [(a, e) for a in range(n)]
        cases["A_eq_B_eq_id"] = [(e, e)]
        comm = [(a, b) for a, b in allpairs if M[a][b] == M[b][a] and a != b and a != e and b != e]
        rng.shuffle(comm)
        cases["A_B_commute_nontrivial"] = comm[:60]
        cases["A_B_commute_all_count"] = len([1 for a, b in allpairs if M[a][b] == M[b][a]])
        cases["A_order2"] = [(a, rng.randrange(n)) for a in range(n) if T.order[a] == 2]
        cases["A_conj_B_notequal"] = []
        for a in range(n):
            g = rng.randrange(n)
            b = M[M[inv[g]][a]][g]
            if b != a:
                cases["A_conj_B_notequal"].append((a, b))
        return cases

    cases = ab_case_lists()
    B["A_B_commute_all_count_m5"] = cases.pop("A_B_commute_all_count")
    stats = {}
    for name, lst in cases.items():
        agg = Counter()
        for a, b in lst:
            agg.update(T.preimage_counts(a, b))
        st = count_stats(list(agg.elements()), n)
        st["n_AB_instances"] = len(lst)
        # normalise pmf by number of instances
        st["pmf_random_target"] = {k: v for k, v in st["pmf_random_target"].items()}
        st["planted_mean_solutions"] = st["planted_mean_solutions"] / len(lst)
        st["random_target_mean_solutions"] = 1.0  # exact: sum of counts = |A_m|^2 = #targets
        st["fraction_targets_with_solution"] = st["targets_with_solution"] / (n * n * len(lst))
        st["targets_with_solution"] = st["targets_with_solution"]
        st.pop("targets")
        sn = solver_nodes(solve_round, T, lst[:15], rng)
        if sn:
            st["solver_nodes_planted_m5"] = sn
        stats[name] = st
    B["stats_by_case_m5"] = stats

    # ---- reductions
    red = {}
    # (1) A=B=id : full cube-root reduction.  U'=UVU, V'=VUV  => (UV)^3 = U'V'
    counts = T.preimage_counts(e, e)
    mism = 0
    for t in range(n * n):
        up, vp = divmod(t, n)
        if counts[t] != len(T.cube_roots.get(M[up][vp], [])):
            mism += 1
    red["A_eq_B_eq_id"] = {"reduction": "cube roots: (U,V) = (W^-1 U', U'^-1 W^2), W^3 = U'V'",
                           "mismatches_over_3600_targets": mism, "polynomial": True}
    # (2) A=id : X = U'^-1 U satisfies X^2 C X = D with C = V'U', D = U'^-1 B U' (3 occurrences)
    mism = 0
    tested = 0
    for b in rng.sample(range(n), 8):
        counts = T.preimage_counts(e, b)
        for t in range(n * n):
            up, vp = divmod(t, n)
            C = M[vp][up]
            D = M[M[inv[up]][b]][up]
            k = sum(1 for x in range(n) if M[M[M[x][x]][C]][x] == D)
            tested += 1
            if k != counts[t]:
                mism += 1
    red["A_eq_id"] = {"reduction": "single unknown X=U'^-1 U with 3 occurrences: X^2 (V'U') X = U'^-1 B U'  (equivalently X^3 (V'U') = X (U'^-1 B U') X^-1); not a pure power unless V'U' = id",
                      "targets_tested": tested, "count_mismatches": mism,
                      "polynomial": "unknown (twisted cube root); yes when V'U'=id or B=id"}
    # (3) conjugation normal form  (A,B;U',V') ~ (g^-1 A g, B; U'g, g^-1 V')  -> A~B reduces to A=B
    mism = 0
    for _ in range(6):
        a, b, g = rng.randrange(n), rng.randrange(n), rng.randrange(n)
        c1 = T.preimage_counts(a, b)
        a2 = M[M[inv[g]][a]][g]
        c2 = T.preimage_counts(a2, b)
        for t in range(n * n):
            up, vp = divmod(t, n)
            if c1[t] != c2[M[up][g] * n + M[inv[g]][vp]]:
                mism += 1
    red["conjugation_normal_form"] = {
        "statement": "(U,V) solves R(U',V';A,B) iff (Ug, g^-1V) solves R(U'g, g^-1V'; g^-1Ag, B); iff (h^-1U, Vh) solves R(h^-1U', V'h; A, h^-1Bh). Hence the difficulty depends on (A,B) only through the conjugacy classes of A and of B (in A_m; in S_m if g,h odd, moving unknowns to the odd coset).",
        "count_mismatches_6_instances_x_3600_targets": mism,
        "consequence": "A~B reduces to A=B; A~B^-1 (true in S_m for all A) reduces to A=B^-1; A=B^-1 reduces to A=B whenever A is A_m-conjugate to A^-1"}
    # A ~_{A_m} A^-1 ?
    selfinv = sum(1 for a in range(n) if T.am_class[a] == T.am_class[inv[a]])
    red["conjugation_normal_form"]["fraction_A_in_A5_with_A_Am_conjugate_to_Ainv"] = selfinv / n
    # (4) A=B : swap symmetry, locus fraction, diagonal solutions when U'=V'
    swap_mism = 0
    locus_sol = 0
    tot_sol = 0
    diag_mism = 0
    for a in rng.sample(range(n), 10):
        counts = T.preimage_counts(a, a)
        for t in range(n * n):
            up, vp = divmod(t, n)
            if counts[t] != counts[vp * n + up]:
                swap_mism += 1
            if counts[t]:
                sols = T.solutions(up, vp, a, a)
                tot_sol += len(sols)
                locus_sol += sum(1 for u, v in sols if M[v][a] == M[a][v])
        # U'=V' : diagonal solutions U=V <-> U A U^2 = U'
        for up in range(n):
            sols = T.solutions(up, up, a, a)
            diag = {u for u, v in sols if u == v}
            pred = {u for u in range(n) if M[M[M[u][a]][u]][u] == up}
            if diag != pred:
                diag_mism += 1
    red["A_eq_B"] = {"swap_symmetry_count_mismatches": swap_mism,
                     "fraction_solutions_with_V_in_centralizer_of_A": locus_sol / max(tot_sol, 1),
                     "note": "solutions with V in C(A) are exactly the cube-root family (obstruction locus); the rest are not reachable by the reduction",
                     "diagonal_solutions_when_Up_eq_Vp_mismatches": diag_mism,
                     "polynomial": "no reduction found (swap symmetry only)"}
    # (5) frequencies of the cases for random instances at m=5,6,7
    freq = {}
    for m in (5, 6, 7):
        if m in tables:
            Tm = tables[m]
            nm = Tm.n
            pc = Counter(Tm.ctype)
            sum_p2 = sum((v / nm) ** 2 for v in pc.values())
            k_am = Tm.n_am_classes
            invol = sum(1 for o in Tm.order if o == 2)
            am_p2 = sum((v / nm) ** 2 for v in Counter(Tm.am_class).values())
        else:
            G = alternating_group(m)
            nm = len(G)
            pc = Counter(cycle_type(p) for p in G)
            sum_p2 = sum((v / nm) ** 2 for v in pc.values())
            k_am = None
            invol = sum(1 for p in G if order(p) == 2)
            am_p2 = None
        freq[str(m)] = {
            "|A_m|": nm,
            "P(A=B)": 1 / nm, "P(A=B^-1)": 1 / nm, "P(A=id)": 1 / nm, "P(B=id)": 1 / nm,
            "P(A~B in S_m)": sum_p2, "P(A~B in A_m)": am_p2,
            "P(A,B commute)=k(A_m)/|A_m|": (k_am / nm) if k_am else None,
            "P(ord A = 2)": invol / nm,
            "P(U'=V')": 1 / nm, "P(U',V' commute)": (k_am / nm) if k_am else None}
    B["case_frequencies"] = freq
    # empirical commuting probability m=5
    B["case_frequencies"]["5"]["P(A,B commute) empirical"] = B["A_B_commute_all_count_m5"] / (n * n)
    # (6) U'=V' and U',V' commuting: stats only (generic A,B)
    misc = {}
    for name, cond in (("Up_eq_Vp", lambda up, vp: up == vp),
                       ("Up_Vp_commute", lambda up, vp: M[up][vp] == M[vp][up])):
        agg = Counter()
        n_t = 0
        for a, b in cases["generic_sample"][:20]:
            counts = T.preimage_counts(a, b)
            for t in range(n * n):
                up, vp = divmod(t, n)
                if cond(up, vp):
                    agg[counts[t]] += 1
                    n_t += 1
        misc[name] = {"targets": n_t, "pmf": {str(k): v / n_t for k, v in sorted(agg.items())},
                      "mean_solutions": sum(k * v for k, v in agg.items()) / n_t,
                      "max": max(agg)}
    B["target_special_cases_m5"] = misc
    # (7) m=6 spot checks: A=id reduction, normal form, A=B=id cube roots, swap symmetry (random targets)
    T6 = tables[6]
    n6, M6, inv6, e6 = T6.n, T6.mul, T6.inv, T6.e
    chk = Counter()
    for _ in range(300):
        up, vp, b, g = (rng.randrange(n6) for _ in range(4))
        cnt = len(T6.solutions(up, vp, e6, b))
        C_, D_ = M6[vp][up], M6[M6[inv6[up]][b]][up]
        k = sum(1 for x in range(n6) if M6[M6[M6[x][x]][C_]][x] == D_)
        chk["A_eq_id_mismatch"] += (k != cnt)
        a = rng.randrange(n6)
        c1 = len(T6.solutions(up, vp, a, b))
        c2 = len(T6.solutions(M6[up][g], M6[inv6[g]][vp], M6[M6[inv6[g]][a]][g], b))
        chk["normal_form_mismatch"] += (c1 != c2)
        chk["cube_root_mismatch"] += (len(T6.solutions(up, vp, e6, e6)) != len(T6.cube_roots.get(M6[up][vp], [])))
        chk["swap_mismatch"] += (len(T6.solutions(up, vp, a, a)) != len(T6.solutions(vp, up, a, a)))
        chk["trials"] += 1
    B["m6_spot_checks"] = dict(chk)
    # (8) reference-solver node counts (full enumeration) per case at m=6 and m=7
    if solve_round is not None:
        comp = {}
        for m in (6, 7):
            row = {}
            for name, gen in (("generic", lambda: (random_even(m, rng), random_even(m, rng))),
                              ("A_eq_B", lambda: (lambda a: (a, a))(random_even(m, rng))),
                              ("A_eq_id", lambda: (identity(m), random_even(m, rng))),
                              ("A_eq_B_eq_id", lambda: (identity(m), identity(m))),
                              ("A_B_commute", lambda: (lambda a: (a, power(a, rng.randrange(1, order(a) + 1))))(random_even(m, rng)))):
                nodes = []
                for _ in range(25):
                    A_, B_ = gen()
                    U, V = random_even(m, rng), random_even(m, rng)
                    Up, Vp = mul(U, A_, V, U), mul(V, B_, U, V)
                    _, st = solve_round(Up, Vp, A_, B_, max_solutions=10 ** 6, node_budget=2_000_000)
                    nodes.append(st["nodes"])
                nodes.sort()
                row[name] = {"instances": len(nodes), "mean_nodes": sum(nodes) / len(nodes),
                             "median_nodes": nodes[len(nodes) // 2], "max_nodes": nodes[-1]}
            comp[str(m)] = row
        B["solver_nodes_by_case_m6_m7"] = comp
    B["reductions"] = red
    B["elapsed_s"] = time.time() - t0
    results["B"] = B
    return B


# --------------------------------------------------------------------------
# Part C: solution statistics
# --------------------------------------------------------------------------

def part_C(tables, rng, results, m6_systems=10, m7_systems=0, solve_round=None):
    t0 = time.time()
    C = {}
    T = tables[5]
    n = T.n
    agg = Counter()
    max_by_ab = []
    for a in range(n):
        for b in range(n):
            counts = T.preimage_counts(a, b)
            agg.update(counts)
            max_by_ab.append(max(counts))
    pmf_r = {k: v / sum(agg.values()) for k, v in sorted(agg.items())}
    tot_pl = sum(k * v for k, v in agg.items())
    pmf_p = {k: k * v / tot_pl for k, v in sorted(agg.items()) if k > 0}
    kmax = max(agg)
    P1 = poisson_pmf(1.0, kmax)
    SB = size_biased_poisson_pmf(1.0, kmax)
    C["m5_all_3600_AB"] = {
        "systems": n * n, "targets_per_system": n * n,
        "pmf_random_target": {str(k): v for k, v in pmf_r.items()},
        "pmf_planted_target": {str(k): v for k, v in pmf_p.items()},
        "poisson1_pmf": {str(k): v for k, v in P1.items()},
        "size_biased_poisson1_pmf": {str(k): v for k, v in SB.items()},
        "random_mean_var": moments(pmf_r), "planted_mean_var": moments(pmf_p),
        "poisson1_mean_var": (1.0, 1.0), "size_biased_mean_var": (2.0, 1.0),
        "tv_random_vs_poisson1": tv_distance(pmf_r, P1),
        "tv_planted_vs_size_biased_poisson1": tv_distance(pmf_p, SB),
        "max_solutions_observed": kmax,
        "max_over_AB_histogram": {str(k): v for k, v in sorted(Counter(max_by_ab).items())},
        "fraction_targets_with_solution": 1 - pmf_r.get(0, 0.0),
        "fraction_planted_with_unique_solution": pmf_p.get(1, 0.0)}
    # which (A,B) attain the max?
    worst = []
    for a in range(n):
        for b in range(n):
            pass
    C["m5_all_3600_AB"]["note_max"] = "see m5_max_instances"
    # locate a few max instances
    mx = []
    for a in range(0, n, 7):
        for b in range(0, n, 11):
            counts = T.preimage_counts(a, b)
            k = max(counts)
            if k >= kmax - 2:
                t = counts.index(k)
                mx.append({"A": T.elems[a], "B": T.elems[b], "Up": T.elems[t // n], "Vp": T.elems[t % n], "solutions": k})
    C["m5_max_instances_sample"] = mx[:5]

    # m=6
    T6 = tables[6]
    n6 = T6.n
    agg6 = Counter()
    per_sys = []
    for _ in range(m6_systems):
        a, b = rng.randrange(n6), rng.randrange(n6)
        counts = T6.preimage_counts(a, b)
        c = Counter(counts)
        agg6.update(c)
        per_sys.append({"A": T6.elems[a], "B": T6.elems[b], "max": max(c),
                        "frac_with_solution": 1 - c[0] / (n6 * n6),
                        "planted_mean": sum(k * k * v for k, v in c.items()) / (n6 * n6)})
    pmf_r6 = {k: v / sum(agg6.values()) for k, v in sorted(agg6.items())}
    tot6 = sum(k * v for k, v in agg6.items())
    pmf_p6 = {k: k * v / tot6 for k, v in sorted(agg6.items()) if k > 0}
    k6 = max(agg6)
    C["m6_random_AB"] = {
        "systems": m6_systems, "targets_per_system": n6 * n6,
        "pmf_random_target": {str(k): v for k, v in pmf_r6.items()},
        "pmf_planted_target": {str(k): v for k, v in pmf_p6.items()},
        "random_mean_var": moments(pmf_r6), "planted_mean_var": moments(pmf_p6),
        "tv_random_vs_poisson1": tv_distance(pmf_r6, poisson_pmf(1.0, k6)),
        "tv_planted_vs_size_biased_poisson1": tv_distance(pmf_p6, size_biased_poisson_pmf(1.0, k6)),
        "max_solutions_observed": k6, "per_system": per_sys,
        "fraction_planted_with_unique_solution": pmf_p6.get(1, 0.0)}
    # m=7 exhaustive (table) if available
    if 7 in tables:
        T7 = tables[7]
        n7 = T7.n
        agg7 = Counter()
        per7 = []
        for _ in range(m7_systems):
            a, b = rng.randrange(n7), rng.randrange(n7)
            counts = T7.preimage_counts(a, b)
            c = Counter(counts)
            agg7.update(c)
            per7.append({"max": max(c), "frac_with_solution": 1 - c[0] / (n7 * n7),
                         "planted_mean": sum(k * k * v for k, v in c.items()) / (n7 * n7)})
        pmf_r7 = {k: v / sum(agg7.values()) for k, v in sorted(agg7.items())}
        tot7 = sum(k * v for k, v in agg7.items())
        pmf_p7 = {k: k * v / tot7 for k, v in sorted(agg7.items()) if k > 0}
        k7 = max(agg7)
        C["m7_random_AB_exhaustive"] = {
            "systems": m7_systems, "targets_per_system": n7 * n7,
            "pmf_random_target": {str(k): v for k, v in pmf_r7.items()},
            "pmf_planted_target": {str(k): v for k, v in pmf_p7.items()},
            "random_mean_var": moments(pmf_r7), "planted_mean_var": moments(pmf_p7),
            "tv_random_vs_poisson1": tv_distance(pmf_r7, poisson_pmf(1.0, k7)),
            "tv_planted_vs_size_biased_poisson1": tv_distance(pmf_p7, size_biased_poisson_pmf(1.0, k7)),
            "max_solutions_observed": k7, "per_system": per7,
            "fraction_planted_with_unique_solution": pmf_p7.get(1, 0.0)}
    # m=7 random planted instances via the reference solver (exhaustive enumeration per instance)
    if m7_systems and solve_round is not None:
        cnts = []
        nodes = []
        for _ in range(m7_systems):
            U, V, A, Bp = (random_even(7, rng) for _ in range(4))
            Up, Vp = mul(U, A, V, U), mul(V, Bp, U, V)
            sols, st = solve_round(Up, Vp, A, Bp, max_solutions=10 ** 6, node_budget=2_000_000)
            cnts.append(len(sols) if not st["exhausted"] else -1)
            nodes.append(st["nodes"])
        C["m7_planted_via_solver"] = {"systems": m7_systems,
                                     "solution_count_histogram": {str(k): v for k, v in sorted(Counter(cnts).items())},
                                     "mean_nodes": sum(nodes) / len(nodes), "note": "-1 = budget exhausted"}
    C["elapsed_s"] = time.time() - t0
    results["C"] = C
    return C



# --------------------------------------------------------------------------
# Part K: wreath-product normal form  x a x^2 = c  in K = A_m wr C_2 < S_{2m}
# --------------------------------------------------------------------------

def embed(U, V):
    """(U,V) in A_m x A_m -> element of S_{2m} acting by U on {0..m-1}, V on {m..2m-1}."""
    m = len(U)
    return tuple(U) + tuple(v + m for v in V)


def sigma(m):
    return tuple((j + m) % (2 * m) for j in range(2 * m))


def wreath(U, V):
    """x = (U,V)·sigma  (apply sigma first, then (U,V))."""
    return compose(embed(U, V), sigma(len(U)))


def unwreath(x):
    """inverse of wreath: x = (U,V)·sigma  ->  (U,V); returns None if x is not in the sigma-coset of the base group."""
    m2 = len(x)
    m = m2 // 2
    if any(x[j] < m for j in range(m)) or any(x[j] >= m for j in range(m, m2)):
        return None
    V = tuple(x[j] - m for j in range(m))
    U = tuple(x[j] for j in range(m, m2))
    return U, V


def check_wreath_identity(U, V, A, B):
    m = len(U)
    x = wreath(U, V)
    a = embed(B, A)
    Up, Vp = mul(U, A, V, U), mul(V, B, U, V)
    c = wreath(Up, Vp)
    ok = {}
    ok["x_a_x2_eq_c"] = mul(x, a, x, x) == c
    # x^3 = (UVU, VUV)·sigma ;  a x^3 = x^-1 c x  ;  ctype(a x^3) = ctype(c)
    x3 = mul(x, x, x)
    ok["x3_form"] = x3 == wreath(mul(U, V, U), mul(V, U, V))
    ok["a_x3_eq_conj_c"] = mul(a, x3) == mul(inverse(x), c, x)
    ok["ctype_a_x3_eq_ctype_c"] = cycle_type(mul(a, x3)) == cycle_type(c)
    # cycle type of (P,Q)·sigma in S_{2m} is the doubled cycle type of PQ
    P, Q = random_even(m, random), random_even(m, random)
    ok["ctype_wreath_doubles_PQ"] = cycle_type(wreath(P, Q)) == tuple(sorted((2 * l for l in cycle_type(mul(P, Q))), reverse=True))
    # base-group conjugation g=(g1,g2): x -> g x g^-1 realises (U,V) -> (g1 U g2^-1, g2 V g1^-1), a -> (g1 B g1^-1, g2 A g2^-1) (Lemma A4)
    g1, g2 = random_even(m, random), random_even(m, random)
    g = embed(g1, g2)
    ok["base_conjugation_is_normal_form_A4"] = (mul(g, x, inverse(g)) == wreath(mul(g1, U, inverse(g2)), mul(g2, V, inverse(g1)))
                                                 and mul(g, a, inverse(g)) == embed(mul(g1, B, inverse(g1)), mul(g2, A, inverse(g2))))
    # a = 1  <=>  A = B = id : x^3 = c
    return ok


def part_K(tables, rng, results, exhaustive_AB=20, random_trials=3000, s2m_instances=3):
    t0 = time.time()
    K = {}
    T5 = tables[5]
    n = T5.n
    # exhaustive over all 3600 (U,V) at m=5 for exhaustive_AB random (A,B) (+ (id,id))
    fails = Counter()
    checks = 0
    ab = [(rng.randrange(n), rng.randrange(n)) for _ in range(exhaustive_AB)] + [(T5.e, T5.e)]
    for a_, b_ in ab:
        for U in T5.elems:
            for V in T5.elems:
                for k, v in check_wreath_identity(U, V, T5.elems[a_], T5.elems[b_]).items():
                    checks += 1
                    if not v:
                        fails[k] += 1
    K["identity_exhaustive_m5"] = {"n_AB": len(ab), "UV_pairs_each": n * n, "checks": checks, "failures": dict(fails)}
    rand = {}
    for m in (6, 7):
        f = Counter()
        for _ in range(random_trials):
            U, V, A, B = (random_even(m, rng) for _ in range(4))
            for k, v in check_wreath_identity(U, V, A, B).items():
                if not v:
                    f[k] += 1
        rand[str(m)] = {"trials": random_trials, "failures": dict(f)}
    K["identity_random_m6_m7"] = rand
    # converse at m=5: for ALL 3600 targets of 2 systems, the solutions x of x a x^2 = c in K (base group and sigma-coset)
    # coincide with wreath(round solutions); base-group elements never solve (sigma-exponent parity)
    conv = {"systems": 0, "targets": 0, "mismatch": 0, "base_solutions": 0}
    for a_, b_ in ab[:2]:
        A, B = T5.elems[a_], T5.elems[b_]
        a = embed(B, A)
        coset = {wreath(U, V): (U, V) for U in T5.elems for V in T5.elems}
        base = [embed(U, V) for U in T5.elems for V in T5.elems]
        # x a x^2 for every coset element -> image map
        img = defaultdict(set)
        for x, (U, V) in coset.items():
            img[mul(x, a, x, x)].add((U, V))
        counts = T5.preimage_counts(a_, b_)
        conv["systems"] += 1
        for t in range(n * n):
            up, vp = divmod(t, n)
            c = wreath(T5.elems[up], T5.elems[vp])
            sol_round = set((T5.elems[u], T5.elems[v]) for u, v in T5.solutions(up, vp, a_, b_))
            conv["targets"] += 1
            if img.get(c, set()) != sol_round or len(sol_round) != counts[t]:
                conv["mismatch"] += 1
        # base group: check 50 random targets exhaustively over the 3600 base elements
        for _ in range(50):
            c = wreath(T5.elems[rng.randrange(n)], T5.elems[rng.randrange(n)])
            conv["base_solutions"] += sum(1 for x in base if mul(x, a, x, x) == c)
    K["converse_m5"] = conv
    # solutions in the whole S_{2m} (m=5: 10! = 3,628,800 candidates) for planted instances
    s2m = []
    for _ in range(s2m_instances):
        U, V, A, B = (T5.elems[rng.randrange(n)] for _ in range(4))
        a = embed(B, A)
        c = wreath(mul(U, A, V, U), mul(V, B, U, V))
        in_K = in_SmwrC2 = outside = 0
        for x in itertools.permutations(range(10)):
            if compose(compose(compose(x, a), x), x) == c:
                uv = unwreath(x)
                if uv is None:
                    # in S_m wr C_2 \ K ?  (blocks preserved or swapped, with odd components)
                    blocks_ok = (all(x[j] >= 5 for j in range(5)) and all(x[j] < 5 for j in range(5, 10))) or \
                                (all(x[j] < 5 for j in range(5)) and all(x[j] >= 5 for j in range(5, 10)))
                    if blocks_ok:
                        in_SmwrC2 += 1
                    else:
                        outside += 1
                else:
                    in_K += 1
        s2m.append({"solutions_in_K_sigma_coset": in_K, "solutions_in_SmwrC2_not_K": in_SmwrC2,
                    "solutions_outside_SmwrC2": outside, "round_solutions": len(T5.solutions(T5.idx[mul(U, A, V, U)], T5.idx[mul(V, B, U, V)], T5.idx[A], T5.idx[B]))})
    K["solutions_in_full_S2m_m5_planted"] = s2m
    # pruning power of ctype(a x^3) = ctype(c) over the sigma-coset (pairs (U,V)) at m=5,6
    prune = {}
    for m in (5, 6):
        T = tables[m]
        nm, M = T.n, T.mul
        fr = []
        for _ in range(40 if m == 5 else 15):
            a_, b_, u, v = (rng.randrange(nm) for _ in range(4))
            up, vp = T.round(u, v, a_, b_)
            target_ct = T.ctype[M[up][vp]]  # ctype(c) <-> ctype(U'V')
            # ctype(a x^3) <-> ctype(B UVU A VUV)
            passing = 0
            for x in range(nm):
                for y in range(nm):
                    if T.ctype[T.prod(b_, x, y, x, a_, y, x, y)] == target_ct:
                        passing += 1
            fr.append(passing / (nm * nm))
        prune[str(m)] = {"trials": len(fr), "mean_fraction_UV_passing": sum(fr) / len(fr)}
    K["ctype_a_x3_pruning"] = prune
    K["statement"] = ("K = (A_m x A_m) ⋊ C_2 = A_m wr C_2 < S_{2m}; (U,V) acts by U on {0..m-1}, V on {m..2m-1}; sigma: j -> j+m mod 2m. "
                      "x = (U,V)·sigma, a = (B,A), c = (U',V')·sigma.  Then x a x^2 = c  <=>  (U'=UAVU and V'=VBUV). "
                      "Every solution in K lies in the sigma-coset (sigma-exponent of x a x^2 is that of x). "
                      "a = 1 <=> A = B = id: cube root. Also a x^3 = x^-1 c x, hence ctype(a x^3) = ctype(c) and x^3 = (UVU, VUV)·sigma.")
    K["elapsed_s"] = time.time() - t0
    results["K"] = K
    return K


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="round_algebra_results.json")
    ap.add_argument("--seed", type=int, default=20260918)
    ap.add_argument("--m6-systems", type=int, default=50)
    ap.add_argument("--m7-systems", type=int, default=20)
    ap.add_argument("--parts", default="AKBC")
    ap.add_argument("--m7-table", action="store_true", help="build the A_7 table (2520 elements) for exhaustive m=7 statistics")
    args = ap.parse_args()
    rng = random.Random(args.seed)
    results = {"meta": {"seed": args.seed, "convention": "(pq)(j)=p(q(j)); round U'=UAVU, V'=VBUV",
                        "python": sys.version.split()[0]}}
    if os.path.exists(args.out):
        try:
            results = json.load(open(args.out))
        except Exception:
            pass
    t0 = time.time()
    tables = {5: GroupTable(5), 6: GroupTable(6)}
    if args.m7_table:
        tables[7] = GroupTable(7)
    results["meta"]["table_build_s"] = time.time() - t0
    solve_round = load_solver()
    results["meta"]["reference_solver_available"] = solve_round is not None

    def dump():
        json.dump(results, open(args.out, "w"), indent=1, default=str)

    if "A" in args.parts:
        part_A(tables, rng, results)
        dump()
        print("A done", round(results["A"]["elapsed_s"], 1), "s", flush=True)
    if "K" in args.parts:
        part_K(tables, rng, results)
        dump()
        print("K done", round(results["K"]["elapsed_s"], 1), "s", flush=True)
    if "B" in args.parts:
        part_B(tables, rng, results, solve_round)
        dump()
        print("B done", round(results["B"]["elapsed_s"], 1), "s", flush=True)
    if "C" in args.parts:
        part_C(tables, rng, results, args.m6_systems, args.m7_systems, solve_round)
        dump()
        print("C done", round(results["C"]["elapsed_s"], 1), "s", flush=True)


if __name__ == "__main__":
    main()
