"""round_solver_study.py — profiling and variants of the propagation solver for the
CPWI single-round equation in wreath normal form,  x a x^2 = c  in  A_m wr C_2 <= S_{2m}.

Conventions: permutations are image tuples, (pq)(j) = p(q(j)), the word x a x x is applied
right to left: chain of point j is  j -> x(j) -> x(x(j)) -> a(x^2(j)) -> x(a x^2(j)) = c(j).
Block structure is enforced: x maps block 0 = {0..m-1} onto block 1 = {m..2m-1} and back.
Standard library only (matplotlib/z3 optional, imported lazily in the stages that need them).

Seeds: SHA-256 of "CPWI-SOLV-instance|m|t"  (same planted instance for every variant).
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, random, statistics, sys, time
from collections import Counter

NA = -1

# ----------------------------------------------------------------------------- primitives
def compose(p, q): return tuple(p[q[j]] for j in range(len(p)))
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
def parity(p): return sum(L - 1 for L in cycle_type(p)) % 2
def random_even(m, rng):
    while True:
        p = list(range(m)); rng.shuffle(p); p = tuple(p)
        if parity(p) == 0: return p
def seed_from(label, *args):
    h = hashlib.sha256("|".join([label] + [str(a) for a in args]).encode()).hexdigest()
    return int(h[:16], 16)

def wreath_instance(m, t):
    """Planted round instance in wreath form. Returns dict with a, c, x, U, V, A, B, up, vp."""
    rng = random.Random(seed_from("CPWI-SOLV-instance", m, t))
    U, V, A, B = [random_even(m, rng) for _ in range(4)]
    emb = lambda P, Q: tuple(list(P) + [q + m for q in Q])
    sig = tuple((j + m) % (2 * m) for j in range(2 * m))
    x = compose(emb(U, V), sig); a = emb(B, A)
    c = mul(x, a, x, x)
    up = tuple(c[k + m] for k in range(m)); vp = tuple(c[j] - m for j in range(m))
    return dict(m=m, t=t, a=a, c=c, x=x, U=U, V=V, A=A, B=B, up=up, vp=vp)

def check_solution(a, c, x):
    return mul(x, a, x, x) == c

# ----------------------------------------------------------------------------- the solver
class WordSolver:
    """Iterative DFS with trail-based undo for x a x^2 = c, block structure enforced.

    heuristic: 'first'  — first free point (baseline, = word_solver.solve_word order)
               'mcf'    — most-constrained point first, domain measured by forward checking
                          (assign + propagate for every candidate; singleton domains forced)
               'chains' — point occurring in the most partially-known chains
    ctype:     apply the cycle-type relaxation ctype(a x^3) = ctype(c) on closed cycles and
               open segments of the partial a x^3 after every propagation.
    """
    def __init__(self, a, c, m, heuristic="first", ctype=False, node_budget=3_000_000,
                 max_solutions=1, planted=None, profile=False):
        self.a, self.c, self.m, self.N = a, c, m, len(c)
        self.ai = inverse(a)
        self.heuristic, self.use_ctype = heuristic, ctype
        self.node_budget, self.max_solutions = node_budget, max_solutions
        self.planted = planted
        self.profile = profile
        self.c_ctype = Counter(cycle_type(c))
        self.c_len_sorted = sorted(cycle_type(c), reverse=True)

    # -- state helpers
    def assign(self, s, v):
        x, xi = self.x, self.xi
        if x[s] != NA: return x[s] == v
        if xi[v] != NA: return False
        if (s < self.m) == (v < self.m): return False        # block structure
        x[s] = v; xi[v] = s; self.trail.append(s)
        return True
    def undo(self, mark):
        x, xi, tr = self.x, self.xi, self.trail
        while len(tr) > mark:
            s = tr.pop(); xi[x[s]] = NA; x[s] = NA

    def propagate(self):
        """Force every x-application that is the single unknown of its chain. Returns ok."""
        x, xi, a, ai, c, N = self.x, self.xi, self.a, self.ai, self.c, self.N
        self.prop_passes += 1
        changed = True
        while changed:
            changed = False
            for j in range(N):
                p1 = x[j]
                if p1 != NA:
                    p2 = x[p1]
                    if p2 != NA:
                        p3 = a[p2]; cur = x[p3]
                        if cur == NA:
                            if not self.assign(p3, c[j]): return False
                            changed = True
                        elif cur != c[j]: return False
                    else:
                        q3 = xi[c[j]]
                        if q3 != NA:
                            if not self.assign(p1, ai[q3]): return False
                            changed = True
                else:
                    q3 = xi[c[j]]
                    if q3 != NA:
                        q1 = xi[ai[q3]]
                        if q1 != NA:
                            if not self.assign(j, q1): return False
                            changed = True
        if self.use_ctype and not self.ctype_ok():
            self.ctype_prunes += 1; return False
        return True

    def ctype_ok(self):
        """Relaxation of ctype(a x^3) = ctype(c): closed cycles of the partial map a x^3 must
        form a sub-multiset of ctype(c); the longest open segment must fit in a remaining cycle."""
        x, a, N = self.x, self.a, self.N
        f = [NA] * N                                   # partial a x^3
        for j in range(N):
            p1 = x[j]
            if p1 == NA: continue
            p2 = x[p1]
            if p2 == NA: continue
            p3 = x[p2]
            if p3 == NA: continue
            f[j] = a[p3]
        indeg = [0] * N
        for j in range(N):
            if f[j] != NA: indeg[f[j]] += 1
        seen = [False] * N
        closed = Counter(); longest_open = 0
        # open segments start at points with indegree 0 (within the partial map)
        for j in range(N):
            if indeg[j] == 0 and f[j] != NA:
                L = 0; k = j
                while k != NA and not seen[k]:
                    seen[k] = True; k = f[k]; L += 1
                # L points on the segment => cycle length >= L (>= L+1 if it does not close)
                longest_open = max(longest_open, L + 1)
        for j in range(N):
            if not seen[j] and f[j] != NA:
                L = 0; k = j
                while not seen[k]:
                    seen[k] = True; k = f[k]; L += 1
                closed[L] += 1
        rem = self.c_ctype.copy()
        for L, k in closed.items():
            if rem[L] < k: return False
            rem[L] -= k
        if longest_open:
            if not any(L >= longest_open and k > 0 for L, k in rem.items()): return False
        return True

    # -- branching heuristics
    def candidates(self, j):
        xi, m = self.xi, self.m
        lo, hi = (m, 2 * m) if j < m else (0, m)
        return [v for v in range(lo, hi) if xi[v] == NA]

    def choose_first(self, free):
        return free[0], self.candidates(free[0])

    def choose_mcf(self, free):
        best = None
        for j in free:
            surv = []
            for v in self.candidates(j):
                mark = len(self.trail)
                ok = self.assign(j, v) and self.propagate()
                self.undo(mark)
                self.lookahead += 1
                if ok: surv.append(v)
            if best is None or len(surv) < len(best[1]):
                best = (j, surv)
                if len(surv) <= 1: break
        return best

    def choose_chains(self, free):
        """Score each free point j by the chains in which the missing image x(j) is one of two
        unknown letters whose source point is already known (the chain becomes forced as soon
        as one more image is fixed); 3-unknown chains give a small weight."""
        x, xi, a, ai, c, N = self.x, self.xi, self.a, self.ai, self.c, self.N
        score = [0.0] * N
        for j in range(N):
            p1 = x[j]
            if p1 != NA:
                if x[p1] == NA:
                    q3 = xi[c[j]]
                    if q3 == NA:                   # unknown x(p1) and x(p3), p3 unknown
                        score[p1] += 1.0
                    # q3 known -> chain was forced by propagation
            else:
                q3 = xi[c[j]]
                if q3 != NA:
                    if xi[ai[q3]] == NA:           # unknown x(j) and x(?) = q2
                        score[j] += 1.0
                else:
                    score[j] += 0.34               # three unknowns
        j = max(free, key=lambda k: (score[k], -k))
        return j, self.candidates(j)

    # -- search
    def solve(self):
        N = self.N
        self.x = [NA] * N; self.xi = [NA] * N; self.trail = []
        self.nodes = 0; self.decisions = 0; self.prop_passes = 0; self.lookahead = 0
        self.ctype_prunes = 0; self.exhausted = False; self.sols = []
        self.nodes_first = None; self.nodes_planted = None
        self.depth_nodes = Counter()               # surviving nodes per decision depth
        self.decision_log = []                     # (depth, domain, surviving, forced, on_planted)
        self.planted_path = []                     # (depth, domain, index of planted value, forced)
        self.wasted = Counter()                    # nodes in off-planted subtrees, by depth
        t0 = time.perf_counter()
        self.nodes += 1
        if self.propagate():
            self._dfs()
        self.seconds = time.perf_counter() - t0
        return self.sols

    def on_planted(self):
        x, p = self.x, self.planted
        return all(x[j] == NA or x[j] == p[j] for j in range(self.N))

    def _dfs(self):
        N = self.N
        choose = {"first": self.choose_first, "mcf": self.choose_mcf,
                  "chains": self.choose_chains}[self.heuristic]
        stack = []
        def open_node(depth):
            free = [j for j in range(N) if self.x[j] == NA]
            if not free:
                xt = tuple(self.x)
                if check_solution(self.a, self.c, xt):
                    self.sols.append(xt)
                    if len(self.sols) == 1: self.nodes_first = self.nodes
                    if xt == self.planted: self.nodes_planted = self.nodes
                return False
            j, cands = choose(free)
            self.decisions += 1
            onp = self.on_planted() if self.profile else False
            stack.append(dict(j=j, cands=cands, idx=0, mark=len(self.trail), depth=depth,
                              forced=0, onp=onp, surv=0, off_start=None))
            return True
        if not open_node(0): return
        while stack:
            fr = stack[-1]
            if self.profile and fr["off_start"] is not None:
                self.wasted[fr["depth"]] += self.nodes - fr["off_start"]; fr["off_start"] = None
            if fr["idx"] >= len(fr["cands"]) or len(self.sols) >= self.max_solutions or self.exhausted:
                stack.pop()
                if self.profile:
                    self.decision_log.append((fr["depth"], len(fr["cands"]), fr["surv"], fr["forced"], fr["onp"]))
                self.undo(fr["mark"])
                continue
            v = fr["cands"][fr["idx"]]; fr["idx"] += 1
            self.undo(fr["mark"])
            before = len(self.trail)
            ok = self.assign(fr["j"], v) and self.propagate()
            if not ok:
                self.undo(before); continue
            self.nodes += 1; fr["surv"] += 1
            nforced = len(self.trail) - before - 1
            fr["forced"] += nforced
            self.depth_nodes[fr["depth"] + 1] += 1
            if self.nodes > self.node_budget:
                self.exhausted = True; continue
            if self.profile and fr["onp"]:
                if self.on_planted():
                    self.planted_path.append((fr["depth"], len(fr["cands"]), fr["idx"] - 1, nforced))
                else:
                    fr["off_start"] = self.nodes - 1
            open_node(fr["depth"] + 1)


def solve_instance(inst, heuristic="first", ctype=False, node_budget=3_000_000,
                   max_solutions=1, profile=False):
    s = WordSolver(inst["a"], inst["c"], inst["m"], heuristic=heuristic, ctype=ctype,
                   node_budget=node_budget, max_solutions=max_solutions,
                   planted=inst["x"] if profile else None, profile=profile)
    sols = s.solve()
    assert all(check_solution(inst["a"], inst["c"], x) for x in sols)
    for x in sols:   # block structure
        assert all((j < inst["m"]) != (x[j] < inst["m"]) for j in range(len(x)))
    return s, sols


# ----------------------------------------------------------------------------- reference solver
def load_reference(zip_dir):
    sys.path.insert(0, zip_dir)
    import round_attack
    return round_attack

def run_reference(inst, round_attack, node_budget):
    t0 = time.perf_counter()
    sols, st = round_attack.solve_round(inst["up"], inst["vp"], inst["A"], inst["B"],
                                        max_solutions=1, node_budget=node_budget)
    dt = time.perf_counter() - t0
    for u, v in sols:
        assert round_attack.round_map(u, v, inst["A"], inst["B"]) == (inst["up"], inst["vp"])
    return dict(nodes=st["nodes"], exhausted=bool(st["exhausted"]), found=len(sols), seconds=dt)

# ----------------------------------------------------------------------------- conjugated form (d)
def cube_roots(z):
    """All x with x^3 = z (polynomial-time enumeration by cycle structure)."""
    N = len(z); seen = [False] * N; cycles = {}
    for i in range(N):
        if not seen[i]:
            cyc = []; j = i
            while not seen[j]: seen[j] = True; cyc.append(j); j = z[j]
            cycles.setdefault(len(cyc), []).append(cyc)
    partials = [[]]
    import itertools
    for L, cs in cycles.items():
        options = []   # list of x-fragments (dict) for this length class
        def build(remaining, acc):
            if not remaining:
                options.append(acc); return
            first = remaining[0]; rest = remaining[1:]
            if L % 3 != 0:                       # single cycle: x = cycle^k, 3k = 1 mod L
                k = pow(3, -1, L)
                frag = {first[i]: first[(i + k) % L] for i in range(L)}
                build(rest, acc + [frag])
            for i2 in range(len(rest)):          # triple: interleave first, B, C
                for i3 in range(i2 + 1, len(rest)):
                    B, C = rest[i2], rest[i3]
                    rest2 = rest[:i2] + rest[i2 + 1:i3] + rest[i3 + 1:]
                    for (P, Q) in ((B, C), (C, B)):
                        for s1 in range(L):
                            for s2 in range(L):
                                frag = {}
                                for i in range(L):
                                    frag[first[i]] = P[(i + s1) % L]
                                    frag[P[(i + s1) % L]] = Q[(i + s2) % L]
                                    frag[Q[(i + s2) % L]] = first[(i + 1) % L]
                                build(rest2, acc + [frag])
        build(cs, [])
        partials = [p + o for p in partials for o in options]
    out = []
    for p in partials:
        x = [NA] * N
        for frag in p:
            for k, v in frag.items(): x[k] = v
        x = tuple(x)
        if mul(x, x, x) == z: out.append(x)
    return out

def conjugated_form_search(inst, budget_y=200_000):
    """Variant (d): x a x^2 = c  <=>  x^3 = a^{-1} y with y = x^{-1} c x in the class of c.
    Branch on the conjugator: y ranges over the S_N-class of c (representatives g^{-1} c g over
    coset representatives of the centraliser), propagate by cube-root extraction of a^{-1} y,
    keep the roots x with x^{-1} c x = y and block structure. Cost = size of the class."""
    a, c, m = inst["a"], inst["c"], inst["m"]; N = len(c); ai = inverse(a)
    # class of c: enumerate y by relabelling; use the canonical cycle form to avoid repeats
    import itertools
    ct = cycle_type(c)
    class_size = math.factorial(N)
    for L, k in Counter(ct).items(): class_size //= (L ** k) * math.factorial(k)
    if class_size > budget_y:
        return dict(class_size=class_size, skipped=True)
    # build all y in the class: y = g c g^{-1} for g in S_N, dedupe via set (only small N)
    ys = set(); found = []; roots_tested = 0
    t0 = time.perf_counter()
    for g in itertools.permutations(range(N)):
        gi = inverse(g); y = mul(g, c, gi)
        if y in ys: continue
        ys.add(y)
        z = compose(ai, y)
        for x in cube_roots(z):
            roots_tested += 1
            if mul(inverse(x), c, x) == y and all((j < m) != (x[j] < m) for j in range(N)):
                found.append(x)
        if len(ys) >= class_size: break
    dt = time.perf_counter() - t0
    assert all(check_solution(a, c, x) for x in found)
    return dict(class_size=class_size, skipped=False, y_enumerated=len(ys),
                roots_tested=roots_tested, solutions=len(found), planted_found=inst["x"] in found,
                seconds=dt)

# ----------------------------------------------------------------------------- Z3 (e)
def z3_solve(inst, timeout_s=300):
    import z3
    a, c, m = inst["a"], inst["c"], inst["m"]; N = len(c)
    opp = lambda j: range(m, N) if j < m else range(0, m)
    P = {(j, v): z3.Bool(f"p_{j}_{v}") for j in range(N) for v in opp(j)}
    s = z3.Solver(); s.set("timeout", int(timeout_s * 1000))
    for j in range(N):
        s.add(z3.PbEq([(P[(j, v)], 1) for v in opp(j)], 1))
    for v in range(N):
        s.add(z3.PbEq([(P[(j, v)], 1) for j in opp(v)], 1))
    for j in range(N):
        for i1 in opp(j):
            for i2 in opp(i1):
                s.add(z3.Implies(z3.And(P[(j, i1)], P[(i1, i2)]), P[(a[i2], c[j])]))
    # parity of U and V (both even): omitted — the propagation solvers do not use it either
    t0 = time.perf_counter(); r = s.check(); dt = time.perf_counter() - t0
    out = dict(status=str(r), seconds=dt, timeout_s=timeout_s)
    if r == z3.sat:
        mdl = s.model(); x = [NA] * N
        for (j, v), b in P.items():
            if z3.is_true(mdl.eval(b, model_completion=True)): x[j] = v
        out["verified"] = check_solution(a, c, tuple(x))
        stats = s.statistics()
        try: out["z3_decisions"] = stats.get_key_value("decisions")
        except Exception: pass
        try: out["z3_conflicts"] = stats.get_key_value("conflicts")
        except Exception: pass
    return out

# ----------------------------------------------------------------------------- fits
def fit_line(xs, ys):
    n = len(xs); xb = sum(xs) / n; yb = sum(ys) / n
    sxx = sum((x - xb) ** 2 for x in xs)
    slope = sum((x - xb) * (y - yb) for x, y in zip(xs, ys)) / sxx if sxx else float("nan")
    return slope, yb - slope * xb

# ----------------------------------------------------------------------------- stages
def stage_variants(args):
    variants = {
        "a_baseline_first": dict(heuristic="first", ctype=False),
        "b1_mcf_lookahead": dict(heuristic="mcf", ctype=False),
        "b2_chains": dict(heuristic="chains", ctype=False),
        "c_first_ctype": dict(heuristic="first", ctype=True),
        "c_chains_ctype": dict(heuristic="chains", ctype=True),
    }
    sel = args.variants.split(",") if args.variants else list(variants)
    rows = []
    for name in sel:
        cfg = variants[name]
        for m in range(args.m_min, args.m_max + 1):
            tm = []
            for t in range(args.trials):
                inst = wreath_instance(m, t)
                s, sols = solve_instance(inst, node_budget=args.budget, **cfg)
                rows.append(dict(variant=name, m=m, t=t, nodes=s.nodes, decisions=s.decisions,
                                 nodes_first=s.nodes_first,
                                 prop_passes=s.prop_passes, lookahead=s.lookahead,
                                 ctype_prunes=s.ctype_prunes, found=len(sols),
                                 exhausted=s.exhausted, seconds=round(s.seconds, 4)))
                tm.append(s.seconds)
                print(f"[{name}] m={m} t={t} nodes={s.nodes} found={len(sols)} "
                      f"exh={s.exhausted} {s.seconds:.1f}s", flush=True)
            if sum(tm) > args.time_cap_per_m: 
                print(f"[{name}] time cap reached at m={m}", flush=True); break
    json.dump(rows, open(args.out, "w"), indent=1)

def stage_reference(args):
    ra = load_reference(args.zip_dir)
    rows = []
    for m in range(args.m_min, args.m_max + 1):
        tm = []
        for t in range(args.trials):
            inst = wreath_instance(m, t)
            r = run_reference(inst, ra, args.budget)
            rows.append(dict(variant="f_reference_solve_round", m=m, t=t, **r))
            tm.append(r["seconds"])
            print(f"[ref] m={m} t={t} nodes={r['nodes']} exh={r['exhausted']} {r['seconds']:.1f}s", flush=True)
        if sum(tm) > args.time_cap_per_m: break
    json.dump(rows, open(args.out, "w"), indent=1)

def stage_z3(args):
    rows = []
    for m in range(args.m_min, args.m_max + 1):
        n_to = 0
        for t in range(args.trials):
            inst = wreath_instance(m, t)
            r = z3_solve(inst, timeout_s=args.timeout)
            rows.append(dict(variant="e_z3_onehot", m=m, t=t, **r))
            print(f"[z3] m={m} t={t} {r}", flush=True)
            json.dump(rows, open(args.out, "w"), indent=1)
            if r["status"] != "sat": n_to += 1
        if n_to == args.trials: break
    json.dump(rows, open(args.out, "w"), indent=1)

def stage_conjugated(args):
    rows = []
    for m in range(args.m_min, args.m_max + 1):
        for t in range(args.trials):
            inst = wreath_instance(m, t)
            r = conjugated_form_search(inst, budget_y=args.budget)
            rows.append(dict(variant="d_conjugated_cuberoot", m=m, t=t, **r))
            print(f"[conj] m={m} t={t} {r}", flush=True)
            if r.get("skipped"): break
    json.dump(rows, open(args.out, "w"), indent=1)

def stage_profile(args):
    out = {}
    for heur in args.variants.split(","):
        per_m = {}
        for m in range(args.m_min, args.m_max + 1):
            recs = []
            for t in range(args.trials):
                inst = wreath_instance(m, t)
                s, sols = solve_instance(inst, heuristic=heur, node_budget=args.budget,
                                         max_solutions=10**9, profile=True)
                internal = [d for d in s.decision_log]        # (depth, dom, surv, forced, onp)
                surv = [d[2] for d in internal]
                dom = [d[1] for d in internal]
                forced = [d[3] for d in internal]
                pp = s.planted_path
                recs.append(dict(
                    t=t, nodes=s.nodes, nodes_first=s.nodes_first, nodes_planted=s.nodes_planted,
                    decisions=s.decisions, found=len(sols), exhausted=s.exhausted,
                    seconds=round(s.seconds, 4),
                    mean_surviving_children=statistics.mean(surv) if surv else None,
                    mean_domain=statistics.mean(dom) if dom else None,
                    mean_forced_per_node=statistics.mean(forced) / max(1, statistics.mean(surv)) if surv else None,
                    planted_depth=len(pp),
                    planted_domains=[p[1] for p in pp],
                    planted_index=[p[2] for p in pp],
                    planted_forced=[p[3] for p in pp],
                    depth_nodes=dict(sorted(s.depth_nodes.items())),
                    wasted_by_depth=dict(sorted(s.wasted.items())),
                    surv_hist=dict(sorted(Counter(surv).items())),
                ))
                print(f"[profile {heur}] m={m} t={t} nodes={s.nodes} planted_depth={len(pp)} "
                      f"surv={recs[-1]['mean_surviving_children']:.2f} {s.seconds:.1f}s", flush=True)
            per_m[str(m)] = recs
        out[heur] = per_m
    json.dump(out, open(args.out, "w"), indent=1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["variants", "reference", "z3", "conjugated", "profile"])
    ap.add_argument("--m-min", type=int, default=6); ap.add_argument("--m-max", type=int, default=13)
    ap.add_argument("--trials", type=int, default=5); ap.add_argument("--budget", type=int, default=3_000_000)
    ap.add_argument("--variants", default=""); ap.add_argument("--timeout", type=float, default=300)
    ap.add_argument("--time-cap-per-m", type=float, default=600)
    ap.add_argument("--zip-dir", default="zip/cpwi-post-quantum-dynamical-inversion")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    dict(variants=stage_variants, reference=stage_reference, z3=stage_z3,
         conjugated=stage_conjugated, profile=stage_profile)[args.stage](args)

if __name__ == "__main__":
    main()
