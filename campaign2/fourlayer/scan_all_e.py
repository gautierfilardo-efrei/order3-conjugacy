import sys, json, itertools
sys.path.insert(0, '.')
import fourlayer_track as FT
L, M = 5, 4
v = (1,2,4,3)
adm = [e for e in itertools.product(range(1,L), repeat=M) if sum(e) % L and sum(v[r]*e[r] for r in range(M)) % L]
assert len(adm) == 160
# orbits under rotation (relabel layers) and scaling by units (conjugation by z -> u z blockwise)
seen = set(); reps = []
for e in adm:
    if e in seen: continue
    orb = set()
    for u in range(1, L):
        f = tuple(u*t % L for t in e)
        for s in range(M): orb.add(f[s:]+f[:s])
    assert orb <= set(adm)
    seen |= orb; reps.append((e, len(orb)))
print("orbit representatives:", len(reps), reps, flush=True)
out = {"admissible": 160, "orbits": []}
for e, osz in reps:
    row = {"e": e, "orbit_size": osz}
    for m in (2,3,4):
        a, c = FT.affine_instance(L, (m,)); ap, cp = FT.pad(a, c, e)
        sols, st = FT.P3Solver(ap, cp).solve(time_cap=300)
        row[f"m{m}"] = {"n": len(sols), "complete": st["complete"], "ctype_a": FT.ctype(ap), "ctype_c": FT.ctype(cp),
                        "ctypes": sorted({str(FT.ctype(x)) for x in sols}), "example": list(sols[0]) if sols else None}
    print(row["e"], osz, [(m, row[f"m{m}"]["n"], row[f"m{m}"]["complete"]) for m in (2,3,4)], flush=True)
    out["orbits"].append(row); json.dump(out, open("scan_all_e_k1.json","w"), indent=1)
