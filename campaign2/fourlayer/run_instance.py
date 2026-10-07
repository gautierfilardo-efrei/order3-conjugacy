import sys, json, time
sys.path.insert(0, '.')
import fourlayer_track as FT
L = int(sys.argv[1]); m = tuple(int(t) for t in sys.argv[2].split(',')); e = tuple(int(t) for t in sys.argv[3].split(','))
cap = float(sys.argv[4]); out = sys.argv[5]
a, c = FT.affine_instance(L, m); ap, cp = FT.pad(a, c, e)
k = len(m); M = len(e)
wit = FT.conj3_affine_witnesses(L, m)
t0 = time.time()
sols, st = FT.P3Solver(ap, cp).solve(time_cap=cap)
rec = {"L": L, "m": list(m), "e": list(e), "N": len(ap), "ctype_a": FT.ctype(ap), "ctype_c": FT.ctype(cp),
       "conj3_witnesses": len(wit), "conj3_yes": len(wit) > 0,
       "n_solutions": len(sols), "complete": st["complete"], "nodes": st["nodes"], "time_s": round(st["time"],1),
       "all_verified": all(FT.is_p3_solution(ap, cp, x) for x in sols),
       "ctypes": sorted({str(FT.ctype(x)) for x in sols}),
       "planted_lifts_found": sum(1 for y in wit if FT.lift(y, M) in set(sols)),
       "layer_affine_count": sum(1 for x in sols if (FT.is_layer_affine(x, L, k, M) or {}).get("affine", False)),
       "solutions": [list(x) for x in sols]}
json.dump(rec, open(out, "w"))
print(out, rec["n_solutions"], rec["complete"], rec["nodes"], rec["time_s"])
