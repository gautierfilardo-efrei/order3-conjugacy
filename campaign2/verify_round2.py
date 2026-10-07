"""Independent re-verification (round 2) of the publication-campaign results, written from the statements.
Conventions: image tuples, (pq)(j)=p(q(j))."""
import itertools, random, json, math, time
from collections import Counter, defaultdict
def comp(p,q): return tuple(p[q[j]] for j in range(len(p)))
def inv(p):
    r=[0]*len(p)
    for i,v in enumerate(p): r[v]=i
    return tuple(r)
def power(p,k):
    r=tuple(range(len(p)))
    for _ in range(k): r=comp(r,p)
    return r
def cycles(p):
    seen=[False]*len(p); out=[]
    for i in range(len(p)):
        if not seen[i]:
            c=[]; j=i
            while not seen[j]: seen[j]=True; c.append(j); j=p[j]
            out.append(c)
    return out
def ctype(p): return tuple(sorted((len(c) for c in cycles(p)),reverse=True))
def order(p):
    o=1
    for c in cycles(p): o=o*len(c)//math.gcd(o,len(c))
    return o
rng=random.Random(7); out={}

# ---------- 1. lem:cuberoots count formula vs brute force (N<=7)
def cube_count_formula(t):
    cnt=Counter(len(c) for c in cycles(t)); tot=1
    for l,n in cnt.items():
        s=0
        for k in range(n//3+1):
            sig = 1 if l%3 else 0
            term=math.factorial(n)//(math.factorial(k)*6**k*math.factorial(n-3*k))*(2*l*l)**k*(sig**(n-3*k) if not (sig==0 and n-3*k==0) else 1)
            s+=term
        tot*=s
    return tot
ok=tot=0
for N in range(1,8):
    reps={}
    for t in itertools.permutations(range(N)): reps.setdefault(ctype(t),t)
    for t in reps.values():
        bf=sum(1 for x in itertools.permutations(range(N)) if power(x,3)==t)
        ok+=(bf==cube_count_formula(t)); tot+=1
out['cuberoots']=dict(classes=tot,agree=ok)

# ---------- 2. lem:cent criterion (i) and count (ii) vs brute force, all a (one per class), all c in C(a), N<=7
def centraliser(a):
    N=len(a); return [x for x in itertools.permutations(range(N)) if comp(x,a)==comp(a,x)]
def cent_formula(a,t):
    # classes (L, ell, R)
    acyc=cycles(a); cyc_of={}
    for idx,c in enumerate(acyc):
        for j in c: cyc_of[j]=idx
    # rho: permutation of a-cycles induced by t (t in C(a) maps cycles to cycles)
    rho={idx: cyc_of[t[c[0]]] for idx,c in enumerate(acyc)}
    seen=set(); classes=Counter()
    for idx in range(len(acyc)):
        if idx in seen: continue
        orb=[]; j=idx
        while j not in seen: seen.add(j); orb.append(j); j=rho[j]
        ell=len(orb); L=len(acyc[idx]); p0=acyc[idx][0]
        q=p0
        for _ in range(ell): q=t[q]
        # q = a^R p0
        R=0; z=p0
        while z!=q: z=a[z]; R+=1
        classes[(L,ell,R%L)]+=1
    total=1
    for (L,ell,R),n in classes.items():
        free = (ell%3!=0) and (L%3!=0 or R%3==0)
        sig = 0 if not free else (1 if L%3 else 3)
        s=0
        for k in range(n//3+1):
            term=math.factorial(n)//(math.factorial(k)*6**k*math.factorial(n-3*k))*(2*ell*ell*L*L)**k
            rest=n-3*k
            term*= (sig**rest if rest>0 else 1)
            s+=term
        total*=s
    return total
ok=tot=pos=0
for N in range(2,8):
    reps={}
    for a in itertools.permutations(range(N)): reps.setdefault(ctype(a),a)
    for a in reps.values():
        C=centraliser(a); Cset=set(C)
        for t in C:
            bf=sum(1 for x in C if power(x,3)==t)
            f=cent_formula(a,t); ok+=(bf==f); tot+=1; pos+=(bf>0)
out['cent']=dict(pairs=tot,agree=ok,positive=pos)
# plus: commuting solutions of x a x^2 = c equal cube roots of a^{-1}c in C(a), brute force N<=6 random c
ok=tot=0
for _ in range(200):
    N=6; a=tuple(rng.sample(range(N),N)); c=tuple(rng.sample(range(N),N))
    bf=sum(1 for x in itertools.permutations(range(N)) if comp(x,comp(a,comp(x,x)))==c and comp(x,a)==comp(a,x))
    t=comp(inv(a),c); f = cent_formula(a,t) if comp(t,a)==comp(a,t) else 0
    ok+=(bf==f); tot+=1
out['cent_equation']=dict(checks=tot,agree=ok)

# ---------- 3. Lemma A (affine criterion for x^n = 1, any n) vs exhaustive conjugators, composite n
def build(L,ms):
    k=len(ms); a=tuple(i*L+(z+1)%L for i in range(k) for z in range(L)); c=tuple(i*L+(z+ms[i])%L for i in range(k) for z in range(L)); return a,c
def conjugators(a,c,L,k):
    N=L*k
    for s in itertools.permutations(range(k)):
        for base in itertools.product(range(L),repeat=k):
            x=[0]*N
            for i in range(k):
                p=i*L; q=s[i]*L+base[i]
                for _ in range(L): x[p]=q; p=a[p]; q=c[q]
            yield tuple(x)
def critA(L,ms,n):
    k=len(ms)
    for s in itertools.permutations(range(k)):
        if order(s)==0 or n%order(s): continue
        ok=True
        for cyc in cycles(s):
            pr=1
            for i in cyc: pr=pr*ms[i]%L
            if pow(pr,n//len(cyc),L)!=1: ok=False; break
        if ok: return True
    return False
res=[]
for L,k in ((5,3),(7,3),(4,3),(9,2),(5,4),(8,2),(6,3)):
    units=[m for m in range(1,L) if math.gcd(m,L)==1]
    for _ in range(8):
        ms=[rng.choice(units) for _ in range(k)]; a,c=build(L,ms); idn=tuple(range(L*k))
        orders=set(order(x) for x in conjugators(a,c,L,k))
        for n in (2,3,4,6,8,9,12):
            ex=any(n%o==0 for o in orders); cr=critA(L,ms,n); res.append((L,k,n,ex,cr))
out['lemmaA']=dict(instances=len(res),agree=sum(1 for r in res if r[3]==r[4]),positives=sum(1 for r in res if r[3]))

# ---------- 4. Proposition E witnesses: conjugator orders of (12345)->a^2 and 7-cycle -> a^3
for N,u in ((5,2),(7,3)):
    a=tuple((j+1)%N for j in range(N)); c=power(a,u)
    ords=sorted(set(order(x) for x in itertools.permutations(range(N)) if comp(x,comp(a,inv(x)))==c))
    out[f'propE_N{N}_u{u}']=ords

# ---------- 5. Four-layer witness (k=1, L=5, e=(1,2,3,3), c=a^{-1}): verify the explicit x, and count all solutions by brute force structure? N=20 too big for brute force; verify the witness and the Conj_3-NO claim.
ap=(6,7,8,9,5,12,13,14,10,11,18,19,15,16,17,3,4,0,1,2); cp=(9,5,6,7,8,13,14,10,11,12,17,18,19,15,16,2,3,4,0,1)
x=(0,13,1,7,10,14,9,15,5,4,8,3,17,2,18,16,6,11,12,19)
# rebuild a', c' from the definition to confirm the encoding (j,r)->5r+j
a5=tuple((j+1)%5 for j in range(5)); c5=inv(a5); e=(1,2,3,3)
def padded(base):
    P=[0]*20
    for r in range(4):
        for j in range(5):
            P[5*r+j]=5*((r+1)%4)+power(base,e[r])[j]
    return tuple(P)
out['fourlayer']=dict(a_prime_matches=padded(a5)==ap, c_prime_matches=padded(c5)==cp, witness_solves=comp(x,comp(ap,comp(x,x)))==cp, witness_ctype=ctype(x),
   conj3_no=not any(power(y,3)==tuple(range(5)) and comp(y,comp(a5,inv(y)))==c5 for y in itertools.permutations(range(5))),
   unpadded_P3_count=sum(1 for y in itertools.permutations(range(5)) if comp(y,comp(a5,comp(y,y)))==c5))

# ---------- 6. Audit: rigid-pair census (solvable pairs all of whose solutions have order exactly 3), weighted by class size, N=5,6
rig={}
for N in (5,6):
    reps=Counter()
    for a in itertools.permutations(range(N)): reps[ctype(a)]+=1
    rep_of={}
    for a in itertools.permutations(range(N)): rep_of.setdefault(ctype(a),a)
    solv=rigid=0
    for ct,a in rep_of.items():
        w=reps[ct]
        for c in itertools.permutations(range(N)):
            sols=[y for y in itertools.permutations(range(N)) if comp(y,comp(a,comp(y,y)))==c]
            if sols:
                solv+=w; rigid+= w*all(order(y)==3 for y in sols)
    rig[N]=dict(solvable_pairs=solv, rigid=rigid, fraction=rigid/solv)
out['rigid_census_N5_6']=rig
json.dump(out,open('verify_round2.json','w'),indent=1,default=str); print(json.dumps(out,default=str))
