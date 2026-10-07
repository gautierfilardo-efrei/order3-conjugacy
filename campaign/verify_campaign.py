"""Independent re-verification (main session) of the load-bearing claims of the P3 campaign sub-agents.
Everything here is written from the THEOREM STATEMENTS, not from the sub-agents' code.
Conventions: permutations are image tuples, (pq)(j)=p(q(j))."""
import itertools, random, json, math, time
from collections import Counter
def comp(p,q): return tuple(p[q[j]] for j in range(len(p)))
def inv(p):
    r=[0]*len(p)
    for i,v in enumerate(p): r[v]=i
    return tuple(r)
def power(p,k):
    r=tuple(range(len(p)))
    for _ in range(k): r=comp(r,p)
    return r
def ctype(p):
    seen=[False]*len(p); t=[]
    for i in range(len(p)):
        if not seen[i]:
            l=0; j=i
            while not seen[j]: seen[j]=True; j=p[j]; l+=1
            t.append(l)
    return tuple(sorted(t,reverse=True))
def nsol(a,c):
    N=len(a); return sum(1 for x in itertools.permutations(range(N)) if comp(x,comp(a,comp(x,x)))==c)
out={}

# ---------- 1. Conj_2 (Theorem 2): own implementation of the criterion vs brute force over all involutions
def orbits(a,c):
    N=len(a); seen=[False]*N; orbs=[]
    for s in range(N):
        if seen[s]: continue
        st=[s]; seen[s]=True; o=[]
        while st:
            u=st.pop(); o.append(u)
            for v in (a[u],c[u]):
                if not seen[v]: seen[v]=True; st.append(v)
        orbs.append(sorted(o))
    return orbs
def iso(a,c,O,O2,j0,v,twisted):
    """unique bijection O->O2 with j0->v commuting with (a,c) [plain] or intertwining a<->c [twisted]; None if inconsistent"""
    A,C=(c,a) if twisted else (a,c)
    phi={j0:v}; st=[j0]
    while st:
        u=st.pop()
        for g,h in ((a,A),(c,C)):
            u2=g[u]; v2=h[phi[u]]
            if u2 in phi:
                if phi[u2]!=v2: return None
            else: phi[u2]=v2; st.append(u2)
    if len(phi)!=len(O) or len(set(phi.values()))!=len(O) or set(phi.values())!=set(O2): return None
    return phi
def exists_iso(a,c,O,O2,twisted):
    if len(O)!=len(O2): return None
    for v in O2:
        phi=iso(a,c,O,O2,O[0],v,twisted)
        if phi is not None: return phi
    return None
def conj2_criterion(a,c):
    orbs=orbits(a,c)
    # plain-isomorphism classes
    classes=[]
    for O in orbs:
        for cl in classes:
            if exists_iso(a,c,cl[0],O,False) is not None: cl.append(O); break
        else: classes.append([O])
    for cl in classes:
        # twisted partner class
        partner=None
        for cl2 in classes:
            if exists_iso(a,c,cl[0],cl2[0],True) is not None: partner=cl2; break
        if partner is None or len(partner)!=len(cl): return False
        if partner is cl and len(cl)%2==1:
            O=cl[0]; ok=False
            for v in O:
                phi=iso(a,c,O,O,O[0],v,True)
                if phi is not None and phi[phi[O[0]]]==O[0]: ok=True; break
            if not ok: return False
    return True
def conj2_brute(a,c):
    N=len(a)
    for x in itertools.permutations(range(N)):
        if comp(x,x)==tuple(range(N)) and comp(x,comp(a,inv(x)))==c: return True
    return False
rng=random.Random(20261007); t0=time.time(); agree=0; tot=0; yes=0
for N in range(3,8):
    for _ in range(120 if N<7 else 80):
        a=tuple(rng.sample(range(N),N)); g=tuple(rng.sample(range(N),N)); c=comp(g,comp(a,inv(g)))
        b=conj2_brute(a,c); k=conj2_criterion(a,c); agree+=(b==k); tot+=1; yes+=b
# plus ALL pairs at N=5 (a one representative per class, all conjugates c)
for N in (4,5):
    reps={}
    for a in itertools.permutations(range(N)): reps.setdefault(ctype(a),a)
    for a in reps.values():
        for g in itertools.permutations(range(N)):
            c=comp(g,comp(a,inv(g))); b=conj2_brute(a,c); k=conj2_criterion(a,c); agree+=(b==k); tot+=1; yes+=b
out['conj2']=dict(pairs=tot,agree=agree,yes=yes,seconds=round(time.time()-t0,1))

# ---------- 2. Affine criterion for x^p=1 (Lemma 2) at p=5 and p=3, vs exhaustive conjugator enumeration
def build(L,ms):
    k=len(ms); a=tuple(i*L+(z+1)%L for i in range(k) for z in range(L)); c=tuple(i*L+(z+ms[i])%L for i in range(k) for z in range(L)); return a,c
def all_conjugators(a,c,L,k):
    N=L*k
    for s in itertools.permutations(range(k)):
        for base in itertools.product(range(L),repeat=k):
            x=[0]*N
            for i in range(k):
                p=i*L; q=s[i]*L+base[i]
                for _ in range(L): x[p]=q; p=a[p]; q=c[q]
            yield tuple(x)
def crit(L,ms,p):
    k=len(ms)
    for s in itertools.permutations(range(k)):
        if any(power_perm(s,p)[i]!=i for i in range(k)): continue
        ok=True; seen=set()
        for i in range(k):
            if i in seen: continue
            cyc=[]; j=i
            while j not in seen: seen.add(j); cyc.append(j); j=s[j]
            if len(cyc)==1:
                if pow(ms[i],p,L)!=1: ok=False; break
            else:
                pr=1
                for j in cyc: pr=pr*ms[j]%L
                if pr!=1: ok=False; break
        if ok: return True
    return False
def power_perm(s,p):
    r=tuple(range(len(s)))
    for _ in range(p): r=tuple(s[r[i]] for i in range(len(s)))
    return r
res2=[]; t0=time.time()
for p,L,k,trials in ((5,7,5,12),(5,4,5,12),(3,7,3,20),(3,5,4,20),(5,7,4,10)):
    units=[m for m in range(1,L) if math.gcd(m,L)==1]
    for t in range(trials):
        ms=[rng.choice(units) for _ in range(k)]
        if t<trials//2: ms=[rng.choice([u for u in units if u!=1]) for _ in range(k)]
        a,c=build(L,ms); N=L*k; idn=tuple(range(N))
        sols=[x for x in all_conjugators(a,c,L,k) if power(x,p)==idn]
        ex=bool(sols); cr=crit(L,ms,p)
        hyp = all(m!=1 for m in ms) and not any(pow(m,p,L)==1 and m%L!=1 for m in range(1,L) if math.gcd(m,L)==1)
        types_ok = all(ctype(x)==tuple([p]*(N//p)) for x in sols) if hyp else None
        res2.append(dict(p=p,L=L,k=k,ms=ms,exhaustive=ex,criterion=cr,agree=ex==cr,all_mi_ne1=all(m!=1 for m in ms),fixed_point_free_type=types_ok))
out['affine_criterion']=dict(instances=len(res2),agree=sum(r['agree'] for r in res2),positives=sum(r['exhaustive'] for r in res2),
   type_claim_checked=sum(1 for r in res2 if r['fixed_point_free_type'] is not None and r['exhaustive']),type_claim_hyp_violated_instances=sum(1 for r in res2 if r['fixed_point_free_type'] is None and r['all_mi_ne1']),type_claim_ok=all(r['fixed_point_free_type'] in (None,True) for r in res2),seconds=round(time.time()-t0,1))

# ---------- 3. Theorem C.1 (N=7 PSL(2,7) obstruction)
a=(1,2,3,4,5,6,0); c1=(0,4,5,1,3,6,2); c2=(0,2,5,6,3,1,4)
n1,n2=nsol(a,c1),nsol(a,c2)
def group_with_words(a,c):
    idn=tuple(range(7)); G={idn:()}; frontier=[idn]
    while frontier:
        nf=[]
        for g in frontier:
            for sym,h in (('a',a),('c',c)):
                g2=comp(h,g)
                if g2 not in G: G[g2]=(sym,)+G[g]; nf.append(g2)
        frontier=nf
    return G
G1=group_with_words(a,c1); G2=group_with_words(a,c2)
def evalword(w,a,c):
    g=tuple(range(7))
    for sym in reversed(w): g=comp(a if sym=='a' else c,g)
    return g
# phi: g=w(a,c1) -> w(a,c2); well-defined iff w(a,c1)=w'(a,c1) implies w(a,c2)=w'(a,c2): check via all words stored + check it's a bijection and homomorphism
phi={g:evalword(w,a,c2) for g,w in G1.items()}
well_defined=True
# test well-definedness on products: phi(gh) == phi(g)phi(h) for all g,h (homomorphism implies well-defined on the word representation)
hom=all(phi[comp(g,h)]==comp(phi[g],phi[h]) for g in G1 for h in G1)
bij=len(set(phi.values()))==len(G1) and set(phi.values())==set(G2)
ctypes_ok=all(ctype(g)==ctype(phi[g]) for g in G1)
out['theorem_C1']=dict(order_G1=len(G1),order_G2=len(G2),nsol_c1=n1,nsol_c2=n2,ctype_c1=ctype(c1),ctype_c2=ctype(c2),phi_is_isomorphism=hom and bij,phi_preserves_ctypes=ctypes_ok,
   transitive=(len(orbits(a,c1))==1 and len(orbits(a,c2))==1))

# ---------- 4. Theorem N6 (transpositions, a = N-cycle)
def eps(m):
    if m==1: return 0
    return (-pow(3,-1,m))%m if math.gcd(3,m)==1 else None
def N6(N,d):
    cnt=0
    for l in range(1,N):
        if l%3==0 or (N-l)%3==0: continue
        if (l+eps(N-l)-eps(l))%N==d%N: cnt+=1
    return cnt
n6=[]
for N in range(4,10):
    a=tuple((j+1)%N for j in range(N))
    for d in range(1,N):
        c=list(range(N)); c[0],c[d]=d,0; c=tuple(c)
        n6.append((N,d,nsol(a,c),N6(N,d)))
out['theorem_N6']=dict(cases=len(n6),agree=sum(1 for r in n6 if r[2]==r[3]),solvable=sum(1 for r in n6 if r[2]>0))

# ---------- 5. Lemma P2: {x a x^2 = c, x a^2 x^2 = c^2} <=> {x^3=1, x a x^-1 = c}  (random check N=7, plus algebraic proof in notes)
cnt=0; tot=0
for _ in range(300):
    N=7; a=tuple(rng.sample(range(N),N)); x=tuple(rng.sample(range(N),N)); c=comp(x,comp(a,comp(x,x)))
    lhs = comp(x,comp(power(a,2),comp(x,x)))==comp(c,c)
    rhs = power(x,3)==tuple(range(N)) and comp(x,comp(a,inv(x)))==c
    cnt+=(lhs==rhs); tot+=1
# also planted order-3 x
for _ in range(100):
    N=6; 
    while True:
        x=tuple(rng.sample(range(N),N))
        if power(x,3)==tuple(range(N)): break
    a=tuple(rng.sample(range(N),N)); c=comp(x,comp(a,inv(x)))
    lhs = comp(x,comp(a,comp(x,x)))==c and comp(x,comp(power(a,2),comp(x,x)))==comp(c,c)
    cnt+=lhs; tot+=1
out['lemma_P2']=dict(checks=tot,agree=cnt)
json.dump(out,open('verify_campaign.json','w'),indent=1,default=str)
print(json.dumps(out,default=str))
