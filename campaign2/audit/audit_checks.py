"""audit_checks.py -- independent numerical checks performed during the referee audit of p3_article.tex v6.
Standard library + numpy.  Each block prints the number of instances checked and the number of discrepancies.
Checks: (1) Theorem transposition, N=2..8, all d;  (2) Lemma cent criterion, (L,n) in {(2,4),(4,3),(5,3),(7,3),(4,4),(5,4),(2,6)};
(3) Theorem conj2 criterion vs brute force, exhaustive N<=6 (one a per class, every c~a) + 300 random N=7;
(4) Theorem psl certificate;  (5) rho_3 formula N<=7 and closed-form first moment N<=6;  (6) Lemma affinep for composite and
prime L (incl. p=2) vs enumeration of all conjugators;  (7) ten-case analysis of Theorem main, B<=12;  (8) admissible p-sets of
Theorem conjp for p=3,5,7 (arithmetic);  (9) solvable fractions / solution-count law N=4..9;  (10) complete bounded-support
algorithm (replacement proof of Lemma support) vs brute force;  (11) rigid-pair census N=5..8.
"""
import itertools, math, random, json, time
from collections import Counter, defaultdict
import numpy as np
random.seed(1)
PERMS={n:list(itertools.permutations(range(n))) for n in range(1,9)}
def comp(p,q): return tuple(p[q[j]] for j in range(len(p)))


def inv(p):
    r=[0]*len(p)
    for i,v in enumerate(p): r[v]=i
    return tuple(r)


def ctype(p):
    n=len(p); seen=[False]*n; cyc=[]
    for i in range(n):
        if not seen[i]:
            l=0;j=i
            while not seen[j]: seen[j]=True; j=p[j]; l+=1
            cyc.append(l)
    return tuple(sorted(cyc,reverse=True))


def cycles(p):
    n=len(p); seen=[False]*n; out=[]
    for i in range(n):
        if not seen[i]:
            c=[];j=i
            while not seen[j]: seen[j]=True; c.append(j); j=p[j]
            out.append(c)
    return out


def P3count(a,c,perms):
    return sum(1 for x in perms if comp(x,comp(a,comp(x,x)))==c)


def eps(m):
    if m==1: return 0
    return (-pow(3,-1,m))%m


def formula(N,d):
    return sum(1 for l in range(1,N) if l%3 and (N-l)%3 and (l+eps(N-l)-eps(l)-d)%N==0)


def gen_group(gens):
    e=tuple(range(len(gens[0]))); seen={e:()}; frontier=[e]
    while frontier:
        nf=[]
        for g in frontier:
            for i,s in enumerate(gens):
                h=comp(s,g)
                if h not in seen: seen[h]=seen[g]+(i,); nf.append(h)
        frontier=nf
    return seen


def evalword(w,gens):
    g=tuple(range(7))
    for i in w: g=comp(gens[i],g)
    return g


def rho3_formula(nu):
    m=Counter(nu); r=1
    for L,mL in m.items():
        s=0
        for k in range(0,mL//3+1):
            if L%3==0 and 3*k!=mL: continue
            s+=math.factorial(mL)//(math.factorial(mL-3*k)*math.factorial(k)*6**k)*(2*L*L)**k
        r*=s
    return r


def wreath_elements(L,n):
    # elements of Z_L wr Sym(n) as permutations of n*L points: block i -> block rho(i), z -> z + r_i
    out=[]
    for rho in itertools.permutations(range(n)):
        for r in itertools.product(range(L),repeat=n):
            p=[0]*(n*L)
            for i in range(n):
                for z in range(L):
                    p[i*L+z]=rho[i]*L+(z+r[i])%L
            out.append((tuple(p),rho,r))
    return out


def cent_criterion(L,n,rho,r):
    # rho-cycles of length l divisible by 3, rotation sum class counts divisible by 3
    cnt=Counter()
    for cyc in cycles(rho):
        if len(cyc)%3==0:
            cnt[(len(cyc),sum(r[i] for i in cyc)%L)]+=1
    return all(v%3==0 for v in cnt.values())


def orbits_of(gens,N):
    seen=[-1]*N; orbs=[]
    for i in range(N):
        if seen[i]<0:
            o=[i]; seen[i]=len(orbs); k=0
            while k<len(o):
                u=o[k]; k+=1
                for g in gens:
                    v=g[u]
                    if seen[v]<0: seen[v]=len(orbs); o.append(v)
            orbs.append(o)
    return orbs


def iso_test(O,Op,j0,v,a,c,twisted):
    # try to close x(j0)=v with x(a u)=a' x(u), x(c u)=c' x(u) where (a',c')=(c,a) if twisted else (a,c)
    ap,cp=(c,a) if twisted else (a,c)
    if len(O)!=len(Op): return None
    x={j0:v}; stack=[j0]
    while stack:
        u=stack.pop()
        for g,gp in ((a,ap),(c,cp)):
            w=g[u]; img=gp[x[u]]
            if w in x:
                if x[w]!=img: return None
            else:
                x[w]=img; stack.append(w)
    if len(set(x.values()))!=len(O) or set(x.values())!=set(Op): return None
    return x


def conj2_criterion(a,c):
    N=len(a); orbs=orbits_of([a,c],N); k=len(orbs)
    def plain_iso(i,j): return any(iso_test(orbs[i],orbs[j],orbs[i][0],v,a,c,False) for v in orbs[j])
    def tw_iso(i,j): return any(iso_test(orbs[i],orbs[j],orbs[i][0],v,a,c,True) for v in orbs[j])
    cls=[-1]*k; classes=[]
    for i in range(k):
        if cls[i]<0:
            cls[i]=len(classes); members=[i]
            for j in range(i+1,k):
                if cls[j]<0 and plain_iso(i,j): cls[j]=cls[i]; members.append(j)
            classes.append(members)
    for C in classes:
        i=C[0]
        partner=[cl for cl,D in enumerate(classes) if tw_iso(i,D[0])]
        if not partner: return False
        D=classes[partner[0]]
        if len(D)!=len(C): return False
        if D==C and len(C)%2==1:
            O=orbs[i]; ok=False
            for v in O:
                x=iso_test(O,O,O[0],v,a,c,True)
                if x and x[x[O[0]]]==O[0]: ok=True;break
            if not ok: return False
    return True


def blockwise(L,k,m):
    a=[0]*(k*L); c=[0]*(k*L)
    for i in range(k):
        for z in range(L):
            a[i*L+z]=i*L+(z+1)%L; c[i*L+z]=i*L+(z+m[i])%L
    return tuple(a),tuple(c)


def conjugators(L,k,m):
    for sigma in itertools.permutations(range(k)):
        for t in itertools.product(range(L),repeat=k):
            x=[0]*(k*L)
            for i in range(k):
                for z in range(L):
                    x[i*L+z]=sigma[i]*L+(m[sigma[i]]*z+t[i])%L
            yield tuple(x),sigma,t


def crit(L,k,m,p):
    for sigma in itertools.permutations(range(k)):
        ok=True
        for cyc in cycles(sigma):
            if len(cyc)==1:
                if pow(m[cyc[0]],p,L)!=1: ok=False
            elif len(cyc)==p:
                if math.prod(m[i] for i in cyc)%L!=1: ok=False
            else: ok=False
        if ok: return True
    return False


def perm_power(x,p):
    e=tuple(range(len(x)))
    for _ in range(p): e=comp(x,e)
    return e


def main_case(B,T):
    M=3*B
    for al,be,ga in [(x,y,3-x-y) for x in range(4) for y in range(4-x)]:
        for S in range(3,3*B+1):
            z=(S+be*M+ga*(T-M-B))%T==0
            # admissible iff (1,1,1) and S==B ; but S must be realisable by al,be,ga sizes in [1,B]: always realisable if 3<=S<=3B
            if z!=((al,be,ga)==(1,1,1) and S==B): return (al,be,ga,S)
    return None


def conjp_arith(p,Bp):
    M=p*Bp+1; w=[(p+1)**(j) for j in range(p-1)]; W=sum(w); T0=M*(p+1)**(p-1)
    T=1
    while T<T0: T*=2
    # enumerate multisets of (class j, size s) of size p
    items=[(j,s) for j in range(p) for s in range(1,Bp+1)]
    bad=0; tot=0
    for S in itertools.combinations_with_replacement(items,p):
        e=sum((s+M*w[j]) if j<p-1 else (s-M*W-Bp) for j,s in S)
        beta=Counter(j for j,s in S); Sigma=sum(s for j,s in S)
        adm=(e%T==0); match=(all(beta[j]==1 for j in range(p)) and Sigma==Bp); tot+=1; bad+=adm!=match
    return tot,bad


def cube_root_with_constraints(t,f):
    N=len(t); cyc=cycles(t); cid={}
    for i,C in enumerate(cyc):
        for z in C: cid[z]=i
    Ln=[len(C) for C in cyc]
    tpow=lambda z,i: [z:=t[z] for _ in range(i)][-1] if i>0 else z
    xmap={}; pi={}
    def set_cycle(Ci,z0,w0):
        # define x on cycle Ci by x(t^i z0)=t^i w0 ; return False if inconsistent
        z,w=z0,w0
        for _ in range(Ln[Ci]):
            if z in xmap and xmap[z]!=w: return False
            xmap[z]=w; z=t[z]; w=t[w]
        return True
    for p,q in f.items():
        if Ln[cid[p]]!=Ln[cid[q]]: return False
        if not set_cycle(cid[p],p,q): return False
        if cid[p] in pi and pi[cid[p]]!=cid[q]: return False
        pi[cid[p]]=cid[q]
    if len(set(pi.values()))!=len(pi): return False
    closed=set()
    def rec():
        # find an open key
        for C in list(pi):
            if C in closed: continue
            Cp=pi[C]
            if Cp==C:
                if Ln[C]%3==0: return False
                k=pow(3,-1,Ln[C]); z=cyc[C][0]
                if xmap[z]!=tpow(z,k): return False
                closed.add(C); r=rec(); closed.discard(C); return r
            inv_pi={v:k_ for k_,v in pi.items()}
            if Cp in pi:
                Cpp=pi[Cp]
                if Cpp==Cp or Cpp==C: return False
                if Cpp in pi:
                    if pi[Cpp]!=C: return False
                    z=cyc[C][0]
                    if xmap[xmap[xmap[z]]]!=t[z]: return False
                    closed.update([C,Cp,Cpp]); r=rec(); closed.difference_update([C,Cp,Cpp]); return r
                else:
                    if C in inv_pi: return False
                    z=cyc[C][0]; w=xmap[xmap[z]]
                    saved=dict(xmap); pi[Cpp]=C
                    ok=set_cycle(Cpp,w,t[z]) and rec()
                    del pi[Cpp]; xmap.clear(); xmap.update(saved); return ok
            else:
                # Cp image unknown: branch
                z=cyc[C][0]; w=xmap[z]
                if C in inv_pi:
                    D=inv_pi[C]
                    if D==Cp or Ln[D]!=Ln[C]: return False
                    # x on Cp forced: x(w)=x_D^{-1}(t z): need inverse of x on D
                    xinvD={xmap[u]:u for u in cyc[D]}
                    if t[z] not in xinvD: return False
                    saved=dict(xmap); pi[Cp]=D
                    ok=set_cycle(Cp,w,xinvD[t[z]]) and rec()
                    del pi[Cp]; xmap.clear(); xmap.update(saved); return ok
                cands=[D for D in range(len(cyc)) if D not in pi and D not in inv_pi and D!=C and D!=Cp and Ln[D]==Ln[C]]
                for D in cands:
                    for w2 in cyc[D]:
                        saved=dict(xmap); pi[Cp]=D; pi[D]=C
                        ok=set_cycle(Cp,w,w2) and set_cycle(D,w2,t[z]) and rec()
                        del pi[Cp]; del pi[D]; xmap.clear(); xmap.update(saved)
                        if ok: return True
                return False
        # all closed: remaining cycles
        rem=Counter(Ln[i] for i in range(len(cyc)) if i not in pi)
        return all(v%3==0 for l,v in rem.items() if l%3==0)
    return rec()


def support_decide(a,c):
    N=len(a); P=[p for p in range(N) if a[p]!=p]
    for img in itertools.permutations(range(N),len(P)):
        f=dict(zip(P,img))
        b=list(range(N))
        for p in P: b[f[p]]=f[a[p]]
        b=tuple(b)
        if len(set(b))!=N: continue
        t=comp(inv(b),c)
        if cube_root_with_constraints(t,f): return True
    return False


def stats(N):
    P=np.array(list(itertools.permutations(range(N))),dtype=np.int8)  # rows = x
    byclass=defaultdict(list)
    for p in map(tuple,P): byclass[ctype(p)].append(p)
    XX=np.take_along_axis(P,P,axis=1)  # x(x(j))
    tot_pairs=0; solvable=0; hist=Counter(); fact=math.factorial(N)
    for lam,members in byclass.items():
        a=np.array(members[0]); size=len(members)
        AXX=a[XX]                       # a(x(x(j)))
        C=np.take_along_axis(P,AXX,axis=1)  # x(a(x(x(j))))
        # encode rows to ints
        key=np.zeros(len(C),dtype=np.int64)
        for j in range(N): key=key*N+C[:,j]
        u,cnt=np.unique(key,return_counts=True)
        solvable+=size*len(u); tot_pairs+=size*fact
        h=Counter(cnt.tolist())
        for k,v in h.items(): hist[k]+=size*v
        hist[0]+=size*(fact-len(u))
    return solvable/tot_pairs, {k:hist[k]/tot_pairs for k in (0,1,2,3)}


def rigid_census2(N):
    P=np.array(list(itertools.permutations(range(N))),dtype=np.int8)
    XX=np.take_along_axis(P,P,axis=1); X3=np.take_along_axis(P,XX,axis=1)
    idn=np.arange(N,dtype=np.int8)
    ord3=np.all(X3==idn,axis=1) & ~np.all(P==idn,axis=1)
    byclass=defaultdict(list)
    for p in map(tuple,P): byclass[ctype(p)].append(p)
    solv=rig=uniq=0
    for lam,members in byclass.items():
        a=np.array(members[0]); size=len(members)
        C=np.take_along_axis(P,a[XX],axis=1)
        key=np.zeros(len(C),dtype=np.int64)
        for j in range(N): key=key*N+C[:,j]
        order=np.argsort(key,kind='stable'); ks=key[order]; o3=ord3[order]
        starts=np.flatnonzero(np.r_[True,ks[1:]!=ks[:-1]]); ends=np.r_[starts[1:],len(ks)]
        allo3=np.logical_and.reduceat(o3,starts)
        solv+=size*len(starts); rig+=size*allo3.sum(); uniq+=size*(allo3&((ends-starts)==1)).sum()
    return solv,int(rig),round(rig/solv,4),round(uniq/rig,3)

if __name__=="__main__":
    bad=[];tot=0
    for N in range(2,9):
        a=tuple((j+1)%N for j in range(N))
        for d in range(1,N):
            c=list(range(N)); c[0]=d; c[d]=0; c=tuple(c); tot+=1
            if P3count(a,c,PERMS[N])!=formula(N,d): bad.append((N,d))
    print('transposition',tot,bad)
    for (L,n) in [(2,4),(4,3),(5,3),(7,3),(4,4),(5,4),(2,6)]:
        els=wreath_elements(L,n); cubes={comp(p,comp(p,p)) for p,_,_ in els}
        print('cent',(L,n),len(els),sum((p in cubes)!=cent_criterion(L,n,rho,r) for p,rho,r in els))
    for N in range(1,7):
        perms=PERMS[N]; invols=[x for x in perms if comp(x,x)==tuple(range(N))]; reps={}
        for p in perms: reps.setdefault(ctype(p),p)
        tot=bad=0
        for lam,a_ in reps.items():
            for c_ in {comp(g,comp(a_,inv(g))) for g in perms}:
                bf=any(comp(x,comp(a_,inv(x)))==c_ for x in invols); tot+=1; bad+=bf!=conj2_criterion(a_,c_)
        print('conj2',N,tot,bad)
    a=tuple((j+1)%7 for j in range(7)); c1=(0,4,5,1,3,6,2); c2=(0,2,5,6,3,1,4)
    G1=gen_group([a,c1]); phi={g:evalword(w,[a,c2]) for g,w in G1.items()}
    print('psl',P3count(a,c1,PERMS[7]),P3count(a,c2,PERMS[7]),len(G1),
          all(phi[comp(g,h)]==comp(phi[g],phi[h]) for g in G1 for h in G1),all(ctype(g)==ctype(phi[g]) for g in G1))
    for N in range(1,8):
        cubes=Counter(comp(x,comp(x,x)) for x in PERMS[N]); reps={}
        for y in PERMS[N]: reps.setdefault(ctype(y),y)
        print('rho3',N,sum(cubes[y]!=rho3_formula(nu) for nu,y in reps.items()))
    for (L,k,p) in [(9,3,3),(8,3,3),(7,3,3),(9,2,3),(4,3,2),(5,2,2),(8,2,2),(5,3,3)]:
        units=[u for u in range(1,L) if math.gcd(u,L)==1]; bad=0
        for _ in range(40):
            m=[random.choice(units) for _ in range(k)]; e=tuple(range(k*L))
            sols=[x for x,s,t in conjugators(L,k,m) if perm_power(x,p)==e]
            bad+=bool(sols)!=crit(L,k,m,p)
        print('affinep',(L,k,p),bad)
    print('main cases bad',[(B,T) for B in range(1,13) for T in range(12*B+1,12*B+4) if main_case(B,T)])
    print('conjp arith',[conjp_arith(3,B) for B in (1,2,3,5)],[conjp_arith(5,B) for B in (1,2,3)])
    for N in range(4,10): print('stats',N,stats(N))
    for N in (5,6,7):
        tot=bad=0
        for s in (0,2,3,4):
            for _ in range(30):
                pts=random.sample(range(N),s); a=list(range(N))
                if s:
                    while True:
                        sh=pts[:]; random.shuffle(sh)
                        if all(x!=y for x,y in zip(pts,sh)): break
                    for p,q in zip(pts,sh): a[p]=q
                a=tuple(a); c=random.choice(PERMS[N]); tot+=1
                bad+=(P3count(a,c,PERMS[N])>0)!=support_decide(a,c)
        print('support algorithm',N,tot,bad)
    for N in (5,6,7,8): print('rigid',N,rigid_census2(N))
