import itertools, random, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_lemmaC import build, has_order3_conjugator, criterion
def is_prime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True
def choose_L(B):
    L=12*B+2
    while not (is_prime(L) and L%3==2): L+=1
    return L
worst=max((choose_L(B)/(12*B+2),B) for B in range(1,10001))
print('max L/(12B+2) over B<=10^4:',worst)
ok=True
for B in range(1,41):
    for T in (12*B+1, 12*B+7, 30*B):
        M=3*B
        E={'W':[s for s in range(1,B+1)],'X':[s+M for s in range(1,B+1)],'Y':[s+T-M-B for s in range(1,B+1)]}
        for types in itertools.combinations_with_replacement('WXY',3):
            for vals in itertools.product(*[range(1,B+1) for _ in types]):
                e=sum(E[t][v-1] for t,v in zip(types,vals))
                zero=(e%T==0); should=(types==('W','X','Y') and sum(vals)==B)
                if zero!=should: ok=False; print('MISMATCH',B,T,types,vals)
print('case analysis exhaustive B<=40:',ok)
def reduce_n3dm(W,X,Y,B):
    L=choose_L(B); T=L-1; M=3*B
    g=next(g for g in range(2,L) if len({pow(g,e,L) for e in range(T)})==T)
    es=[s for s in W]+[s+M for s in X]+[s+T-M-B for s in Y]
    return L,[pow(g,e,L) for e in es]
def n3dm_bruteforce(W,X,Y,B):
    q=len(W)
    return any(all(W[i]+X[px[i]]+Y[py[i]]==B for i in range(q)) for px in itertools.permutations(range(q)) for py in itertools.permutations(range(q)))
rng=random.Random(7); agree=0; tot=0; res=[]
for q in (1,2,3):
    for _ in range(60):
        B=rng.randrange(3,9); W=[rng.randrange(1,B) for _ in range(q)]; X=[rng.randrange(1,B) for _ in range(q)]; Y=[rng.randrange(1,B) for _ in range(q)]
        truth=n3dm_bruteforce(W,X,Y,B); L,ms=reduce_n3dm(W,X,Y,B)
        d0=criterion(L,ms); agree+=(truth==d0); tot+=1
        if q==1 and len(res)<8:
            a,c=build(L,ms); ex=has_order3_conjugator(a,c,L,3); res.append(dict(W=W,X=X,Y=Y,B=B,L=L,N=3*L,n3dm=truth,criterion=d0,exhaustive=ex))
print('reduction vs N3DM brute force (criterion):',agree,'/',tot)
print('q=1 exhaustive conjugator checks (n3dm, criterion, exhaustive, N):',[(r['n3dm'],r['criterion'],r['exhaustive'],r['N']) for r in res])
json.dump(dict(prime_bound_worst=worst,case_analysis_ok=ok,criterion_agree=agree,total=tot,exhaustive_q1=res),open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'verify_reduction.json'),'w'),indent=1)
# planted YES instances at q=1 and q=2 (q=2: k=6, exhaustive over L^6*720 conjugators infeasible -> criterion only; q=1 exhaustive)
yes_res=[]
for (W,X,Y,B) in (([1],[1],[2],4),([2],[3],[1],6),([1],[2],[2],5)):
    L,ms=reduce_n3dm(W,X,Y,B); a,c=build(L,ms)
    yes_res.append(dict(W=W,X=X,Y=Y,B=B,L=L,N=3*L,n3dm=n3dm_bruteforce(W,X,Y,B),criterion=criterion(L,ms),exhaustive=has_order3_conjugator(a,c,L,3)))
print('planted YES q=1 (n3dm, criterion, exhaustive, N):',[(r['n3dm'],r['criterion'],r['exhaustive'],r['N']) for r in yes_res])
d=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'verify_reduction.json'))); d['exhaustive_q1_planted_yes']=yes_res
json.dump(d,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'verify_reduction.json'),'w'),indent=1)
