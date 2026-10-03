import itertools, random, json, time
def compose(p,q): return tuple(p[q[j]] for j in range(len(p)))
def build(L,ms):
    k=len(ms)
    a=tuple(i*L+(z+1)%L for i in range(k) for z in range(L))
    c=tuple(i*L+(z+ms[i])%L for i in range(k) for z in range(L))
    return a,c
def has_order3_conjugator(a,c,L,k):
    N=L*k; idn=tuple(range(N))
    for s in itertools.permutations(range(k)):
        for base in itertools.product(range(L),repeat=k):
            x=[0]*N
            for i in range(k):
                p=i*L; q=s[i]*L+base[i]
                for _ in range(L):
                    x[p]=q; p=a[p]; q=c[q]
            x=tuple(x)
            if compose(x,compose(x,x))==idn: return True
    return False
def criterion(L,ms):
    k=len(ms)
    for s in itertools.permutations(range(k)):
        if any(s[s[s[i]]]!=i for i in range(k)): continue
        ok=True
        for i in range(k):
            if s[i]==i:
                if pow(ms[i],3,L)!=1: ok=False;break
            elif i<s[i] and i<s[s[i]]:
                if (ms[i]*ms[s[i]]*ms[s[s[i]]])%L!=1: ok=False;break
        if ok: return True
    return False
if __name__=='__main__':
    rng=random.Random(3); agree=tot=pos=0; t0=time.time()
    for L,k in ((5,3),(7,3),(5,4),(11,3)):
        for _ in range(40):
            pass
        for _ in range(40):
            ms=[rng.randrange(1,L) for _ in range(k)]
            a,c=build(L,ms); ex=has_order3_conjugator(a,c,L,k); cr=criterion(L,ms)
            agree+=(ex==cr); tot+=1; pos+=ex
    print('criterion vs exhaustive conjugators:',agree,'/',tot,'positives',pos,'%.0fs'%(time.time()-t0))
    L=11; g=2
    out={}
    for name,es in (('yes',[1,2,7]),('no',[1,2,8]),('no_sum15',[5,5,5]),('no2',[5,5,6])):
        ms=[pow(g,e,L) for e in es]; a,c=build(L,ms)
        ex=has_order3_conjugator(a,c,L,len(ms)); cr=criterion(L,ms); out[name]=(es,ex,cr)
        print(name,es,'D0 exhaustive:',ex,'criterion:',cr)
    json.dump(dict(agree=agree,total=tot,positives=pos,n3dm_toy=out),open('verify_lemmaC.json','w'))
    