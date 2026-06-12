"""
ПАСЬЯНС-3b: H-тест на СИММЕТРИЗОВАННОМ графе (взаимные связи => силы потенциальны,
энергия сохраняется => честная замкнутая система).
"""
import numpy as np
def build_sym(N, deg, seed, box=30.0):
    rng=np.random.RandomState(seed)
    pos=rng.uniform(0,box,(N,3))
    cell=box/8; grid={}
    for i,p in enumerate(pos):
        grid.setdefault(tuple((p//cell).astype(int)),[]).append(i)
    pairs=set()
    for i,p in enumerate(pos):
        key=(p//cell).astype(int); cand=[]
        for a in(-1,0,1):
            for b in(-1,0,1):
                for cc in(-1,0,1):
                    cand+=grid.get((key[0]+a,key[1]+b,key[2]+cc),[])
        cand=[j for j in cand if j!=i]
        d=np.sqrt(((pos[cand]-p)**2).sum(1))
        order=np.argsort(d)[:deg]
        for k in order: pairs.add((min(i,cand[k]),max(i,cand[k])))
    adj=[[] for _ in range(N)]
    for i,j in pairs:
        adj[i].append(j); adj[j].append(i)
    maxd=max(len(a) for a in adj)
    na=np.full((N,maxd),-1)
    for i,a in enumerate(adj): na[i,:len(a)]=a
    valid=na>=0
    return na,valid,list(pairs)

N=2000; deg=6; K=1.0; dt=0.01
na,valid,pairs=build_sym(N,deg,0)
pairs=np.array(pairs)
rng=np.random.RandomState(77)
th=np.zeros(N); v=np.zeros(N)
hot=rng.choice(N,10,replace=False)
v[hot]=3.0
def forces(th):
    tn=np.where(valid,th[na],0.0)
    return K*np.where(valid,np.sin(tn-th[:,None]),0.0).sum(1)   # БЕЗ нормировки на cnt: парные силы симметричны
def total_energy(th,v):
    kin=0.5*(v**2).sum()
    pot=K*(1-np.cos(th[pairs[:,1]]-th[pairs[:,0]])).sum()
    return kin+pot
def node_energy(th,v):
    kin=0.5*v**2
    dE=K*(1-np.cos(th[pairs[:,1]]-th[pairs[:,0]]))
    pot=np.zeros(N)
    np.add.at(pot,pairs[:,0],0.5*dE); np.add.at(pot,pairs[:,1],0.5*dE)
    return kin+pot
def shannon(E):
    p=E/E.sum(); p=p[p>1e-15]
    return -(p*np.log(p)).sum()
E0=total_energy(th,v)
print(f"N={N}, ln N={np.log(N):.3f}, E0={E0:.2f}")
print(f"{'t':>7} {'S':>7} {'E/E0':>8}")
S_hist=[]
# leapfrog
a=forces(th)
for t in range(120001):
    v=v+0.5*dt*a
    th=th+dt*v
    a=forces(th)
    v=v+0.5*dt*a
    if t%8000==0:
        S=shannon(node_energy(th,v))
        S_hist.append(S)
        print(f"{t:>7} {S:>7.3f} {total_energy(th,v)/E0:>8.4f}")
inc=sum(1 for i in range(1,len(S_hist)) if S_hist[i]>=S_hist[i-1]-0.02)
print(f"\nмонотонных (с допуском): {inc}/{len(S_hist)-1}; старт {S_hist[0]:.3f} -> финал {S_hist[-1]:.3f} (ln N={np.log(N):.3f})")
