"""
ПАСЬЯНС-2c: пороги в ХВОСТЕ распределения фона (квантили) — распад как редкое событие.
tau(квантиль) = иерархия времён жизни от высоты барьера.
"""
import numpy as np
def build(N, deg, seed, box=30.0):
    rng=np.random.RandomState(seed)
    pos=rng.uniform(0,box,(N,3))
    cell=box/8; grid={}
    for i,p in enumerate(pos):
        grid.setdefault(tuple((p//cell).astype(int)),[]).append(i)
    nbrs=[]
    for i,p in enumerate(pos):
        key=(p//cell).astype(int); cand=[]
        for a in(-1,0,1):
            for b in(-1,0,1):
                for cc in(-1,0,1):
                    cand+=grid.get((key[0]+a,key[1]+b,key[2]+cc),[])
        cand=[j for j in cand if j!=i]
        d=np.sqrt(((pos[cand]-p)**2).sum(1))
        order=np.argsort(d)[:deg]
        nbrs.append([cand[k] for k in order])
    maxd=max(len(n) for n in nbrs)
    na=np.full((N,maxd),-1)
    for i,n in enumerate(nbrs): na[i,:len(n)]=n
    valid=na>=0; cnt=valid.sum(1); cnt[cnt==0]=1
    return na,valid,cnt

N=3000; deg=6; K=1.0; dt=0.02; M=60
na,valid,cnt=build(N,deg,0)
rng=np.random.RandomState(900)
th=rng.uniform(0,2*np.pi,N)
heavy0=rng.choice(N,M,replace=False).tolist()
om=np.ones(N); om[heavy0]=2.0
def step(th,om):
    tn=np.where(valid,th[na],0.0)
    coup=np.where(valid,np.sin(tn-th[:,None]),0.0).sum(1)/cnt
    return th+dt*(om+K*coup)
def mis(th,h):
    return np.abs(np.where(valid[h],np.sin(th[na[h]]-th[h]),0.0)).sum()/cnt[h]
for t in range(20000): th=step(th,om)
# фон: длинное окно для хвоста
samples=[]
thb=th.copy()
for t in range(6000):
    thb=step(thb,om)
    if t%20==0: samples += [mis(thb,h) for h in heavy0]
samples=np.array(samples)
qs={0.95:np.quantile(samples,0.95), 0.99:np.quantile(samples,0.99), 0.999:np.quantile(samples,0.999)}
print("фон: mean=%.3f, q95=%.4f, q99=%.4f, q999=%.4f, max=%.4f"%(samples.mean(),qs[0.95],qs[0.99],qs[0.999],samples.max()))
print()
print("=== tau от высоты барьера (квантиль фона) ===")
results={}
for q,dE in qs.items():
    omr=np.ones(N); omr[heavy0]=2.0
    thr=th.copy(); heavy=set(heavy0); times=[]
    for t in range(80000):
        thr=step(thr,omr)
        if t%25==0 and heavy:
            for h in list(heavy):
                if mis(thr,h)>dE:
                    omr[h]=1.0; heavy.discard(h); times.append(t*dt)
    dec=len(times)
    med=np.median(times) if dec else float('inf')
    results[q]=(times,dec)
    print(f"q={q}: dE={dE:.4f}, распалось {dec}/{M}, медиана {med:8.1f}")
print()
# экспонента выживания для q с лучшей статистикой
best=max(results, key=lambda q: min(results[q][1], M-1) if results[q][1]<M else M//2)
times,dec=results[best]
if dec>=10:
    times=np.sort(np.array(times)); surv=M-np.arange(1,dec+1)
    mask=surv>0
    A=np.vstack([times[mask],np.ones(mask.sum())]).T
    coef,_,_,_=np.linalg.lstsq(A,np.log(surv[mask]/M),rcond=None)
    tau=-1/coef[0]
    pred=A@coef; lnS=np.log(surv[mask]/M)
    r2=1-((lnS-pred)**2).sum()/((lnS-lnS.mean())**2).sum()
    print(f"кривая выживания (q={best}): tau={tau:.0f}, R²={r2:.4f} ({dec} распадов)")
