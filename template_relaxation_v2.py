"""
ПАСЬЯНС-2b: релаксация с прогревом.
1) Прогрев 20000 шагов БЕЗ распадов. 2) Замер фоновой рассинхронии тяжёлых.
3) Пороги dE = фон + к*сигма_фона. 4) Окно наблюдения распадов.
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

def experiment(kappa_list, M=60, N=3000, deg=6, K=1.0, dt=0.02, warm=20000, obs=80000, seed=0):
    na,valid,cnt=build(N,deg,seed)
    rng=np.random.RandomState(seed+900)
    th0=rng.uniform(0,2*np.pi,N)
    heavy0=rng.choice(N,M,replace=False).tolist()
    def step(th,om):
        tn=np.where(valid,th[na],0.0)
        coup=np.where(valid,np.sin(tn-th[:,None]),0.0).sum(1)/cnt
        return th+dt*(om+K*coup)
    def mis_of(th,h):
        return np.abs(np.where(valid[h],np.sin(th[na[h]]-th[h]),0.0)).sum()/cnt[h]
    # прогрев
    om=np.ones(N); om[heavy0]=2.0
    th=th0.copy()
    for t in range(warm): th=step(th,om)
    # фон рассинхронии тяжёлых (по окну 2000 шагов)
    samples=[]
    th_b=th.copy()
    for t in range(2000):
        th_b=step(th_b,om)
        if t%100==0:
            samples += [mis_of(th_b,h) for h in heavy0]
    mu,sd=np.mean(samples),np.std(samples)
    print(f"фон рассинхронии тяжёлых: {mu:.4f} ± {sd:.4f}")
    results={}
    for kap in kappa_list:
        dE=mu+kap*sd
        omr=np.ones(N); omr[heavy0]=2.0
        thr=th.copy()
        heavy=set(heavy0)
        times=[]
        for t in range(obs):
            thr=step(thr,omr)
            if t%25==0 and heavy:
                for h in list(heavy):
                    if mis_of(thr,h)>dE:
                        omr[h]=1.0; heavy.discard(h); times.append(t*dt)
        results[kap]=(dE,times,len(heavy))
        dec=len(times)
        med=np.median(times) if dec else float('nan')
        print(f"kappa={kap}: dE={dE:.4f}, распалось {dec}/{M}, медианное время {med:8.0f}, выжило {len(heavy)}")
    return results, M

print("=== Релаксация шаблонов: пороги от прогретого фона ===")
res,M=experiment([1.5, 2.5, 3.5])
print()
# Кривая выживания и экспонента для среднего порога
kap=2.5
dE,times,_=res[kap]
if len(times)>=10:
    times=np.sort(np.array(times))
    surv_t=times
    surv_n=M-np.arange(1,len(times)+1)
    mask=surv_n>0
    A=np.vstack([surv_t[mask],np.ones(mask.sum())]).T
    coef,_,_,_=np.linalg.lstsq(A,np.log(surv_n[mask]/M),rcond=None)
    tau=-1/coef[0] if coef[0]<0 else float('inf')
    pred=A@coef
    lnS=np.log(surv_n[mask]/M)
    r2=1-((lnS-pred)**2).sum()/((lnS-lnS.mean())**2).sum()
    print(f"kappa=2.5: выживание ~ exp(-t/tau), tau={tau:.0f}, R²={r2:.4f} (по {mask.sum()} распадам)")
