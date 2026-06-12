"""
ПАСЬЯНС-2: слабое как релаксация шаблонов.
Сеть Курамото; M "тяжёлых" носителей (drive=2.0). Тяжёлый шаблон пересобирается
в лёгкий (drive->1.0), когда локальная рассинхрония вокруг него (|фазовый дефицит|
относительно соседей) случайно превысит порог dE (цена перехода = глубина
несостоятельности, которую надо набрать флуктуацией).
Ожидания: (а) кривая выживания N(t) ~ exp(-t/tau); (б) tau растёт с dE.
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

def run(dE, M=60, N=3000, deg=6, K=1.0, dt=0.02, steps=60000, seed=0):
    na,valid,cnt=build(N,deg,seed)
    rng=np.random.RandomState(seed+900)
    th=rng.uniform(0,2*np.pi,N)
    heavy=set(rng.choice(N,M,replace=False).tolist())
    om=np.ones(N)
    for h in heavy: om[h]=2.0
    decay_times=[]
    alive=[]
    for t in range(steps):
        tn=np.where(valid,th[na],0.0)
        sino=np.sin(tn-th[:,None])
        coup=np.where(valid,sino,0.0).sum(1)/cnt
        th=th+dt*(om+K*coup)
        # локальная рассинхрония носителя: средний |sin(соседи - я)|
        if t%25==0 and heavy:
            for h in list(heavy):
                mis=np.abs(np.where(valid[h],np.sin(th[na[h]]-th[h]),0.0)).sum()/cnt[h]
                if mis>dE:
                    om[h]=1.0; heavy.discard(h)
                    decay_times.append(t*dt)
        if t%2000==0: alive.append((t*dt,len(heavy)))
    return alive, decay_times, M

print("=== (а) Кривая выживания при dE=0.62 ===")
alive,times,M=run(0.62)
for tt,n in alive[::3]: print(f"t={tt:7.0f}: живых {n}/{M}")
# фит экспоненты по кривой выживания
import numpy as np
ts=np.array([a[0] for a in alive if a[1]>0]); ns=np.array([a[1] for a in alive if a[1]>0],float)
if (ns<M).sum()>=4:
    mask=ns>0
    A=np.vstack([ts[mask],np.ones(mask.sum())]).T
    coef,_,_,_=np.linalg.lstsq(A,np.log(ns[mask]/M),rcond=None)
    tau=-1/coef[0] if coef[0]<0 else float('inf')
    pred=A@coef
    r2=1-((np.log(ns[mask]/M)-pred)**2).sum()/((np.log(ns[mask]/M)-np.log(ns[mask]/M).mean())**2).sum()
    print(f"экспонента: tau={tau:.0f}, R²={r2:.4f}")
print()
print("=== (б) tau(dE): иерархия времён жизни ===")
for dE in [0.55,0.62,0.70]:
    alive,times,M=run(dE, steps=60000)
    dec=len(times)
    if dec>=10:
        tau_est=np.mean(times) if dec>M*0.6 else (alive[-1][0]*M/max(dec,1))
        print(f"dE={dE}: распалось {dec}/{M}, средн. время распавшихся {np.mean(times):7.0f}, оценка tau~{tau_est:7.0f}")
    else:
        print(f"dE={dE}: распалось {dec}/{M} — слишком стабильны на этом окне")
