"""
ПАСЬЯНС-3: H-тест — ведёт ли себя распределённость возмущения как энтропия.
ЗАМКНУТАЯ система: фазовая сеть ВТОРОГО порядка (с инерцией) без затухания и накачки —
энергия сохраняется. Старт: вся рассинхрония локализована (горб фазовой скорости на k узлах).
Метрика: энтропия Шеннона S = -sum p_i ln p_i, p_i = E_i/sum(E), E_i = 0.5*v_i^2 + связевая часть.
Гипотезы (свойства энтропии): (1) S растёт монотонно (с флуктуациями) к плато;
(2) плато ~ ln N (равнораспределение); (3) энергия сохранена; (4) спонтанного возврата к старту нет.
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

N=2000; deg=6; K=1.0; dt=0.01
na,valid,cnt=build(N,deg,0)
rng=np.random.RandomState(77)
th=np.zeros(N)            # идеальный синхрон фаз
v=np.zeros(N)             # скорости (отклонения темпа)
hot=rng.choice(N,10,replace=False)
v[hot]=3.0                # вся рассинхрония в 10 узлах из 2000
def forces(th):
    tn=np.where(valid,th[na],0.0)
    return K*np.where(valid,np.sin(tn-th[:,None]),0.0).sum(1)/cnt
def energy_parts(th,v):
    kin=0.5*v**2
    tn=np.where(valid,th[na],0.0)
    pot=K*np.where(valid,(1-np.cos(tn-th[:,None]))*0.5,0.0).sum(1)/cnt  # делим связь пополам
    return kin+pot
def shannon(E):
    p=E/E.sum()
    p=p[p>1e-15]
    return -(p*np.log(p)).sum()
E0=energy_parts(th,v).sum()
print(f"N={N}, старт: вся энергия в 10 узлах; S_max теоретич ~ ln(N)={np.log(N):.2f}")
print(f"{'t':>7} {'S':>7} {'E/E0':>7}")
S_hist=[]
for t in range(120001):
    a=forces(th)
    v=v+dt*a; th=th+dt*v   # симплектический Эйлер: энергия дрейфует слабо
    if t%8000==0:
        E=energy_parts(th,v)
        S=shannon(E)
        S_hist.append(S)
        print(f"{t:>7} {S:>7.3f} {E.sum()/E0:>7.4f}")
# монотонность (по сглаженной истории)
inc=sum(1 for i in range(1,len(S_hist)) if S_hist[i]>=S_hist[i-1]-0.02)
print(f"\nмонотонных шагов (с допуском флуктуаций): {inc}/{len(S_hist)-1}")
print(f"стартовая S={S_hist[0]:.3f}, финальная S={S_hist[-1]:.3f}, ln N={np.log(N):.3f}")
