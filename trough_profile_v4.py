"""
ПАСЬЯНС-1d: ложбина = kink-KINK (одноимённые, SG-отталкивание, хвост e^-r)
+ диффузионное поле задержки phi (гравитация: провал темпа, распределяемый сетью,
хвост e^{-r/lambda}, lambda=sqrt(D/gamma) >> 1).
V_eff = 1 - phi  (где провал темпа — там потенциал слабее: кинку энергетически
выгодно сидеть в провале -> эффективное притяжение к общему провалу).
Ожидание: ОТТАЛК на малых sep (короткий хвост сильнее), ПРИТЯЖ на больших (длинный хвост живёт),
смена знака = дно ложбины r*.
"""
import numpy as np
L=300; dx=0.25; Nx=int(L/dx); c=1.0; dt=0.2*dx/c
x=np.arange(Nx)*dx
D=0.5; gam=0.005; s_src=0.004   # lambda=sqrt(D/gam)=10; CFL diff: D*dt/dx^2=0.4<0.5
def prof(x0,s): return 4*np.arctan(np.exp(s*(x-x0)))
def lapn(f):
    l=np.empty_like(f)
    l[1:-1]=(f[2:]-2*f[1:-1]+f[:-2])/dx**2
    l[0]=(f[1]-2*f[0]+f[1])/dx**2
    l[-1]=(f[-2]-2*f[-1]+f[-2])/dx**2
    return l
def centers_subpix(th):
    g=np.abs(np.gradient(th,dx))
    idx=[i for i in range(2,Nx-2) if g[i]>g[i-1] and g[i]>g[i+1] and g[i]>0.3*g.max()]
    idx=sorted(idx,key=lambda i:-g[i])[:2]
    out=[]
    for i in sorted(idx):
        y0,y1,y2=g[i-1],g[i],g[i+1]
        d=0.5*(y0-y2)/(y0-2*y1+y2+1e-12)
        out.append(x[i]+d*dx)
    return out,len(idx)
def run(sep0, gmode, steps=24000):
    x1=L/2-sep0/2; x2=L/2+sep0/2
    th=prof(x1,+1)+prof(x2,+1)   # KINK-KINK
    tht=np.zeros_like(th)
    phi=np.zeros_like(th)
    s0=None
    for t in range(steps):
        gx=np.gradient(th,dx)
        E=0.5*tht**2+0.5*c**2*gx**2+(1-np.cos(th))
        if gmode:
            phi=phi+dt*(D*lapn(phi)+s_src*E-gam*phi)
            Veff=np.clip(1.0-phi,0.5,1.0)
        else:
            Veff=1.0
        a=c**2*lapn(th)-Veff*np.sin(th)
        tht=tht+dt*a; th=th+dt*tht
        if t==200:
            cs,n=centers_subpix(th); s0=cs[1]-cs[0] if n==2 else None
    cs,n=centers_subpix(th)
    s1=cs[1]-cs[0] if n==2 else float('nan')
    ok=(n==2) and (s0 is not None)
    return s0,s1,ok,phi.max()
print(f"lambda = {np.sqrt(D/gam):.1f}")
print("=== Контроль: kink-kink БЕЗ phi (эталон: отталкивание 14->17) ===")
for sep in [10,14]:
    s0,s1,ok,_=run(sep,False)
    d=s1-s0
    print(f"sep={sep:>3}: {s0:.2f}->{s1:.2f} d={d:+.2f} {'ОТТАЛК' if d>0.3 else 'ПРИТЯЖ' if d<-0.3 else 'нейтр'}")
print()
print("=== kink-kink + диффузионная задержка: ищем ложбину ===")
print(f"{'sep':>4} {'s0':>7} {'s1':>7} {'d':>7} {'phi_max':>8} {'знак':>8}")
for sep in [6,9,12,15,18,22,26,30]:
    s0,s1,ok,pm=run(sep,True)
    if s0 is None: s0=float('nan')
    d=s1-s0
    sign='ОТТАЛК' if d>0.3 else 'ПРИТЯЖ' if d<-0.3 else 'нейтр'
    print(f"{sep:>4} {s0:>7.2f} {s1:>7.2f} {d:>+7.2f} {pm:>8.3f} {sign if ok else 'БРАК':>8}")
