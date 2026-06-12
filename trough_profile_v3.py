"""
ПАСЬЯНС-1c: ложбина, мягкий зонд.
g=0.12, клиппинг V_eff >= 0.6 (кинки выживают), контроль целостности:
- E/E0 (бракуем при |E/E0-1|>0.05)
- счёт кинков (должно быть ровно 2 пика градиента > 0.3*max)
Скан sep=5..26 шаг ~2, длинное окно.
"""
import numpy as np
L=300; dx=0.25; Nx=int(L/dx); c=1.0; dt=0.2*dx/c
x=np.arange(Nx)*dx
def prof(x0,s): return 4*np.arctan(np.exp(s*(x-x0)))
def lap_neumann(th):
    l=np.empty_like(th)
    l[1:-1]=(th[2:]-2*th[1:-1]+th[:-2])/dx**2
    l[0]=(th[1]-2*th[0]+th[1])/dx**2
    l[-1]=(th[-2]-2*th[-1]+th[-2])/dx**2
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
    return out, len(idx)
def energy(th,tht):
    gx=np.gradient(th,dx)
    return (0.5*tht**2+0.5*c**2*gx**2+(1-np.cos(th))).sum()*dx
def run(sep0,g,steps=24000):
    x1=L/2-sep0/2; x2=L/2+sep0/2
    th=prof(x1,+1)+prof(x2,-1); tht=np.zeros_like(th)
    E0=energy(th,tht); s0=None
    for t in range(steps):
        gx=np.gradient(th,dx)
        E=0.5*tht**2+0.5*c**2*gx**2+(1-np.cos(th))
        Veff=np.clip(1.0-g*E,0.6,1.0) if g>0 else np.ones_like(th)
        a=c**2*lap_neumann(th)-Veff*np.sin(th)
        tht=tht+dt*a; th=th+dt*tht
        if t==200:
            cs,n=centers_subpix(th); s0=cs[1]-cs[0] if n==2 else None
    cs,n=centers_subpix(th)
    Ef=energy(th,tht)
    ok = (n==2) and (s0 is not None) and abs(Ef/E0-1)<0.05
    s1=cs[1]-cs[0] if n==2 else float('nan')
    return s0,s1,Ef/E0,ok
print("=== g=0.12, клиппинг>=0.6 ===")
print(f"{'sep':>4} {'s0':>7} {'s1':>7} {'d':>7} {'E/E0':>7} {'целост':>7} {'знак':>8}")
for sep in [5,7,9,11,13,15,17,19,21,24,26]:
    s0,s1,er,ok=run(sep,0.12)
    if s0 is None: s0=float('nan')
    d=s1-s0
    sign='ПРИТЯЖ' if d<-0.3 else 'ОТТАЛК' if d>0.3 else 'нейтр'
    print(f"{sep:>4} {s0:>7.2f} {s1:>7.2f} {d:>+7.2f} {er:>7.3f} {'OK' if ok else 'БРАК':>7} {sign if ok else '--':>8}")
