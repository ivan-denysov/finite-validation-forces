"""
ПАСЬЯНС-1b: ложбина — честное окно (как в kink_force_static_v2: 24000 шагов)
и субпиксельные центры (параболическая интерполяция пика градиента).
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
    return out
def run(sep0,g,steps=24000):
    x1=L/2-sep0/2; x2=L/2+sep0/2
    th=prof(x1,+1)+prof(x2,-1); tht=np.zeros_like(th)
    s0=None
    for t in range(steps):
        gx=np.gradient(th,dx)
        E=0.5*tht**2+0.5*c**2*gx**2+(1-np.cos(th))
        Veff=np.clip(1.0-g*E,0.2,1.0) if g>0 else 1.0
        a=c**2*lap_neumann(th)-Veff*np.sin(th)
        tht=tht+dt*a; th=th+dt*tht
        if t==200:
            cs=centers_subpix(th); s0=cs[1]-cs[0] if len(cs)==2 else sep0
    cs=centers_subpix(th)
    s1=cs[1]-cs[0] if len(cs)==2 else 0.0
    return s0,s1
print("=== Контроль g=0 (ожидаем притяжение, эталон 14->9.5) ===")
for sep in [10,14,18]:
    s0,s1=run(sep,0.0)
    d=s1-s0
    print(f"sep={sep:>3}: {s0:.2f} -> {s1:.2f}  d={d:+.2f}  {'ПРИТЯЖ' if d<-0.3 else 'ОТТАЛК' if d>0.3 else 'нейтр'}")
print()
print("=== g=0.35: профиль знака ===")
for sep in [5,7,9,12,15,19,24]:
    s0,s1=run(sep,0.35)
    d=s1-s0
    print(f"sep={sep:>3}: {s0:.2f} -> {s1:.2f}  d={d:+.2f}  {'ПРИТЯЖ' if d<-0.3 else 'ОТТАЛК' if d>0.3 else 'нейтр'}")
