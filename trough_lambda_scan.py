"""
ПАСЬЯНС-1e: формула ложбины — r*(lambda).
Скан lambda = 5, 7, 10, 14, 20 (gam = D/lambda^2, D=0.5 фикс, CFL ok).
r* мерим сходимостью с двух сторон: старт sep=12 (сверху) и sep=5 (снизу),
равновесие = среднее s1 двух прогонов (контроль согласия < 0.3).
Гипотеза для фита: r* ~ a*ln(lambda)+b ИЛИ r* ~ c*lambda^p — увидим по данным.
"""
import numpy as np
L=300; dx=0.25; Nx=int(L/dx); c=1.0; dt=0.2*dx/c
x=np.arange(Nx)*dx
D=0.5; s_src=0.004
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
def run(sep0, gam, steps=24000):
    x1=L/2-sep0/2; x2=L/2+sep0/2
    th=prof(x1,+1)+prof(x2,+1); tht=np.zeros_like(th); phi=np.zeros_like(th)
    for t in range(steps):
        gx=np.gradient(th,dx)
        E=0.5*tht**2+0.5*c**2*gx**2+(1-np.cos(th))
        phi=phi+dt*(D*lapn(phi)+s_src*E-gam*phi)
        Veff=np.clip(1.0-phi,0.5,1.0)
        a=c**2*lapn(th)-Veff*np.sin(th)
        tht=tht+dt*a; th=th+dt*tht
    cs,n=centers_subpix(th)
    return cs[1]-cs[0] if n==2 else float('nan')
print(f"{'lambda':>7} {'r*(сверху)':>11} {'r*(снизу)':>10} {'r*':>7} {'согласие':>9}")
res=[]
for lam in [5.0,7.0,10.0,14.0,20.0]:
    gam=D/lam**2
    r_hi=run(12,gam); r_lo=run(5,gam)
    agree=abs(r_hi-r_lo)
    rstar=0.5*(r_hi+r_lo)
    res.append((lam,rstar,agree))
    print(f"{lam:>7.1f} {r_hi:>11.2f} {r_lo:>10.2f} {rstar:>7.2f} {agree:>9.2f}")
good=[(l,r) for l,r,a in res if a<0.5 and not np.isnan(r)]
if len(good)>=3:
    ls=np.array([g[0] for g in good]); rs=np.array([g[1] for g in good])
    # два фита
    A=np.vstack([np.log(ls),np.ones(len(ls))]).T
    c1,_,_,_=np.linalg.lstsq(A,rs,rcond=None)
    pred=A@c1; r2log=1-((rs-pred)**2).sum()/((rs-rs.mean())**2).sum()
    B=np.vstack([np.log(ls),np.ones(len(ls))]).T
    c2,_,_,_=np.linalg.lstsq(B,np.log(rs),rcond=None)
    pred2=B@c2; r2pow=1-((np.log(rs)-pred2)**2).sum()/((np.log(rs)-np.log(rs).mean())**2).sum()
    print(f"\nлог-фит:   r* = {c1[0]:.2f}*ln(lambda) + {c1[1]:.2f}   (R²={r2log:.4f})")
    print(f"степенной: r* ~ lambda^{c2[0]:.2f}                  (R²={r2pow:.4f})")
