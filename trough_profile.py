"""
ПАСЬЯНС-1: ложбина как комбинация двух градиентов.
Кинк-антикинк в sine-Gordon с ДОБАВКОЙ задержки: локальный провал темпа
(коэффициент при потенциале снижен там, где энергия поля высока — узлы с массой тикают медленнее:
V_eff = 1 - g*E_local, клиппинг снизу). Это вводит вторую, медленно спадающую компоненту.
Мерим знак НАЧАЛЬНОГО ускорения пары из покоя при разных sep:
  чистый SG: притяжение на всех sep (монотонная яма в нуле -> аннигиляция/бризер).
  с задержкой: если возникает ЛОЖБИНА на конечном r*, знак силы меняется:
  sep < r* -> отталкивание (выталкивает ко дну), sep > r* -> притяжение (затягивает ко дну).
Контроль: g=0 (чистый SG) — притяжение везде.
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

def local_energy(th,tht):
    gx=np.gradient(th,dx)
    return 0.5*tht**2+0.5*c**2*gx**2+(1-np.cos(th))

def initial_accel_sep(sep0, g, probe_steps=600):
    x1=L/2-sep0/2; x2=L/2+sep0/2
    th=prof(x1,+1)+prof(x2,-1)   # кинк-антикинк
    tht=np.zeros_like(th)
    # короткая эволюция; знак силы = знак d(sep)/dt после probe_steps
    def centers(th):
        gr=np.abs(np.gradient(th,dx))
        idx=[i for i in range(2,Nx-2) if gr[i]>gr[i-1] and gr[i]>gr[i+1] and gr[i]>0.3*gr.max()]
        idx=sorted(idx,key=lambda i:-gr[i])[:2]
        return sorted(x[i] for i in idx)
    s_start=None
    for t in range(probe_steps):
        E=local_energy(th,tht)
        Veff=np.clip(1.0-g*E, 0.2, 1.0)   # задержка: высокая энергия -> медленный тик -> слабее потенциал
        a=c**2*lap_neumann(th)-Veff*np.sin(th)
        tht=tht+dt*a; th=th+dt*tht
        if t==50:
            cs=centers(th); s_start=cs[1]-cs[0] if len(cs)==2 else sep0
    cs=centers(th)
    s_end=cs[1]-cs[0] if len(cs)==2 else 0.0
    return s_end-s_start  # >0 разошлись (отталкивание), <0 сошлись (притяжение)

print("=== Контроль: чистый SG (g=0) — ожидаем притяжение на всех sep ===")
print(f"{'sep':>5} {'d(sep)':>9} {'знак':>14}")
for sep in [6,10,14,18,24]:
    d=initial_accel_sep(sep, 0.0)
    print(f"{sep:>5} {d:>9.3f} {'ПРИТЯЖЕНИЕ' if d<-0.05 else 'ОТТАЛКИВАНИЕ' if d>0.05 else 'нейтрально':>14}")
print()
print("=== С задержкой (g=0.35): ищем смену знака = дно ложбины ===")
for sep in [4,6,8,10,14,18,24]:
    d=initial_accel_sep(sep, 0.35)
    print(f"{sep:>5} {d:>9.3f} {'ПРИТЯЖЕНИЕ' if d<-0.05 else 'ОТТАЛКИВАНИЕ' if d>0.05 else 'нейтрально':>14}")
print()
print("Ложбина: ОТТАЛКИВАНИЕ на малых sep + ПРИТЯЖЕНИЕ на больших = минимум между ними")
