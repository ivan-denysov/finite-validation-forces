# -*- coding: utf-8 -*-
"""
E.1 STRESS TEST (reviewer request): does the network-solves-Poisson result
survive at large frequency deficit, where sin saturates and the linearised
identity K*psi'' = Omega - omega no longer holds trivially?

At small d_omega the Gauss jump = d_omega/K is guaranteed by linearisation
(the run only confirms the linearisation is valid). The nontrivial question:
where does it break? Scan d_omega and report, per value:
  - locking (std of instantaneous frequencies; large => not locked)
  - Gauss jump vs the linear prediction d_omega/K
  - linear-vs-exponential tail fit
"""
import numpy as np

N = 801; K = 1.0; dt = 0.05; IC = N // 2

def relax(dw, T=15000):
    w = np.zeros(N); w[IC] = -dw; th = np.zeros(N)
    for _ in range(int(T / dt)):
        s = np.sin(np.diff(th))
        f = np.zeros(N); f[:-1] += K*s; f[1:] -= K*s
        th += dt * (w + f)
    # locking
    s = np.sin(np.diff(th)); f = np.zeros(N); f[:-1] += K*s; f[1:] -= K*s
    freq = w + f
    g = np.gradient(th - th[IC])
    jump = g[IC+2] - g[IC-2]
    r = np.arange(10, 121); gr = g[IC+10:IC+121]
    def r2(y, fp): return 1 - np.sum((y-fp)**2)/np.sum((y-np.mean(y))**2)
    r2lin = r2(gr, np.polyval(np.polyfit(r, gr, 1), r))
    ce = np.polyfit(r, np.log(np.abs(gr)), 1)
    r2exp = r2(gr, np.sign(gr[0])*np.exp(np.polyval(ce, r)))
    return freq.std(), jump, r2lin, r2exp

print(f"N={N}, K={K}, T=15000. Linear prediction: Gauss jump = d_omega/K")
print(f"{'d_omega':>8} {'jump':>9} {'d_om/K':>8} {'err%':>7} {'std(freq)':>10} "
      f"{'R2_lin':>7} {'R2_exp':>7} {'locked?':>8}")
for dw in [0.05, 0.2, 0.5, 1.0, 1.5, 2.0]:
    std, jump, r2l, r2e = relax(dw)
    err = 100 * abs(abs(jump) - dw/K) / (dw/K)
    locked = "yes" if std < 1e-2 else "NO"
    print(f"{dw:>8.2f} {jump:>9.4f} {dw/K:>8.4f} {err:>7.1f} {std:>10.2e} "
          f"{r2l:>7.4f} {r2e:>7.4f} {locked:>8}")
print()
print("Reading: while locked, jump tracks d_omega/K (Gauss holds); the run's")
print("content is WHERE locking/linearity fails, not that Poisson is 'discovered'.")
