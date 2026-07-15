# -*- coding: utf-8 -*-
"""
RUN 2 — The delay field built by the network itself (no field equation put in).
Completes the open item of paper 1 §4.2 ("the exact fall-off law is not
strictly established by this run") and supplies the bridge for run 1.

Setup ([1] §4.2, tightened)
---------------------------
A 1D Kuramoto chain with free ends,
    dtheta_i/dt = omega_i + K [ sin(theta_{i+1}-theta_i)
                              + sin(theta_{i-1}-theta_i) ],
mass = one node with lowered natural frequency (omega = -d_omega, background 0).

Prediction made before running
------------------------------
In the phase-locked state the stationary lag psi(x) must satisfy
    K * psi'' = Omega - omega(x),      Omega = <omega>,
i.e. the network SOLVES THE POISSON EQUATION BY ITSELF — the same unscreened
equation that run 1 imposes on the trough. Concretely:
  (a) the lag-gradient jump across the defect equals d_omega / K
      (the 1D Gauss theorem, from pure synchronisation dynamics);
  (b) the tail carries NO intrinsic screening length: any apparent decay
      scale is the causal front (how far the field has propagated),
      growing with relaxation time, not a fixed lambda.

Falsifier: an exponential fit g(r) = A e^{-r/lambda} with a
time-INDEPENDENT lambda beating the linear (Poisson) fit.

Results obtained with this script (for the record)
--------------------------------------------------
- locking: std(theta_dot) = 2.9e-5; Omega = -3.12e-5 = theory to 3 digits;
- gradient jump at the defect: 0.0497 vs theory 0.0500 (0.6%);
- tail fit r = 10..200 at T = 30000: linear R^2 = 0.9982 beats
  exponential R^2 = 0.9949;
- decisive control: the fitted "lambda" GROWS with relaxation time,
  lambda(T=7500) = 86 -> lambda(T=30000) = 224 (ratio 2.6; a screening
  length would give 1.0). The apparent scale is the causal front —
  retardation, not screening.

Status: confirmed in-model (1D). The field the network builds around mass is
the unscreened Poisson field; feeding it to run 1 replaces the v4 proxy phi
with the picture's own delay field — the 1D bridge of debt 6.1.

Usage
-----
  python3 field_from_dynamics_v1.py            # main run, T = 30000
  python3 field_from_dynamics_v1.py front      # lambda(T) front control
"""
import numpy as np
import sys

N = 1601; K = 1.0; dt = 0.05
IC = N // 2                    # defect (mass) at the centre
D_OMEGA = 0.05                 # frequency deficit of the heavy node

def relax(T):
    w = np.zeros(N); w[IC] = -D_OMEGA
    th = np.zeros(N)
    for _ in range(int(T / dt)):
        s = np.sin(np.diff(th))
        f = np.zeros(N); f[:-1] += K * s; f[1:] -= K * s
        th += dt * (w + f)
    return th, w

def r2(y, fit):
    return 1.0 - np.sum((y - fit)**2) / np.sum((y - np.mean(y))**2)

def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'front':
        print("Front control: fitted 'lambda' must grow with T if unscreened")
        lams = []
        for T in [7500, 30000]:
            th, _ = relax(T)
            g = np.gradient(th - th[IC])
            r = np.arange(10, 201); gr = g[IC+10:IC+201]
            ce = np.polyfit(r, np.log(np.abs(gr)), 1)
            lams.append(-1.0 / ce[0])
            print(f"  T={T:>6}: lambda_fit = {lams[-1]:.0f}")
        print(f"  ratio = {lams[1]/lams[0]:.2f}  "
              f"(front ~ sqrt(4) = 2.0; true screening = 1.0)")
        return

    T = 30000
    th, w = relax(T)
    # locking check
    s = np.sin(np.diff(th))
    f = np.zeros(N); f[:-1] += K * s; f[1:] -= K * s
    freq = w + f
    print(f"locking: std(theta_dot) = {freq.std():.2e}; "
          f"Omega = {freq.mean():+.2e} (theory {-D_OMEGA/N:+.2e})")
    # Gauss jump at the defect
    g = np.gradient(th - th[IC])
    print(f"gradient jump at defect: {g[IC+2]-g[IC-2]:+.4f} "
          f"(theory {D_OMEGA/K:+.4f})")
    # tail: linear (Poisson) vs exponential (screened)
    r = np.arange(10, 201); gr = g[IC+10:IC+201]
    cl = np.polyfit(r, gr, 1)
    r2_lin = r2(gr, np.polyval(cl, r))
    ce = np.polyfit(r, np.log(np.abs(gr)), 1)
    r2_exp = r2(gr, np.sign(gr[0]) * np.exp(np.polyval(ce, r)))
    print(f"tail fit r=10..200: LINEAR R^2 = {r2_lin:.4f}  |  "
          f"EXP R^2 = {r2_exp:.4f} (lambda_fit = {-1/ce[0]:.0f})")
    print("profile (r, gradient, lag):")
    psi = th - th[IC]
    for rr in [2, 5, 10, 20, 50, 100, 150, 200, 300]:
        print(f"  r={rr:>4}: g = {g[IC+rr]:+.3e}   psi = {psi[IC+rr]:+.4f}")

if __name__ == '__main__':
    main()
