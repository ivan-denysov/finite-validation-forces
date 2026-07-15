# -*- coding: utf-8 -*-
"""
RUN 1 — The trough with the picture's own (unscreened) delay field.
Closes the direction of debt 6.1 of "The Solitaire of Forces" in 1D.

Basis and claim
---------------
trough_profile_v4.py (paper 4, §3.1) used a phenomenological diffusive field
phi with a decay term gamma*phi, giving a screened tail e^{-r/lambda}. But the
delay field of paper 1 §4 is UNSCREENED by construction: the tick field is a
solution of the Laplace/Poisson equation ("one cannot screen against the tempo
of time"), and the composition law selected by Mercury is path-compounded
exp(phi), not the linear truncation 1-phi. This script removes the screening
(gamma = 0), solves the stationary Poisson equation for phi, and uses
V_eff = exp(-phi).

Setup
-----
Sine-Gordon kink-kink (like charges, repulsive, tail ~ 32*e^{-r}) on a 1D
lattice, plus the delay field phi obtained quasi-statically from
    d2phi/dx2 = -(s_src/D) * E(x),   phi(0) = phi(L) = 0,
where E is the local field energy density. Grounded box ends are an artefact
of the finite 1D box (image-force toward the centre ~ 2e-5, two orders below
the mutual force; stated, not hidden). Optional weak velocity damping eta is a
numerical convergence tool, NOT part of the picture; the attractor must not
depend on it (checked below).

Predictions made before running
-------------------------------
1. In 1D the unscreened potential is LINEAR in r  =>  constant attractive
   force at all separations  =>  confinement: the pair cannot be torn apart.
2. An attractor r* exists where the kink repulsion 32*e^{-r} balances the
   constant pull; naive estimate r* ~ ln(32/F_att) ~ 10 for the parameters
   below.
3. No screening length lambda exists, so the v4 failure to find a smooth
   r*(lambda) law dissolves: there is nothing for r* to be a function of.

Falsifier: no sign change of the drift (pure repulsion or pure attraction),
or an attractor that moves with the damping eta.

Results obtained with this script (for the record)
--------------------------------------------------
- r* = 10.6 +/- 0.1 from starts 6, 18, 30 at eta = 0.006;
  eta = 0.003 (under-damped) gives consistent 10.3-11.5. Attractor confirmed.
- Undamped trajectories: bound orbits through the trough; a pair kicked from
  sep = 6 out to sep ~ 85 RETURNS (linear tail = confinement).
- Source scan s_src = 1e-5 / 4e-5 / 1.6e-4 (cap not binding):
  r* = 11.6 / 10.5 / 10.3 — a plateau pinned to the scale of the short force,
  sub-logarithmic in the source over x16. This IS the s_src control named in
  debt 6.3: the plateau does not move with source strength.
- phi tail measured linear, slope = rho/2 as predicted (Gauss in 1D).

Status: confirmed in-model (1D). Open: 3D, and the bridge from the Kuramoto
network to this lattice (see run 2, field_from_dynamics_v1.py).

Usage
-----
  python3 trough_unscreened_v1.py attractor        # eta-damped scan -> r*
  python3 trough_unscreened_v1.py orbits           # undamped trajectories
  python3 trough_unscreened_v1.py source_scan      # r* vs s_src (plateau)
  python3 trough_unscreened_v1.py control          # no-phi benchmark (14->17)
"""
import numpy as np
import sys

# --- lattice and field parameters -------------------------------------------
L = 300.0; dx = 0.25; Nx = int(L / dx); c = 1.0; dt = 0.2 * dx / c
x = np.arange(Nx) * dx
D = 0.5                 # kept only to preserve the v4 source normalisation s/D
S_SRC = 4e-5            # delay-field source strength (calibrated: tent ~ 0.1)
PHI_CAP = 0.25          # safety cap; must NOT bind in a valid run (reported)

def kink(x0, s):
    return 4.0 * np.arctan(np.exp(s * (x - x0)))

def lap(f):
    l = np.empty_like(f)
    l[1:-1] = (f[2:] - 2.0 * f[1:-1] + f[:-2]) / dx**2
    l[0]  = (f[1]  - 2.0 * f[0]  + f[1])  / dx**2
    l[-1] = (f[-2] - 2.0 * f[-1] + f[-2]) / dx**2
    return l

def poisson_1d(rho):
    """Solve phi'' = -rho with phi(0)=phi(L)=0 (double cumsum + linear part)."""
    g = np.cumsum(rho) * dx
    phi = -np.cumsum(g) * dx
    phi -= phi[0] + (phi[-1] - phi[0]) * (x - x[0]) / (x[-1] - x[0])
    return phi

def kink_centers(th):
    """Sub-pixel kink positions from the two largest |grad theta| peaks."""
    g = np.abs(np.gradient(th, dx))
    idx = [i for i in range(2, Nx - 2)
           if g[i] > g[i-1] and g[i] > g[i+1] and g[i] > 0.3 * g.max()]
    idx = sorted(idx, key=lambda i: -g[i])[:2]
    out = []
    for i in sorted(idx):
        y0, y1, y2 = g[i-1], g[i], g[i+1]
        d = 0.5 * (y0 - y2) / (y0 - 2.0 * y1 + y2 + 1e-12)
        out.append(x[i] + d * dx)
    return out, len(idx)

def run(sep0, mode, eta=0.0, s_src=S_SRC, phi_cap=PHI_CAP,
        steps=30000, traj_every=0):
    """mode: 'poisson' (unscreened field) or 'none' (control, pure repulsion)."""
    th = kink(L/2 - sep0/2, +1) + kink(L/2 + sep0/2, +1)   # like charges
    tht = np.zeros_like(th)
    phi = np.zeros_like(th)
    traj = []
    for t in range(steps):
        if mode == 'poisson' and t % 10 == 0:              # quasi-static field
            gx = np.gradient(th, dx)
            E = 0.5*tht**2 + 0.5*c**2*gx**2 + (1.0 - np.cos(th))
            phi = np.clip(poisson_1d((s_src / D) * E), 0.0, phi_cap)
        Veff = np.exp(-phi) if mode == 'poisson' else 1.0  # composition law [1]
        a = c**2 * lap(th) - Veff * np.sin(th) - eta * tht
        tht += dt * a; th += dt * tht
        if traj_every and t % traj_every == 200:
            cs, n = kink_centers(th)
            traj.append(cs[1] - cs[0] if n == 2 else float('nan'))
    cs, n = kink_centers(th)
    r_fin = (cs[1] - cs[0]) if n == 2 else float('nan')
    return r_fin, phi.max(), traj

def main():
    scen = sys.argv[1] if len(sys.argv) > 1 else 'attractor'
    if scen == 'control':
        print("Control: no field — pure kink-kink repulsion (v4 benchmark 14->17)")
        for sep in [10, 14]:
            r, _, tr = run(sep, 'none', steps=24000, traj_every=24000-201)
            print(f"  sep={sep:>3}: {tr[0]:.2f} -> {r:.2f}  (repulsion expected)")
    elif scen == 'orbits':
        print("Undamped trajectories (no eta): bound orbits, return = confinement")
        for sep in [6, 12, 18, 26]:
            _, pm, tr = run(sep, 'poisson', steps=24000, traj_every=3000)
            ts = " ".join(f"{v:6.2f}" for v in tr)
            print(f"  sep={sep:>3}: {ts}   phi_max={pm:.3f}")
    elif scen == 'attractor':
        print("Damped scan (eta is a numerical tool; r* must not depend on it)")
        for eta in [0.003, 0.006]:
            for sep in [6, 18, 30]:
                r, pm, _ = run(sep, 'poisson', eta=eta)
                print(f"  eta={eta}: sep={sep:>3} -> r* = {r:6.2f}  phi_max={pm:.3f}")
    elif scen == 'source_scan':
        print("r* vs source strength (the s_src control of debt 6.3):")
        print("prediction: plateau pinned to the short-force scale")
        for s in [1e-5, 4e-5, 1.6e-4]:
            cap = 0.6 if s > 1e-4 else PHI_CAP
            r, pm, _ = run(18, 'poisson', eta=0.006, s_src=s, phi_cap=cap)
            note = "  (cap binding — INVALID)" if pm >= cap - 1e-6 else ""
            print(f"  s_src={s:.0e}: r* = {r:6.2f}  phi_max={pm:.4f}{note}")
    else:
        print(__doc__)

if __name__ == '__main__':
    main()
