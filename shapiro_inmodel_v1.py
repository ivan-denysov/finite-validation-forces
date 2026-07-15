# -*- coding: utf-8 -*-
"""
RUN 4 — The metric mark: a massless carrier is slowed by the delay field
(in-model Shapiro delay).

Basis and claim
---------------
In the picture the carrier of electromagnetism (the photon: confirmation
transfer along edges, [2] §4) couples to no charge of the delay field — yet
it cannot help propagating at the tempo the delay field sets, because the
tempo of ticks IS the medium ([1] §4: gravity is a state of the network).
This is the cleanest "mark of gravitational origin" for EM, where no
compositional mark exists: every carrier lives on the metric that gravity is.

Setup
-----
A massless wave packet on a 1D lattice with local propagation speed set by
the delay field phi (same Poisson machinery and normalisation as run 1,
NO retuning):
    c_loc(x) = c * exp(-2*phi(x)),
the factor 2 imported from [1] §4: two equal channels (time + space), their
equality derived there from the invariance of c = L/tau. phi is the tent
built by a point mass m at the centre (Poisson, grounded box ends).

Predictions made before running
-------------------------------
1. Transit lag vs the empty run: dt = integral (n - 1) dx with n = e^{2 phi}
   along the traversed path — the Shapiro form (the lag is set by the field
   integral along the path, not by local gradients).
2. Linearity in the mass: m -> 2m gives dt -> 2*dt (weak field).

Falsifier: a lag that does not match the path integral, or scales
non-linearly with m in the weak field.

Note on register: the factor 2 (two channels) is an input from [1], not
re-derived here; what the run tests is that ONE field, with the SAME
normalisation that binds kinks in run 1, also retards a chargeless carrier
by exactly its path integral. One deficit, two faces.

Usage:  python3 shapiro_inmodel_v1.py
"""
import numpy as np

L = 1000.0; dx = 0.25; Nx = int(L / dx); c = 1.0; dt = 0.2 * dx / c
x = np.arange(Nx) * dx

def poisson_1d(rho):
    g = np.cumsum(rho) * dx
    phi = -np.cumsum(g) * dx
    phi -= phi[0] + (phi[-1] - phi[0]) * (x - x[0]) / (x[-1] - x[0])
    return phi

def phi_of_mass(m):
    """Delay field of a point mass at the centre (same machinery as run 1)."""
    rho = np.exp(-((x - L/2)/4.0)**2)
    rho *= m / (rho.sum() * dx)                 # integral rho dx = m
    return poisson_1d(rho)

def transit(phi, T=700.0, x0=100.0, w=8.0):
    """Right-moving packet; returns energy-centroid position at time T."""
    n2 = np.exp(-4.0 * phi)                     # c_loc^2 = e^{-4 phi}
    th = np.exp(-((x - x0)/w)**2)
    tht = -c * np.gradient(th, dx)              # d'Alembert right-mover
    for _ in range(int(T / dt)):
        lapl = np.empty_like(th)
        lapl[1:-1] = (th[2:] - 2*th[1:-1] + th[:-2]) / dx**2
        lapl[0] = lapl[-1] = 0.0
        tht += dt * n2 * lapl
        th += dt * tht
    E = tht**2 + np.gradient(th, dx)**2
    return (x * E).sum() / E.sum()

def predicted_lag(phi, x_from, x_to):
    sel = (x >= x_from) & (x <= x_to)
    return np.trapezoid(np.exp(2.0 * phi[sel]) - 1.0, dx=dx)

print("empty run (control):")
xc0 = transit(np.zeros(Nx))
print(f"  centroid at T=700: x = {xc0:.2f}  (ballistic: 800.00)")

for m in [2e-4, 4e-4]:
    phi = phi_of_mass(m)
    xc = transit(phi)
    lag_meas = (xc0 - xc) / c
    lag_pred = predicted_lag(phi, 100.0, xc)
    print(f"m = {m:.0e}: phi_max = {phi.max():.4f}, centroid x = {xc:.2f}")
    print(f"  measured lag dt = {lag_meas:6.2f}   predicted path integral = "
          f"{lag_pred:6.2f}   ratio = {lag_meas/lag_pred:.3f}")
