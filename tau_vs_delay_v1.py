# -*- coding: utf-8 -*-
"""
RUN 3 — Lifetime vs template delay: the gravitational mark on the weak channel.

Basis and claim
---------------
template_relaxation_v3.py (paper 4, §3.2) confirmed the lifetime hierarchy in
the barrier: a third-decimal shift of the transition cost multiplies the
median lifetime severalfold. This run scans the OTHER axis: the barrier is
held fixed (one transition cost, set by the environment) and the TEMPLATE
DELAY is varied — the tick-raise d_omega of the heavy nodes, which in the
picture IS the template's mass ([1] §4: mass = tick delay; [3] §2).

If the weak channel carries the mark of its gravitational origin, the
lifetime must be set by the template's delay measured against a fixed
environmental cost: heavier template -> larger local desynchronisation ->
the tail of the background distribution is reached sooner -> shorter life.
(Structurally: the rough mass-lifetime trend of real particle physics.)

Protocol
--------
3D random geometric network (N nodes, degree ~6), Kuramoto dynamics, M heavy
nodes with omega = 1 + d_omega. The barrier dE is fixed ONCE, from the
background mismatch statistics of the reference ensemble d_omega = 1.0
(the v3 q99 protocol), and reused unchanged for every d_omega.
A heavy node decays (omega -> 1) when its local mismatch exceeds dE.
Median lifetime and decayed fraction are reported per d_omega.

Prediction made before running: median tau falls monotonically — and steeply
(the hierarchy axis) — with d_omega at fixed dE; light templates
(d_omega well below the reference) sit below the barrier and are stable
(the in-model analogue of stable species).

Falsifier: tau flat or rising with d_omega.

Status criteria: monotone hierarchy = confirmed in-model; the clean single
exponential of the survival curve is NOT claimed here (debt 6.3, saturating
|sin| metric named suspect in paper 4).

Usage:  python3 tau_vs_delay_v1.py <d_omega> [<d_omega> ...]
        (one value per call keeps runs short; barrier is recomputed
         deterministically from the reference ensemble each call, seed-fixed)
"""
import numpy as np
import sys

N = 3000; DEG = 6; K = 1.0; DT = 0.02; M = 60
T_RELAX = 20000; T_BG = 6000; T_DECAY = 80000
Q = 0.99            # reference quantile defining the fixed environmental cost
DW_REF = 1.0        # reference tick-raise (v3: omega_heavy = 2.0)

def build(n, deg, seed, box=30.0):
    rng = np.random.RandomState(seed)
    pos = rng.uniform(0, box, (n, 3))
    cell = box / 8; grid = {}
    for i, p in enumerate(pos):
        grid.setdefault(tuple((p // cell).astype(int)), []).append(i)
    nbrs = []
    for i, p in enumerate(pos):
        key = (p // cell).astype(int); cand = []
        for a in (-1, 0, 1):
            for b in (-1, 0, 1):
                for c in (-1, 0, 1):
                    cand += grid.get((key[0]+a, key[1]+b, key[2]+c), [])
        cand = [j for j in cand if j != i]
        d = np.sqrt(((pos[cand] - p)**2).sum(1))
        nbrs.append([cand[k] for k in np.argsort(d)[:deg]])
    maxd = max(len(nn) for nn in nbrs)
    na = np.full((n, maxd), -1)
    for i, nn in enumerate(nbrs):
        na[i, :len(nn)] = nn
    valid = na >= 0; cnt = valid.sum(1); cnt[cnt == 0] = 1
    return na, valid, cnt

NA, VALID, CNT = build(N, DEG, 0)

def step(th, om):
    tn = np.where(VALID, th[NA], 0.0)
    coup = np.where(VALID, np.sin(tn - th[:, None]), 0.0).sum(1) / CNT
    return th + DT * (om + K * coup)

def mismatch(th, h):
    return np.abs(np.where(VALID[h], np.sin(th[NA[h]] - th[h]), 0.0)).sum() / CNT[h]

def relax_and_background(dw):
    rng = np.random.RandomState(900)
    th = rng.uniform(0, 2*np.pi, N)
    heavy = rng.choice(N, M, replace=False).tolist()
    om = np.ones(N); om[heavy] = 1.0 + dw
    for _ in range(T_RELAX):
        th = step(th, om)
    samples = []
    thb = th.copy()
    for t in range(T_BG):
        thb = step(thb, om)
        if t % 20 == 0:
            samples += [mismatch(thb, h) for h in heavy]
    return th, heavy, np.array(samples)

def lifetime_scan(dw, dE):
    th, heavy0, _ = relax_and_background(dw)
    om = np.ones(N); om[heavy0] = 1.0 + dw
    heavy = set(heavy0); times = []
    for t in range(T_DECAY):
        th = step(th, om)
        if t % 25 == 0 and heavy:
            for h in list(heavy):
                if mismatch(th, h) > dE:
                    om[h] = 1.0; heavy.discard(h); times.append(t * DT)
    med = np.median(times) if times else float('inf')
    return len(times), med

def main():
    dws = [float(v) for v in sys.argv[1:]] or [DW_REF]
    # fixed environmental cost from the reference ensemble (computed once,
    # deterministic: same seeds every call)
    _, _, bg = relax_and_background(DW_REF)
    dE = np.quantile(bg, Q)
    print(f"fixed barrier dE = {dE:.4f}  (q{int(Q*100)} of reference "
          f"d_omega = {DW_REF} background; frozen for all templates)")
    for dw in dws:
        dec, med = lifetime_scan(dw, dE)
        print(f"d_omega = {dw:4.2f}: decayed {dec:>2}/{M}, "
              f"median tau = {med:8.1f}")

if __name__ == '__main__':
    main()
