# -*- coding: utf-8 -*-
"""
5.2 MULTI-SEED (reviewer request): the v1 table is one realisation (seed 900).
Re-run tau(d_omega) across several network+state seeds; report the per-seed
median lifetime and the spread, so the "x30 for +15% mass" claim carries an
error bar and the non-monotonic decayed-count is explained as seed noise.

Reduced cost (N, T shortened) so absolute tau differ from Table 5.2; the TEST
is the monotonicity of the median and its spread, not the absolute values.
Barrier is frozen per seed from that seed's reference (d_omega=1.0) q99.

Usage: python3 tau_multiseed.py <d_omega> <seed1> <seed2> ...
"""
import numpy as np, sys

N = 1500; DEG = 6; K = 1.0; DT = 0.02; M = 60
T_RELAX = 10000; T_BG = 4000; T_DECAY = 40000
Q = 0.99; DW_REF = 1.0

def build(n, deg, seed, box=30.0):
    rng = np.random.RandomState(seed)
    pos = rng.uniform(0, box, (n, 3)); cell = box/8; grid = {}
    for i, p in enumerate(pos):
        grid.setdefault(tuple((p//cell).astype(int)), []).append(i)
    nbrs = []
    for i, p in enumerate(pos):
        key = (p//cell).astype(int); cand = []
        for a in (-1,0,1):
            for b in (-1,0,1):
                for c in (-1,0,1):
                    cand += grid.get((key[0]+a, key[1]+b, key[2]+c), [])
        cand = [j for j in cand if j != i]
        d = np.sqrt(((pos[cand]-p)**2).sum(1))
        nbrs.append([cand[k] for k in np.argsort(d)[:deg]])
    maxd = max(len(nn) for nn in nbrs); na = np.full((n, maxd), -1)
    for i, nn in enumerate(nbrs): na[i, :len(nn)] = nn
    valid = na >= 0; cnt = valid.sum(1); cnt[cnt == 0] = 1
    return na, valid, cnt

def make_step(na, valid, cnt):
    def step(th, om):
        tn = np.where(valid, th[na], 0.0)
        coup = np.where(valid, np.sin(tn - th[:, None]), 0.0).sum(1)/cnt
        return th + DT*(om + K*coup)
    def mis(th, h):
        return np.abs(np.where(valid[h], np.sin(th[na[h]]-th[h]), 0.0)).sum()/cnt[h]
    return step, mis

def relax_bg(dw, seed, step, mis):
    rng = np.random.RandomState(seed)
    th = rng.uniform(0, 2*np.pi, N)
    heavy = rng.choice(N, M, replace=False).tolist()
    om = np.ones(N); om[heavy] = 1.0 + dw
    for _ in range(T_RELAX): th = step(th, om)
    samp = []; thb = th.copy()
    for t in range(T_BG):
        thb = step(thb, om)
        if t % 20 == 0: samp += [mis(thb, h) for h in heavy]
    return th, heavy, np.array(samp)

def lifetimes(dw, seed, dE, step, mis):
    th, heavy0, _ = relax_bg(dw, seed, step, mis)
    om = np.ones(N); om[heavy0] = 1.0 + dw
    heavy = set(heavy0); times = []
    for t in range(T_DECAY):
        th = step(th, om)
        if t % 25 == 0 and heavy:
            for h in list(heavy):
                if mis(th, h) > dE:
                    om[h] = 1.0; heavy.discard(h); times.append(t*DT)
    return times

dw = float(sys.argv[1]); seeds = [int(s) for s in sys.argv[2:]]
print(f"d_omega={dw}: per-seed (network+state seed shared)")
meds = []
for seed in seeds:
    na, valid, cnt = build(N, DEG, seed)
    step, mis = make_step(na, valid, cnt)
    _, _, bg = relax_bg(DW_REF, seed, step, mis)   # frozen barrier per seed
    dE = np.quantile(bg, Q)
    t = lifetimes(dw, seed, dE, step, mis)
    med = np.median(t) if t else float('inf')
    meds.append(med)
    print(f"  seed={seed}: decayed {len(t):>2}/{M}, dE={dE:.4f}, median tau={med:8.1f}")
finite = [m for m in meds if np.isfinite(m)]
if finite:
    print(f"  => median-of-medians={np.median(finite):.1f}, "
          f"spread [{min(finite):.1f}, {max(finite):.1f}]")
