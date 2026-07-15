# -*- coding: utf-8 -*-
"""
FIGURE for paper 4 v2, Section 5 (Gravity as the Proto-Force).
Four panels, one per run. Panel B recomputes the delay-field profile live
(cheap, reduced N for the figure and noted); panels A, C, D use the measured
values from the protocol scripts (baked in, provenance in each panel title).
Panel C now carries the multi-seed spread (B.4), so the figure is graded to
the same honesty as the v2 text. Panel titles carry the mark GRADE.
Output: fig_protoforce.png (PNG, thin lines on light bg).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE   = "#1f5fbf"   # gravity / field
RED    = "#c0392b"   # strong / trough
GOLD   = "#b8860b"   # this-work overlay
GREEN  = "#2e7d32"   # weak
PURPLE = "#6a3fb5"   # EM / metric
GREY   = "#888888"

fig, ax = plt.subplots(2, 2, figsize=(11, 8.4))
fig.suptitle("Gravity as the Proto-Force — four in-model marks (1D), of unequal strength",
             fontsize=13.5, fontweight="bold")

# ---------------------------------------------------------------- Panel A
# Run 1: the trough is an attractor. Final vs starting separation.
a = ax[0, 0]
starts = np.array([6, 9, 12, 15, 18, 22, 26, 30])
r_unscr = np.array([10.6]*len(starts))
r_scr   = np.array([7.17]*len(starts))
a.plot([5, 31], [5, 31], "--", color=GREY, lw=1, label="no motion (r_fin=r_0)")
a.scatter(starts, r_unscr, color=RED, s=45, zorder=3,
          label="unscreened field  r*=10.6 (this work)")
a.scatter(starts, r_scr, facecolors="none", edgecolors=GOLD, s=45, zorder=3,
          label="screened proxy  r*=7.17 (v1)")
a.axhline(10.6, color=RED, lw=0.8, alpha=0.5)
a.axhline(7.17, color=GOLD, lw=0.8, alpha=0.5)
for s in starts:
    a.annotate("", xy=(s, 10.6), xytext=(s, s),
               arrowprops=dict(arrowstyle="->", color=RED, alpha=0.35, lw=0.8))
a.set_xlabel("starting separation  r_0")
a.set_ylabel("final separation  r_fin")
a.set_title("A. The strong: a trough with an attractor\n"
            "(compositional mark — the substantive one)", fontsize=10)
a.legend(fontsize=7.5, loc="upper left")
a.set_xlim(5, 31); a.set_ylim(5, 31)

# ---------------------------------------------------------------- Panel B
# Run 2: consistency check — the linearised network satisfies Poisson.
b = ax[0, 1]
N = 801; K = 1.0; dt = 0.05; IC = N // 2; DW = 0.05     # reduced N for the figure
def relax(T):
    w = np.zeros(N); w[IC] = -DW; th = np.zeros(N)
    for _ in range(int(T / dt)):
        s = np.sin(np.diff(th))
        f = np.zeros(N); f[:-1] += K*s; f[1:] -= K*s
        th += dt*(w + f)
    return th
th = relax(15000)
g = np.gradient(th - th[IC])
r_axis = np.arange(-120, 121)
b.plot(r_axis, g[IC-120:IC+121], color=BLUE, lw=1.6)
b.axvline(0, color=GREY, lw=0.7, ls=":")
jump = g[IC+2] - g[IC-2]
b.annotate(f"Gauss jump = {jump:.4f}\n(theory {DW/K:.4f})",
           xy=(0, 0), xytext=(28, g.max()*0.55), fontsize=8,
           arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
b.set_xlabel("distance from mass  r")
b.set_ylabel("lag gradient  d(psi)/dx")
b.set_title("B. Gravity: the delay field from the network\n"
            "(consistency check, weak field)", fontsize=10)
ins = b.inset_axes([0.60, 0.14, 0.36, 0.30])
Ts = [3000, 7500, 15000]
lams = []
for T in Ts:
    thh = relax(T); gg = np.gradient(thh - thh[IC])
    rr = np.arange(10, 121); grr = gg[IC+10:IC+121]
    ce = np.polyfit(rr, np.log(np.abs(grr)), 1); lams.append(-1/ce[0])
ins.plot(Ts, lams, "o-", color=GOLD, ms=4, lw=1)
ins.set_title("lambda grows with T (front)", fontsize=6.5)
ins.tick_params(labelsize=6)
ins.set_xlabel("relax time T", fontsize=6.5)
ins.set_ylabel("fitted lambda", fontsize=6.5)

# ---------------------------------------------------------------- Panel C
# Run 3: weak lifetime vs template delay — MULTI-SEED (B.4), with spread.
# median-of-medians and [min,max] over seeds 900/901/902 (tau_multiseed.py).
c = ax[1, 0]
dw  = np.array([0.85, 1.00, 1.15, 1.30])
med = np.array([139.0, 158.0, 24.0, 8.5])
lo  = np.array([112.0, 87.0, 20.0, 3.8])
hi  = np.array([233.0, 173.0, 91.0, 50.0])
yerr = np.vstack([med - lo, hi - med])
c.axvspan(0.60, 0.78, color=GREEN, alpha=0.10)
c.text(0.615, 300, "stable\n(below barrier)", fontsize=7.5, color=GREEN)
c.scatter([0.70], [430], marker="v", color=GREEN, s=55)          # stable point
c.errorbar(dw, med, yerr=yerr, fmt="o-", color=GREEN, ms=6, lw=1.4,
           capsize=4, ecolor=GREEN, elinewidth=1,
           label="median-of-medians ± seed spread")
c.text(1.00, 205, "reference", fontsize=7, ha="center")
c.text(1.28, 40, "×2–5 seed\nscatter; slope\nnot a clean law",
       fontsize=6.5, ha="right", color=GREY)
c.set_yscale("log")
c.set_xlabel("template delay  d_omega  (= template mass)")
c.set_ylabel("median lifetime  tau  (log)")
c.set_title("C. The weak: clock set by template delay\n"
            "(structural mark — partly definitional; multi-seed)", fontsize=10)
c.set_xlim(0.60, 1.38)
c.legend(fontsize=6.8, loc="lower left")

# ---------------------------------------------------------------- Panel D
# Run 4: in-model Shapiro (eikonal check). measured vs predicted.
d = ax[1, 1]
pred = np.array([44.59, 87.69])
meas = np.array([44.61, 87.79])
lim = [0, 100]
d.plot(lim, lim, "--", color=GREY, lw=1, label="identity (measured = predicted)")
d.scatter(pred, meas, color=PURPLE, s=70, zorder=3)
for p, m in zip(pred, meas):
    d.annotate(f"ratio {m/p:.3f}", xy=(p, m), xytext=(p+4, m-6), fontsize=8)
d.text(20, 80, "lag doubles when\nmass doubles", fontsize=8, color=PURPLE)
d.set_xlabel("predicted lag  = eikonal integral of (e^{2phi}-1)")
d.set_ylabel("measured transit lag")
d.set_title("D. Electromagnetism: the metric mark\n"
            "(consistency check, eikonal)", fontsize=10)
d.legend(fontsize=7.5, loc="lower right")
d.set_xlim(0, 100); d.set_ylim(0, 100)

for row in ax:
    for a_ in row:
        a_.grid(True, alpha=0.25, lw=0.5)

fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig("fig_protoforce.png", dpi=170)
print("saved fig_protoforce.png")
print(f"panel B (reduced N={N}, T=15000): Gauss jump {jump:.4f} vs {DW/K:.4f}; "
      f"lambda(T) = {[round(l) for l in lams]}")
