# finite-validation-forces

Simulation code and numerical protocols for:

**I. Denysov, *Physics from Finite Validation: The Solitaire of Forces — Reading the Four Interactions Through Finite Validation*, United Field Initiative, 2026.**
DOI: [to be assigned after Zenodo deposit]

Series: [Paper 1](https://doi.org/10.5281/zenodo.20633926) · [Paper 2](https://doi.org/10.5281/zenodo.20634290) · [Paper 3](https://doi.org/10.5281/zenodo.20640511) · this work.

Requirements: Python 3.10+, numpy. Figure scripts additionally: matplotlib, Pillow. Every script runs standalone; seeds are in the code. `fig_overlay_forces_v2.py` overlays Tier 0 onto `figure1_four_ring_continuum.png` (the original master figure of Paper 0 Part 3, included here for reproducibility) and produces `figure1_with_forces_tier0.png`.

---

# Numerical Companion: Interference Troughs, Template Relaxation, and an Entropy Test in Finite-Validation Models

*Ivan Denysov, Lead Researcher, United Field Initiative (UFI). Supplementary Material to: "Physics from Finite Validation: The Solitaire of Forces" (this record). Code: https://github.com/ivan-denysov/finite-validation-forces*

## Purpose and Status Discipline

This companion contains every numerical run supporting the main article, with full code, parameters, seeds, and per-result status labels. The discipline follows the series: **confirmed in-model**, **observation**, **structural result**, **rejected** (tested and failed), **lesson/artefact** (a result traced to the setup, documented). Runs are exploratory instruments within stated models, not measurements of nature. The programme-wide instrument note holds here too: forms of laws reproduce at R² > 0.98 while coefficients run 0.8–0.85 of analytic predictions — the form is the claim, the coefficient is instrument-limited.

The work uses two model regimes inherited from [DOI 10.5281/zenodo.20634290]: a Kuramoto-type phase network (relaxation, entropy) and a sine-Gordon field (the trough). Their continuum bridge is the programme's load-bearing open debt. All scripts are standalone Python 3 (numpy; figures additionally use matplotlib + Pillow).

---

## Block A — The Strong as an Interference Trough (main article §3.1)

**A.1 Setup lessons [lessons, documented].**
- `trough_profile.py` — first attempt (kink–antikink, instantaneous V_eff(E)): the observation window (600 steps) was ~13× shorter than the sensitivity threshold; centres read off a coarse gradient threshold were noise. Lesson preserved.
- `trough_profile_v2.py` — honest window (24 000 steps), sub-pixel centres: the control is clean (14→9.70, the kink–antikink benchmark), but an instantaneous V_eff(E) coupling is an **active medium** — it injects energy (E/E₀ drifts +18–24%) and kills the long tail. Lesson: gravity must be a *diffusion* of the dip, not an instantaneous response.
- `trough_profile_v3.py` — soft probe (g=0.12, clip ≥0.6): same conclusion — the active coupling damps the tails (everything neutral at sep≥13). The two lessons together motivated the reversal of the setup.

**A.2 The trough as an attractor [confirmed in-model].**
- `trough_profile_v4.py` — the reversal: **kink–kink** (like charges, SG repulsion, e^(−r) tail) **+ a phenomenological diffusive field** φ (a proxy for a long attractive profile, *not* the gravitational field of the picture) (∂φ/∂t = D∇²φ + s·E − γφ, V_eff = 1−φ, λ=√(D/γ)=10; CFL fixed at D=0.5). From all starting separations 9, 12, 15, 18, 22, 26, 30 the pair converges to r* ≈ 7.17±0.01; from a start at 6 it is pushed out to the same r* — an attractor from both sides. The control without φ is clean repulsion (14→17). Structurally the configuration of a nucleon: like-charged constituents bound at a finite scale, repelling under compression. CFL lesson preserved (D=2.0 overflowed the diffusion step).

**A.3 The λ-scan: two regimes, no formula [confirmed in-model + honest failure].**
- `trough_lambda_scan.py` — r*(λ), convergence from both sides (agreement 0.00–0.01 for λ≤14): a plateau r* ≈ 7.1–7.3 for λ=5/7/10 (the bottom pinned to the scale of the short force — robust, like nucleon size), then a step to 13.4 (λ=14) and ~26 (λ=20, agreement degrading to 0.72, the well shallowing). No smooth r*(λ): both a log fit and a power fit fail at R²≈0.58 — reported as a failed fit, not hidden. Reading: two regimes — a hard trough (bottom at the short-force scale, insensitive to the far tail) and a soft well at large λ. Open: the s_src control (the plateau should not move with source strength), the confinement growth law, 1D vs 3D.

---

## Block B — The Weak as Template Relaxation (main article §3.2)

**B.1 Setup lessons [lessons, documented].**
- `template_relaxation.py` — thresholds below the un-warmed phase noise: every carrier decays instantly. Lesson: warm the phase first.
- `template_relaxation_v2.py` — the desynchrony metric is |sin|≤1, the warmed background sits at ~0.82 near the ceiling, so μ+kσ thresholds are unreachable (0 decays). Lesson: a decay threshold must live in the *tail* of the distribution — decay is a rare (Poisson) event.

**B.2 The lifetime hierarchy [confirmed in-model] + two observations.**
- `template_relaxation_v3.py` — thresholds at background quantiles. Shifting the barrier in the third decimal place (q95→q99→q999: 0.9868 → 0.9891 → 0.9912) multiplies the median lifetime fourfold (117 → 213 → 463). The hypersensitivity of τ to the barrier is the profile of the weak (natural lifetimes span ~43 orders of magnitude). **Observation:** at a common barrier the ensemble splits into decayed and stable (27–29 of 60 decay, then stop) — environmental inhomogeneity = inhomogeneity of fate (stable vs unstable isotopes); the survival curve is a mixture of exponentials, R²=0.88.
- `template_relaxation_v4.py` — individual (own-background quantile) thresholds: R²=0.78, 14 survivors. **Observation:** decays are correlated through the medium — each decay clears its neighbours' background, so survivors stabilise (cooperative stabilisation); a lone heavy carrier in a cleared network lives indefinitely. **Open:** a clean single exponential is not reproduced (R²=0.78–0.88 in either normalisation); candidate causes — the saturating metric (|sin| ceiling) vs genuine medium correlation; the discriminator (a non-saturating metric) is a next-session run.

---

## Block C — Relaxation and Entropy (main article §3.3)

**C.1 The H-test [confirmed in-model].**
- `entropy_h_test.py` — first attempt on a non-symmetric kNN graph: the asymmetric coupling makes the forces non-conservative (a hidden pump; E/E₀ ran to ×104). The test does not count — preserved as a lesson (non-symmetric graph = hidden pumping).
- `entropy_h_test_v2.py` — a symmetrised graph (mutual links → pairwise, conservative forces), leapfrog integration. A closed network (no damping, no pumping) started with all disturbance localised in 10 of 2000 nodes: the Shannon entropy of the energy distribution rises monotonically (11/15 windows within fluctuation) from S=2.31 to a plateau ~7.05 (against ln N = 7.60 — 93% of equipartition), energy conserved to E/E₀=1.0001 over 120 000 steps, no return to localisation. A localised disturbance irreversibly spreads to equipartition — H-behaviour from the network's own dynamics, no thermodynamic postulate. Open: the link of Shannon S to thermodynamic S (additivity, temperature, the full second law) — a separate work.

---

## Block D — The Master Figure (main article §4)

- `fig_overlay_forces.py`, `fig_overlay_forces_v2.py` — the Tier-0 overlay on the original [Part 3, Fig. 1]. The original raster is preserved pixel-for-pixel; the right panel is replaced by a "Forces open the rings" column (numbered 1 gravity → 2 strong → 3 EM → 4 weak, bottom-to-top, mirroring the ring anchors), and a Tier-0 row is added at blueprint rank 1. The figure is **schematic**: the entire force sequence sits at ≈13.8 Gyr before present on this axis; the quantitative curves are those of [Part 3, Fig. 1]. PNG (not JPEG — thin lines on a light background).

---

## Reproducibility

Python 3.10+, numpy; figures additionally matplotlib and Pillow. Every script runs standalone. Seeds are written in the code; per-result outputs are printed. Runtime: seconds to ~minutes per script; the longest are the trough scans and the relaxation runs (a few minutes each on a laptop).

## File → article map

| Script | Section | Status |
|---|---|---|
| trough_profile.py | §3.1 | lesson (window too short) |
| trough_profile_v2.py | §3.1 | lesson (active-medium coupling) |
| trough_profile_v3.py | §3.1 | lesson (active coupling damps tails) |
| trough_profile_v4.py | §3.1 | confirmed in-model (a two-profile trough attractor r*≈7.17; φ a phenomenological proxy, not gravity) |
| trough_lambda_scan.py | §3.1 | confirmed in-model (two regimes) + failed fit (no r*(λ) law) |
| template_relaxation.py | §3.2 | lesson (threshold below noise) |
| template_relaxation_v2.py | §3.2 | lesson (saturated metric; threshold must be in the tail) |
| template_relaxation_v3.py | §3.2 | confirmed in-model (lifetime hierarchy) + observation (isotope split) |
| template_relaxation_v4.py | §3.2 | observation (cooperative stabilisation) + open (no clean exponential) |
| entropy_h_test.py | §3.3 | lesson (non-symmetric graph = hidden pump) |
| entropy_h_test_v2.py | §3.3 | confirmed in-model (H-behaviour) |
| fig_overlay_forces.py | §4 | figure (v1) |
| fig_overlay_forces_v2.py | §4 | figure (final overlay) |

## License
Code MIT; text and figures CC BY 4.0. © United Field Initiative. The master figure overlays the original Figure 1 of [Part 3] (same author/programme).

## Research Programme
This document is one output of a broader, staged research programme of the United Field Initiative on cumulative complexity and validation dynamics. Details: ufi.observer/research-programme.
