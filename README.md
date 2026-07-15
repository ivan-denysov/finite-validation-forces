# finite-validation-forces

Simulation code and numerical protocols for:

**I. Denysov, *Physics from Finite Validation: The Solitaire of Forces — Reading the Four Interactions Through Finite Validation, and Gravity as the Proto-Force*, United Field Initiative, 2026 (v2).**
DOI: [10.5281/zenodo.20657843](https://doi.org/10.5281/zenodo.20657843)

Series: [Paper 1](https://doi.org/10.5281/zenodo.21266772) · [Paper 2](https://doi.org/10.5281/zenodo.20634290) · [Paper 3](https://doi.org/10.5281/zenodo.20640511) · this work · [The Bridge](https://doi.org/10.5281/zenodo.21275364).

Requirements: Python 3.10+, numpy. Figure scripts additionally: matplotlib, Pillow. Every script runs standalone; seeds are in the code.

## What's new in v2

v2 pays the load-bearing debt of v1 (the strong's long attractive profile was a screened proxy) with the picture's **own unscreened delay field**, and adds a **proto-force** section (§5) with a graded set of marks — compositional (strong, substantive), structural (weak, partly definitional), metric (EM, near-tautological). Superseded v1 results are re-read, not removed.

Seven new scripts:

| Script | What it does | Status |
|---|---|---|
| `trough_unscreened_v1.py` | trough with the unscreened Poisson field (γ=0, V_eff=e^−φ); attractor r*=10.6, confinement tail, s_src control | confirmed in-model (1D) |
| `field_from_dynamics_v1.py` | Kuramoto chain builds the delay field itself; 1D Gauss jump 0.6% | consistency check (weak field) |
| `e1_stress_domega.py` | where the linearised Gauss law breaks as sin saturates (Δω≳1) | stress test |
| `tau_vs_delay_v1.py` | weak lifetime vs template delay at frozen barrier (single seed) | confirmed in-model (direction) |
| `tau_multiseed.py` | the same across seeds 900/901/902: monotone decayed fraction, large slope scatter | observation (with scatter) |
| `shapiro_inmodel_v1.py` | in-model Shapiro: chargeless carrier retarded by the eikonal path integral | consistency check (eikonal) |
| `fig_protoforce_v1.py` | the four-panel §5 figure → `fig_protoforce.png` (Figure 2) | figure |

**Prediction discipline.** Each new run carries its prediction in the header, committed ahead of the result, so "predictions written before the run" is an auditable git history. See `PREDICTIONS.md`.

## File → article map (v2)

| Script | Section | Status |
|---|---|---|
| trough_profile.py | §3.1 | lesson (window too short) |
| trough_profile_v2.py | §3.1 | lesson (active-medium coupling) |
| trough_profile_v3.py | §3.1 | lesson (active coupling damps tails) |
| trough_profile_v4.py | §3.1 | confirmed in-model (proxy trough r*≈7.17; φ a proxy, stage 1) |
| trough_lambda_scan.py | §3.1 | confirmed in-model (two regimes) + failed fit; **[v2]** re-read as a proxy artefact |
| **trough_unscreened_v1.py** | §3.1, §5 | **[v2]** confirmed in-model (unscreened field, r*=10.6, confinement, s_src control) |
| template_relaxation.py | §3.2 | lesson (threshold below noise) |
| template_relaxation_v2.py | §3.2 | lesson (saturated metric; threshold in the tail) |
| template_relaxation_v3.py | §3.2 | confirmed in-model (barrier hierarchy) + observation (isotope split) |
| template_relaxation_v4.py | §3.2 | observation (cooperative stabilisation) + open (no clean exponential) |
| **tau_vs_delay_v1.py** | §3.2, §5 | **[v2]** confirmed in-model (delay-axis hierarchy, single seed) |
| **tau_multiseed.py** | §5.2 | **[v2]** observation (multi-seed: monotone decayed fraction; slope has large seed scatter) |
| entropy_h_test.py | §3.3 | lesson (non-symmetric graph = hidden pump) |
| entropy_h_test_v2.py | §3.3 | confirmed in-model (H-behaviour) |
| **field_from_dynamics_v1.py** | §5.1 | **[v2]** consistency check (linearised chain = Poisson; closes [1,§4.2] weak-field) |
| **e1_stress_domega.py** | §5.1 | **[v2]** stress test (Gauss law breaks as sin saturates, Δω≳1) |
| **shapiro_inmodel_v1.py** | §5.3 | **[v2]** consistency check (eikonal on the lattice; factor 2 imported) |
| fig_overlay_forces.py | §4 | figure (v1) |
| fig_overlay_forces_v2.py | §4 | figure (final overlay) |
| **fig_protoforce_v1.py** | §5 | **[v2]** figure (four proto-force panels) |

## How to run

Each script is standalone, e.g.:
