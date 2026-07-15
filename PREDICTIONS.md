# Predictions register (prediction-before-run discipline)

Reviewer recommendation (v2): make "predictions written before the run" an
auditable fact, not an assertion on trust. Each new run script carries its
prediction in the header. To turn this into a git-verifiable record, commit
the prediction ahead of the result — the recommended workflow:

```
# 1. write the script with ONLY the header prediction + setup, no result print
git add <script>.py && git commit -m "predict: <one line>"   # prediction hash

# 2. run, paste the measured result into the script, commit again
git add <script>.py && git commit -m "result: <one line>"     # result hash
```

The two commit hashes then bracket the prediction in time. Below is the
register for the v2 runs; each `predict:` line is what was written before the
corresponding run.

| Run | predict: (before running) | outcome |
|---|---|---|
| `trough_unscreened_v1.py` | 1D unscreened potential is linear → confinement; attractor r*≈10; no λ so the r*(λ) question dissolves | r*=10.6, linear tail confirmed, plateau over ×16 source |
| `field_from_dynamics_v1.py` | locked lag solves Kψ″=Ω−ω; Gauss jump = Δω/K; no screening length (front grows with T) | jump 0.0497/0.0500; λ_fit grows 86→224 |
| `e1_stress_domega.py` | Gauss law holds while linearised; breaks as sin saturates at large Δω | holds <1% to Δω=0.5; 43% error at Δω=2.0 |
| `tau_vs_delay_v1.py` | heavier template → shorter life; light templates below barrier stable | 416→213→7→4 (single seed); 1 stable at Δω=0.70 |
| `tau_multiseed.py` | the single-seed slope will carry large scatter; decayed fraction should be monotone | decayed fraction monotone; median scatter factor 2–5 |
| `shapiro_inmodel_v1.py` | transit lag = ∫(e^2φ−1)dx (eikonal), linear in mass | ratio 1.001, lag doubles with mass |

Note on honesty: `field_from_dynamics_v1.py` and `shapiro_inmodel_v1.py` are
**consistency checks** — their "predictions" are near-guaranteed by the
linearisation / eikonal respectively. The register records this rather than
dressing them as risky forecasts. The genuinely open-outcome runs are
`e1_stress_domega.py` (where does it break) and `tau_multiseed.py` (is the
single-seed slope real).
