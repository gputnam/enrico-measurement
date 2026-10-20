# QE-like 3D (pT, p||, ΣTp), ME FHC — reproduction status

**Reference:** Ruterbories et al., PRL 129 021803 (2022), arXiv:2203.08022 (660-bin release with covariances).
**Status: REPRODUCED** in normalization and uncertainties; bin-level agreement within uncertainties
(some residual pulls in sparsely populated corners).

## Inputs
DF tag v1, full ME FHC open data: 1.057e21 data POT, 4.98e21 MC POT, all universes. MnvTune v4.4.1,
10 D'Agostini iterations, 4-sample background fit per (pT, ΣTp) cell with p|| bins as separate fit points.

## Result vs published
| | ours | PRL |
|---|---|---|
| integral (Σ xsec × bin volume) | 5.798e-39 | 5.809e-39 (**ratio 0.998**) |
| bin ratio median | 0.968 | |
| χ² diagonal (ours + ref) / bins | 268 / 660 | |
| χ² (sum of covariances) / bins | 432 / 660 | |
| median frac. stat / flux / total unc. | 7.1% / 4.00% / 18.6% | 8.0% / 4.00% / 18.5% |
| selected samples A / B / C / D | 1.30M / 0.27M / 0.21M / 57k | 1.3M / 0.23M / 0.22M / 57k |

Projections: per-p|| ratios 0.96–0.98; per-pT ratios 0.97–0.99 below 0.85 GeV/c, 0.87–0.90 above;
per-ΣTp ratios 0.88–1.13. Background-fit CV scales (median): signal 1.03, 1π± 1.15, 1π⁰ 1.01, multi-π 0.49.

## Caveats / open items
- χ² with the published total covariance alone is 4.5e3/660: that matrix is near-singular (pinv) and the value is
  dominated by its smallest eigenvalues, so the diagonal and summed-covariance χ² are the meaningful figures.
- Largest pulls: very small cross sections at pT < 0.15 GeV/c and ΣTp > 0.4 GeV.
- Selection details not stated (or stated inverted) in the PRL, chosen because they reproduce the published sample sizes:
  recoil (non-vertex isolated blobs + dispersed) < 500 MeV (from PRD 99 012004); isolated blobs counted only if
  z > 4750 mm and E/nhits **< 10 MeV** (PRL wording is the opposite); reco ΣTp axis = MAT available energy
  1.17 × (tracker + ECAL blob recoil).
- Negative background-subtracted bins are set to 0 before unfolding.
- v1 DFs contain ΣTp = −0.004 MeV for some events (proton-mass rounding); clamped to 0 in the analysis, fixed in makedf
  for future productions.

## Files
- `v1/qe3d_v1{.pkl,_xsec.pkl,.log}` — final (accumulators, extracted xsec+covariances, run log).
- `dev/` — superseded: `s2*` small-scale test, `v1pre*` full-stat iterations before the selection/fit fixes
  (`_b/_c/_d` logs, `nofit` = no background fit, `t*` = test variants).

Reproduce: `venv/bin/python measurements/ccqelike_3d_me_fhc/run.py --load outputs/ccqelike_3d_me_fhc/v1/qe3d_v1.pkl`
(full run: 3.7 min on 16 cores, peak RSS 2.6 GB).
