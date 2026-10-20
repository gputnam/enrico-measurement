# CC1π⁺ (tracked + untracked π), ME FHC — reproduction status

**Reference:** Granados, Messerly et al., arXiv:2605.24224 (code: github.com/MinervaExpt/CC-CH-pip-ana).
**Status: LARGELY REPRODUCED — χ² acceptable in all variables with a published covariance; known residuals at low
Tπ and in high-pT/Q²/θμ tails, consistent with a model ingredient that is not public.**

## Inputs
DF tag v1, full ME FHC open data: 1.057e21 data POT, 4.98e21 MC POT, all universes. Model: MnvTune v4.3.1 composition
(see limitation below). W-sideband background fit: CV scales W 1.4–1.8 GeV 0.70, W > 1.8 GeV 1.04.

## Result vs published
| variable | integral ratio | χ² (ref + ours cov) / bins |
|---|---|---|
| Tπ | 0.982 | 12.9 / 10 |
| θπ | 1.044 | — (no published covariance parsed) |
| pμ | 1.029 | 0.6 / 8 |
| pμ,∥ | 1.030 | 1.0 / 9 |
| pT,μ | 1.088 | 5.8 / 10 |
| θμ | 1.089 | 2.9 / 15 |
| Q² | 1.065 | 4.2 / 12 |

## Known residuals
- Tπ 0–20 and 20–35 MeV bins: ratio 0.55 and 0.71 (1.5–2σ).
- High tails 20–75% above published, within uncertainties: pT > 0.8 GeV/c, Q² > 0.7 GeV², θμ > 12°.
- Our systematic uncertainties are ~20% smaller than published.

## Limitation (not fixable from public inputs)
The MnvTune **v4.7.1** charged-pion (Tπ, Q²) weight and the untracked-pion weight used in the paper
(`GetChargedPionTuneWeight` / `GetUntrackedPionWeight`) are not public. The residuals above are where these
weights act. No ad-hoc substitute has been applied; closing them would require the official weights.

## Files
- `v1/<var>.pkl`, `v1/<var>_xsec.pkl` for var ∈ {tpi, thetapi, pmu, pzmu, ptmu, thetamu, q2}; `v1/ccpip_v1.log` run log
  with per-bin tables.
- `dev/` — superseded: `ccpip_s2` small-scale test, `ccpip_v1pre` full-stat run before the per-universe-migration fix.

Reproduce: `venv/bin/python measurements/ccpip_me_fhc/run.py --load-dir outputs/ccpip_me_fhc/v1 -o <outdir>`
(full run: 6.7 min for all 7 variables, peak RSS 1.6 GB).
