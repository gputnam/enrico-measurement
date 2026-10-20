# CC-inclusive 2D (pT, p||), ME FHC — reproduction status

**Reference:** Ruterbories, Filkins et al., PRD 104 092007 (2021), arXiv:2106.16210 (MINERvA-101 analysis).
**Status: REPRODUCED** (with the published MINERvA-101 signal convention, `--mat-signal`).

## Inputs
DF tag v1, full ME FHC open data: 1.057e21 data POT (219/219 jobs), 4.98e21 MC POT (489 files),
all systematic universes. Flux integral 6.243e-8 /cm²/POT (ν-e constrained, 0–100 GeV; paper quotes 6.32e-8),
3.2353e30 tracker nucleons.

## Result vs published
| | `--mat-signal` (published convention) | strict convention (default) |
|---|---|---|
| integral ratio | 0.998 | 1.001 |
| bin ratio median / RMS | 0.9995 / 0.57% | 1.002 / 2.2% |
| χ² (published total cov) / bins | **4.5 / 205** | 41 / 205 |
| χ² (sum of covariances) | 2.6 / 205 | 10.8 / 205 |
| per-p|| ratios | 0.992–1.005 | 0.989–1.014 |

Fractional uncertainties (median over bins, ours / published):
stat 1.41% / 1.55%, flux 3.55% / 4.01%, total 5.79% / 6.86%.

## Caveats / open items
- **Signal convention.** The published numbers use the MINERvA-101 convention (all νμ CC events in migration and
  efficiency numerator, out-of-phase-space events removed via unfolding/efficiency). The strict convention subtracts them
  as MC background. Both are defensible; the difference is a few-% model-dependent correction at phase-space edges.
- **Muon energy scale band** matches the published `energyscale` covariance (~3.5%) for p|| > 2 GeV/c but is ~2× larger
  in the 1.5–2 GeV/c bin (11% vs 5%). Not yet understood.
- Total uncertainty is ~11% smaller than published (flux band ~12% smaller: our ν-e constrained flux).

## Files
- `v1/ccinc_v1{.pkl,_xsec.pkl,.log}` — final, `--mat-signal` (accumulators, extracted xsec+covariances, run log).
- `v1/ccinc_v1_strict{.pkl,_xsec.pkl,.log}` — final, strict convention.
- `dev/` — superseded development runs (`s1`, `s2`: small-scale tests; `v1pre*`: full-stat runs before the
  per-universe-migration fix). Kept for reference only.

Reproduce: `venv/bin/python measurements/ccinc_2d_me_fhc/run.py --load outputs/ccinc_2d_me_fhc/v1/ccinc_v1.pkl --mat-signal`
(full run: 68 s on 16 cores, peak RSS 2 GB).
