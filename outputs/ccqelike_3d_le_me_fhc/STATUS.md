# QE-like 3D in LE/ME-paper binning, ME FHC — reproduction status

**Reference:** Ruterbories et al., arXiv:2606.00745 (LE vs ME QE-like 3D).
**Status: PARTIAL — ME beam done, consistent with PRL; no direct numerical comparison possible; LE beam pending.**

## Why partial
- The paper has **no numerical data release**; only figures. A quantitative check is made by projecting our result and the
  PRL 129 result (same ME data) onto common coarse bins.
- The paper does not list ΣTp bin edges; the PRL ΣTp binning is used.
- **LE open-data samples are not yet released** (OpenData/LowEnergy_* empty, "coming soon"). The LE configuration
  (playlists, LE flux + ν-e constraint, GENIE 2.8.4-based tune) is prepared as `LE_CONFIG` but not run.

## Result (ME FHC, DF v1, 1.057e21 POT, all universes)
Compared with the PRL result projected onto common bins (pT as in the LE/ME paper, p|| 1.5–3.5 and 3.5–4.5 GeV/c):

| | value |
|---|---|
| integral ratio to PRL (coarse bins) | 0.950 |
| p|| 1.5–3.5 / 3.5–4.5 ratio | 0.957 / 0.943 |
| χ² (sum of covariances) / bins | 215 / 200 |

The ~5% lower normalization relative to the PRL-binned extraction on the same data is attributed to the background fit
being integrated over p|| (one fit point per sample and cell), as the LE/ME paper prescribes.
Qualitative comparison against the paper's figures has not been done (no plots produced yet).

## Files
- `v1/qele_v1{.pkl,_xsec.pkl,.log}` — final ME result (accumulators, extracted xsec+covariances, run log).
- `dev/` — superseded: `qe3d_leme_s2*` small-scale test, `qele_v1pre*` full-stat run before the QE 3D selection fixes.

Reproduce: `venv/bin/python measurements/ccqelike_3d_le_me_fhc/run.py --load outputs/ccqelike_3d_le_me_fhc/v1/qele_v1.pkl`
(full run: 3.1 min on 16 cores, peak RSS 1.7 GB).
