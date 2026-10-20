# plot/: comparison figures

Draw each reproduced cross section (`outputs/<measurement>/v1/*_xsec.pkl`) with the published MINERvA result
overlaid, plus ratio panels.  Output goes to `figures/<measurement>/` (git-ignored), as PDF and PNG.

```
venv/bin/python plot/plot_all.py                 # all measurements
venv/bin/python plot/plot_qe3d.py --format png   # one measurement; --input <xsec pkl>, --outdir <dir>
venv/bin/python plot/plot_ccinc_2d.py --variant strict   # CC-inclusive with the strict signal convention
```

| script | figures |
|---|---|
| `plot_ccinc_2d.py` | `xsec_pt_in_pz` (d²σ/dp_T dp_∥ per p_∥ bin), `xsec_pz_projection`, `frac_unc` |
| `plot_qe3d.py` | `xsec_pz<k>` (d³σ vs ΣT_p per p_T bin, one figure per p_∥ bin), `projections`, `frac_unc` |
| `plot_qe3d_leme.py` | `xsec_coarse_pz<k>`: LE/ME-paper binning vs PRL 129 projected onto common bins (no numerical release for that paper) |
| `plot_ccpip.py` | `xsec_<var>` for the 7 CC1π⁺ variables, `summary` |

Blue points: this reproduction (capped bar: statistical, full bar: total uncertainty).  Orange line and band:
published result with its total uncertainty.  χ² and integral ratios printed on the figures are computed with
the same masks and covariances as each measurement's `run.py`.  Reference parsing is reused from
`measurements/*/run.py` and `measurements/ccpip_me_fhc/reference_tables.py`.
