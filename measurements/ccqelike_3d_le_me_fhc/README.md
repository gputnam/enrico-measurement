# ccqelike_3d_le_me_fhc

```
QE-like triple-differential cross section in the binning of the LE/ME comparison paper.

Reference: D. Ruterbories et al., arXiv:2606.00745 (LE vs ME QE-like 3D).  Same selection and
method as PRL 129, 021803 (see ../ccqelike_3d_me_fhc), with p|| restricted to 1.5-4.5 GeV/c and
the background fit integrated over that range.

Binning: p|| and pT edges are from the paper text / figure panel labels; the SumTp edges are not
given in the paper, so the PRL SumTp binning is used (to be revisited if tables become public).
The LE beam configuration (playlists minerva1/7/9/13A-E, LE flux, GENIE 2.8.4-based tune without
the elastic-FSI fix) is prepared via `LE_CONFIG` but the LE open-data samples are not yet released.
```

## Usage
```
venv/bin/python measurements/ccqelike_3d_le_me_fhc/run.py --mc <DF dir>/me_fhc_mc --data <DF dir>/me_fhc_data -j 16
```

## Results
Full ME FHC open data (DF tag v1): 1.057e21 data POT, 4.98e21 MC POT (489 files), all universes.

No numerical release exists; the result is projected onto bins common with the PRL binning
(pT as in the LE/ME paper, p|| in [1.5, 3.5] and [3.5, 4.5] GeV/c) and compared with the PRL
result projected the same way:

| | value |
|---|---|
| integral ratio to PRL (coarse bins) | 0.950 |
| p|| 1.5-3.5 / 3.5-4.5 ratio | 0.957 / 0.943 |
| chi2 (sum of covariances) | 215 / 200 |

The ~5% lower normalisation relative to the PRL-binned extraction on the same data comes from the
background fit integrated over p|| (one point per sample and cell, 4 parameters), as described in
the LE/ME paper.  LE beam: pending release of the LE open-data samples (`LE_CONFIG` prepared).
