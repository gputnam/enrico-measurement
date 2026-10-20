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
(filled after the full-statistics run)
