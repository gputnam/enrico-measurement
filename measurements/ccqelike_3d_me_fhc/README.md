# ccqelike_3d_me_fhc

```
MINERvA ME FHC QE-like triple-differential cross section d3sigma/(dpT dp|| dSumTp).

Reference: D. Ruterbories et al., Phys. Rev. Lett. 129, 021803 (2022), arXiv:2203.08022.
MnvTune v4.4.1, 10 D'Agostini iterations, background fit to 4 samples per (pT, SumTp) bin.

Signal: numu CC with only muon + nucleons (+ photons < 10 MeV) in the final state, true vertex
in the tracker FV, theta_mu < 20 deg, p|| > 1.5 GeV/c.  SumTp = summed true proton KE.
Reco: MINOS-matched mu- in FV, theta < 20 deg, <= 1 dead discriminator; every non-muon
track proton-like by dE/dx (Q2QE-dependent proton-score cut, CCQENu default).
Samples (Michel x isolated blobs): A (0, <=1) signal, B (>=1, <=1), C (0, >=2), D (>=1, >=2).
Recoil (nonvtx_iso_blobs + dispr_id) < 500 MeV; blobs: z > 4750 mm and < 10 MeV/hit.
Reco SumTp proxy: available energy 1.17 x (blob_recoil_E tracker + ecal).
```

## Usage
```
venv/bin/python measurements/ccqelike_3d_me_fhc/run.py --mc <DF dir>/me_fhc_mc --data <DF dir>/me_fhc_data -j 16
```

## Results
Full ME FHC open data (DF tag v1): 1.057e21 data POT, 4.98e21 MC POT (489 files), all universes.

| | ours | PRL 129 021803 |
|---|---|---|
| integral (sum of xsec x bin volume) | 5.798e-39 | 5.809e-39 (ratio 0.998) |
| bin ratio median | 0.968 | |
| chi2, diagonal (ours + ref) | 268 / 660 | |
| chi2 (sum of covariances) | 432 / 660 | |
| median frac. stat unc. | 7.1% | 8.0% |
| median frac. flux unc. | 4.00% | 4.00% |
| median frac. total unc. | 18.6% | 18.5% |

chi2 with the published total covariance alone is 4.5e3/660; that matrix is near singular (pinv) and
the value is dominated by its smallest eigenvalues, so the summed-covariance and diagonal values are
quoted.  Largest residual pulls: tiny cross sections at pT < 0.15 GeV/c and SumTp > 0.4 GeV.

Details needed for the reproduction that are not in the PRL text (see analysis.py docstring):
recoil < 500 MeV (PRD 99 012004), isolated blobs with < 10 MeV/hit (the PRL wording is inverted;
this reproduces the published sample sizes A/B/C/D = 1.30M/0.27M/0.21M/57k vs 1.3M/0.23M/0.22M/57k),
reco SumTp = MAT available energy 1.17 x (tracker + ECAL blob recoil), background fit with p|| bins
as separate fit points, negative background-subtracted bins set to 0 before unfolding.
Runtime 3.7 min on 16 cores, peak RSS 2.6 GB per process.
