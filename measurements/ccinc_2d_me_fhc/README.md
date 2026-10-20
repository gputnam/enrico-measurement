# ccinc_2d_me_fhc

```
MINERvA ME FHC CC-inclusive double-differential cross section in (pT, p||).

Reference: D. Ruterbories, A. Filkins et al., Phys. Rev. D 104, 092007 (2021), arXiv:2106.16210
(the MINERvA-101 tutorial analysis).  MnvTune v1, 10 D'Agostini iterations.

Reco: tracker FV (z in [5980, 8422] mm, apothem 850 mm), theta_mu < 20 deg, MINOS-matched,
      <= 1 dead discriminator upstream, mu- (MINOS q/p < 0).
Signal: numu CC, true vertex in FV, theta_mu <= 20 deg, p|| >= 1.5 GeV/c.
```

## Usage
```
venv/bin/python measurements/ccinc_2d_me_fhc/run.py --mc <DF dir>/me_fhc_mc --data <DF dir>/me_fhc_data -j 16
```

## Results
(filled after the full-statistics run)
