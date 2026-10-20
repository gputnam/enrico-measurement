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
venv/bin/python measurements/ccinc_2d_me_fhc/run.py --mc <DF dir>/me_fhc_mc --data <DF dir>/me_fhc_data -j 16 --mat-signal   # published convention
```

## Results
Full ME FHC open data (DF tag v1): 1.057e21 data POT, 4.98e21 MC POT (489 files), all universes.

| | ours | PRD 104 092007 |
|---|---|---|
| integrated cross section ratio | 0.998 | |
| bin ratio median / RMS | 0.9995 / 0.57% | |
| chi2 (published total cov) | 4.5 / 205 | |
| chi2 (sum of covariances) | 2.6 / 205 | |
| median frac. stat unc. | 1.41% | 1.55% |
| median frac. flux unc. | 3.55% | 4.01% |
| median frac. total unc. | 5.79% | 6.86% |

Flux integral (nu-e constrained, 0-100 GeV) 6.243e-8 /cm^2/POT (paper quotes 6.32e-8), 3.2353e30 nucleons.

The published result uses the MINERvA-101 signal convention (`--mat-signal`: all numu CC events
enter migration and efficiency numerator, events outside the phase space are removed via the
efficiency/flow-bin unfolding).  With the strict convention (out-of-phase-space events subtracted as
MC background, default) the agreement is chi2 41/205, bin-ratio RMS 2.2%, integral 1.001.  Both are
valid; the difference is a model-dependent correction of a few percent at the phase-space edges.

Muon energy scale band matches the published `energyscale` covariance (~3.5%) for p|| > 2 GeV/c;
it is ~2x larger in the 1.5-2 GeV/c bin (11% vs 5%).  Runtime 68 s on 16 cores, peak RSS 2 GB.
