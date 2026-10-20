# ccpip_me_fhc

```
MINERvA ME FHC numu CC1pi+ differential cross sections (tracked + untracked pions).

Reference: Granados, Messerly et al., arXiv:2605.24224 (code: github.com/MinervaExpt/CC-CH-pip-ana).
Variables: tpi (mixed tracked/untracked), thetapi (mixed, Tpi > 20 MeV signal), pmu, pzmu, ptmu,
thetamu, q2.

Signal (SignalDefinition::OnePi / OnePiThetaPi): numu CC, truth_is_fiducial, true vertex
5991.29 < z < 8400.91 mm, apothem 850 mm, theta_mu < 0.3491 rad, 1.5 < p_mu < 20 GeV/c,
W_exp(true) < 1.4 GeV, exactly one pi+ (Tpi_min < Tpi < 350 MeV), no other mesons, no photon
> 10 MeV, no e+-, any baryons.
Reco: vertex (same volume), MINOS-matched mu-, p_mu, theta_mu, < 2 isolated prongs, and either
  tracked : exactly one hadron with a quality endpoint Michel, LLR > 0 and node cuts, W_exp < 1.4
  untracked: good vertex/cluster Michel (LowRecoilPion), exactly one Michel candidate,
             Tpi_min < Tpi(range) < 350 MeV, pT_mu < 1.8 GeV/c, W_exp(Eavail) <= 1.4
The overlap uses the tracked reconstruction.  Background: W_exp sideband (reco W >= 1.5 GeV,
otherwise selected) fit of the true-W 1.4-1.8 and > 1.8 GeV background normalisations.

Not reproducible from public code: the MnvTune v4.7.1 (Tpi, Q2) charged-pion weight and the
"untracked pion" weight (GetChargedPionTuneWeight / GetUntrackedPionWeight are not public),
so the model is the v4.3.1 composition used by this analysis.
```

## Usage
```
venv/bin/python measurements/ccpip_me_fhc/run.py --mc <DF dir>/me_fhc_mc --data <DF dir>/me_fhc_data -j 16
```

## Results
Full ME FHC open data (DF tag v1): 1.057e21 data POT, 4.98e21 MC POT (489 files), all universes.

| variable | integral ratio | chi2 (ref + ours cov) / bins |
|---|---|---|
| Tpi | 0.982 | 12.9 / 10 |
| theta_pi | 1.044 | (no published covariance parsed) |
| p_mu | 1.029 | 0.6 / 8 |
| p_mu,parallel | 1.030 | 1.0 / 9 |
| pT_mu | 1.088 | 5.8 / 10 |
| theta_mu | 1.089 | 2.9 / 15 |
| Q2 | 1.065 | 4.2 / 12 |

W-sideband background scales (CV): W 1.4-1.8 GeV 0.70, W > 1.8 GeV 1.04.
Known residual differences: the two lowest Tpi bins (0-35 MeV) are 0.55 / 0.71 of the published
values (1.5-2 sigma) and the high pT / Q2 / theta_mu tails are 20-70% high (within uncertainties);
our systematic uncertainties are ~20% smaller than published.  These are consistent with the
missing MnvTune v4.7.1 charged-pion (Tpi, Q2) and untracked-pion weights, which are not public;
no ad-hoc substitute is applied.  Runtime 6.7 min for all 7 variables, peak RSS 1.6 GB.
