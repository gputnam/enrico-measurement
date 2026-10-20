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
(filled after the full-statistics run)
