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
Reco SumTp proxy: calorimetric visible recoil (nonvtx_iso_blobs_energy + dispr_id_energy).
```

## Usage
```
venv/bin/python measurements/ccqelike_3d_me_fhc/run.py --mc <DF dir>/me_fhc_mc --data <DF dir>/me_fhc_data -j 16
```

## Results
(filled after the full-statistics run)
