# MvA-Loop: MINERvA Open-Data cross-section framework

A two-stage framework for neutrino cross-section measurements with the
[MINERvA Open Data Product](https://minerva.fnal.gov/opendata/), built without the MINERvA
analysis software (MAT/PlotUtils), which was used only as a reference for conventions.

```
makedf/        C++ (ROOT only) dataframe producer: AnaTuples -> flat, preselected DF trees
xsec/          python package: weights, tunes, systematic universes, unfolding, extraction
grid/          jobsub_lite scripts to run makedf over the full open data (see grid/README.md)
measurements/  one directory per measurement (analysis.py, run.py, reference/ published data)
```

## Setup
```
source makedf/scripts/setup_env.sh               # spack ROOT 6.28 + gcc 12 (C++ stage, aux converters)
cmake -S makedf -B makedf/build -DCMAKE_CXX_COMPILER=$(which g++) && cmake --build makedf/build -j8
python3 -m venv venv && venv/bin/pip install -r requirements.txt && venv/bin/pip install -e xsec
```
Auxiliary flux/reweight inputs are unpacked from `OpenData/FluxAndReweightFiles` into
`$MVA_AUX` (default `/exp/sbnd/data/users/gputnam/MvA-Loop/aux`) and converted once to npz:
```
source makedf/scripts/setup_env.sh; cd xsec
python3 -m xsec.aux.convert_flux; python3 -m xsec.aux.convert_reweight
```

## Stage 1: makedf
`makedf -c makedf/config -o out.root [-l filelist] files...` reads `MasterAnaDev`, `Truth`, `Meta`
over xrootd, enabling only the needed branches, and writes:
* `reco`  : one row per reconstructed candidate passing a loose preselection (vertex in a loose
  tracker volume, MINOS-matched muon, theta_mu < 25 deg); columns from `config/reco.cols`
  (+ `reco_mc.cols`, `truth.cols` for MC)
* `truth` : MC truth interactions (numu/numubar CC, loose FV, theta_mu < 25 deg)
* `meta`  : per input file POT and bookkeeping
Columns are listed in plain-text files (`name  branch | branch[i] | branch[:K]`), so the content
can change without recompiling.  Quantities that need variable-length inputs are computed in C++:
final-state truth summary (multiplicities, sum of proton KE, leading p/pi), GENIE elastic-FSI
fix weights, Geant hadron reweights, and the low-recoil-pion (Michel) reconstruction summary.
Full ME FHC (11 TB of MC, 1 TB of data) -> 36 GB of DF files (v1).

## Stage 2: xsec (python)
* `flux`: FluxReweighter-equivalent weights, nu-e constraint, universes, integrals
* `reweight`, `tunes`: MnvTune v1 / v4.x pieces (non-res pi, 2p2h fit, RPA, D2 pion tune,
  low-Q2 pion, z-expansion, elastic FSI fix, MINOS efficiency)
* `kinematics`: MINERvA muon reconstruction with lateral shifts
* `universes`: flux (100, nu-e weighted), GENIE, RPA, 2p2h, muon energy (MINERvA/MINOS +
  correlated flux), resolution, beam angle, MINOS efficiency, target mass, ...
* `runner`: parallel event loop over DF files (file-by-file, so memory stays small); MC scaled
  per playlist to the data POT; fills data, signal/background per universe, migration, efficiency
* `extract`: background subtraction (optional sideband fits), D'Agostini unfolding with analytic
  error propagation, efficiency, flux and target normalization, covariance per band
* closure test: `extract(..., closure=True)` unfolds the MC prediction and compares with the
  generator-level cross section from the truth tree

## Measurements
| directory | reference |
|---|---|
| `ccinc_2d_me_fhc` | CC inclusive d2sigma/dpT dp||, PRD 104 092007 (2021) |
| `ccqelike_3d_me_fhc` | QE-like d3sigma/dpT dp|| dSumTp, PRL 129 021803 (2022) |
| `ccqelike_3d_le_me_fhc` | QE-like 3D in the LE/ME paper binning, arXiv:2606.00745 |
| `ccpip_me_fhc` | CC1pi+, arXiv:2605.24224 |

Each `run.py --mc <DF dir> --data <DF dir> -j 16` runs the event loop, extracts the cross
section with all universes and prints comparisons (ratios, chi2 with the published covariance,
uncertainty breakdown) to the published result.
