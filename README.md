# MvA-Loop: MINERvA Open-Data cross-section framework

Two-stage framework for cross-section measurements with the MINERvA Open Data Product.

- `makedf/`  — C++ (ROOT via spack) producer: MasterAnaDev/Truth/Meta tuples → flat, preselected dataframe files (grid-friendly).
- `xsec/`    — python package: selection, weights/tunes, systematic universes, unfolding, cross-section extraction.
- `grid/`    — jobsub_lite submission scripts for running makedf over the full Open Data.
- `measurements/` — one directory per measurement (config, selection, run script, published references).

## Setup
```
source makedf/scripts/setup_env.sh       # spack ROOT (C++ stage, aux converters)
python3 -m venv venv && venv/bin/pip install -r requirements.txt && venv/bin/pip install -e xsec
```
Auxiliary flux/reweight files: `/exp/sbnd/data/users/gputnam/MvA-Loop/aux/` (unpacked from OpenData/FluxAndReweightFiles).
