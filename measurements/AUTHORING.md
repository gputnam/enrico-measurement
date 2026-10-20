# Authoring a new measurement (for the closed loop)

This is the recipe the experimentalist in `orchestrate/` follows.  A *measurement* is a directory
`measurements/<name>/` with

| file | purpose |
|---|---|
| `analysis.py` | a `xsec.measurement.Measurement` subclass: selection, binning, signal definition, truth observables |
| `corpus.json` | one entry per cross section the directory produces (id, class, kwargs, title, result path, motivation) |
| `definition.md` | human/agent-readable signal definition, observables + units, selection, backgrounds, motivation |

Start from `measurements/_template/` (`cp -r measurements/_template measurements/<name>`; rename the
class, `name`, the corpus `id` and `result`).  Existing measurements are worked examples of harder
selections: `ccqelike_3d_me_fhc` (proton score, Michel/blob samples, sideband fit), `ccpip_me_fhc`
(pion tracking + Michel-tagged untracked pions, W sideband), `ccinc_2d_me_fhc` (simplest).

## Interface (see `xsec/xsec/runner.py` docstring)

- `tune_name` / `tune_options`: MnvTune used for the efficiency, backgrounds and the GENIE prediction
  (`v1`, `v4.3.1`, `v4.4.1`, `v4.7.1`, `genie`; see `xsec/xsec/tunes.py`).
- `reco_binning`, `true_binning`: `xsec.binning.Binning(edges_list, names)`; usually identical.
  The first axis varies fastest in the flattened bin order.  Keep the number of bins modest
  (statistics: ~1e6 selected events for the inclusive sample; aim for >= a few hundred events per bin).
- `preselect(ev)`, `reco(ev, shift)`, `truth(ev)`, `signal(ev)`, `category(ev)`, `bkg_categories`,
  `samples`, `signal_sample`, optional `extra(ev, shift)` for sideband histograms and a module-level
  `background_scales(acc, meas)` for a background fit, optional `fv = (zmin, zmax, apothem)`.
- Units in `ev`: MeV, mm, rad.  `ev.muon(shift)` (reco, beam frame), `ev.true_muon()` (truth, beam frame),
  `ev.q0q3_true()`, `ev.reco_vertex_in_fv(...)`, `ev.true_vertex_in_fv(...)`.

## Hard rules

1. **Only existing DF columns.**  The dataframes (`/pnfs/sbnd/scratch/users/gputnam/MvA-Loop/df/v1/me_fhc_{mc,data}`)
   are fixed; new columns would need grid reprocessing of 11 TB of MC, which is out of scope.  Column
   names are the first field of `makedf/config/{reco,reco_mc,truth}.cols`; array columns like
   `vtx[:4]` become `vtx_0..vtx_3`, slot columns like `isoblob_z` become `isoblob_z_0..`.  Columns must
   be listed in `extra_reco_columns` / `extra_truth_columns` / `extra_mc_reco_columns` to be loaded.
2. **The signal definition and truth observables must be computable from generator output.**
   `signal()` and `truth()` are also evaluated on nuflow events by `xsec/xsec/genpred.py`, which provides
   only these truth columns: `mc_incoming, mc_current, mc_intType, mc_targetZ, mc_targetA,
   mc_targetNucleon, mc_Enu, mc_Q2, mc_W, mc_x, mc_y, mc_nu_*, mc_vtx_*, mc_lep_*, t_is_fiducial,
   t_N_pi*, tfs_*` (the final-state summary of `makedf/src/Derived.cxx`: particle counts, `tfs_sumTp`,
   `tfs_sumTn`, `tfs_eavail`, `tfs_maxgammaE`, leading two protons `tfs_p_{T,px,py,pz}_{0,1}` and charged
   pions `tfs_pi_{pdg,T,px,py,pz}_{0,1}`, leading nucleon), `t_*_mu, t_q0, t_q3`.  Define signals in terms
   of *final-state topology and kinematics* (what a detector sees), never in terms of the GENIE
   interaction type (`mc_intType`, `mc_resID`) or GENIE knob weights (`g_*`).  If you need a quantity
   that is a function of these columns (e.g. TKI from the leading proton and the muon) compute it in
   `truth()`.  If you genuinely need another generator-level column, add it to `genpred.load_events`
   with the same definition and units as the DF column.
3. **The measurement must be sound**: run_measurement prints a closure ratio (unfolding the MC
   prediction as data must recover the truth) that must be ~1, an efficiency, and the fractional
   uncertainty.  Bins with total fractional uncertainty >= 100% are ignored by the corpus metrics.
4. **ME FHC only** (numu, 1.057e21 POT).  The flux normalisation is the standard 0-100 GeV integral.

## Running

```bash
cd measurement
# quick development pass: one DF file per playlist (~15 s)
venv/bin/python measurements/run_measurement.py <name> --results-root <results> --max-files 1 -j 16
# full statistics (~1-7 min with 16 workers, peak ~2.5 GB)
venv/bin/python measurements/run_measurement.py <name> --results-root <results> -j 16
# re-extract only (e.g. after changing the background fit or n_iter)
venv/bin/python measurements/run_measurement.py <name> --results-root <results> --load
# add to a corpus and compare with a model
venv/bin/python -m xsec.corpus add --corpus <corpus> --measurement <name> --results-root <results> --origin iteration-<k>
venv/bin/python -m xsec.corpus predict --corpus <corpus> --model <nuflow.root> --model-version v<k> --only <id>
```

The corpus format is documented in `xsec/xsec/corpus.py`.
