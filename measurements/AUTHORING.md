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

- `tune_name` / `tune_options`: the MINERvA simulation tune (MnvTune) used for the detector-level
  ingredients of the extraction: efficiency, migration, backgrounds
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

1. **Only existing DF columns.**  The dataframes (`/pnfs/sbnd/scratch/users/gputnam/MvA-Loop/df/v2/me_fhc_{mc,data}`;
   v2 = v1 + the pi0/photon/shower reconstruction, see "pi0 / photon columns" below) are fixed; new
   columns would need grid reprocessing of 11 TB of MC, which is out of scope.  Column
   names are the first field of `makedf/config/{reco,reco_mc,truth}.cols`; array columns like
   `vtx[:4]` become `vtx_0..vtx_3`, slot columns like `isoblob_z` become `isoblob_z_0..`.  Columns must
   be listed in `extra_reco_columns` / `extra_truth_columns` / `extra_mc_reco_columns` to be loaded.
   **Report schema gaps explicitly.**  If a measurement you want needs a quantity that is not a DF
   column, check `makedf/config/SOURCE_BRANCHES.md` (all source AnaTuple branches *not* in the DFs)
   and say in your report which it is: *blocked by the DF schema* (the branch exists in the source
   tuples; name it, it could be added in a future DF version) or *blocked by the data* (no such
   information in the MINERvA open data).  Do not silently replace the quantity by a proxy.
2. **The signal definition and truth observables must be computable from generator output.**
   `signal()` and `truth()` are also evaluated on nuflow events by `xsec/xsec/genpred.py`, which provides
   only these truth columns: `mc_incoming, mc_current, mc_intType, mc_targetZ, mc_targetA,
   mc_targetNucleon, mc_Enu, mc_Q2, mc_W, mc_x, mc_y, mc_nu_*, mc_vtx_*, mc_lep_*, t_is_fiducial,
   t_N_pi*, tfs_*` (the final-state summary of `makedf/src/Derived.cxx`: particle counts, `tfs_sumTp`,
   `tfs_sumTn`, `tfs_eavail`, `tfs_maxgammaE`, leading two protons `tfs_p_{T,px,py,pz}_{0,1}`, charged
   pions `tfs_pi_{pdg,T,px,py,pz}_{0,1}` and neutral pions `tfs_pi0_{T,px,py,pz}_{0,1}`, leading
   final-state photon `tfs_gamma_{E,px,py,pz}_0`, leading nucleon), `t_*_mu, t_q0, t_q3`.
   The MINERvA-simulation-only truth columns `t_pi0_*`, `t_gamma_*`, `t_N_gamma_pi0` and the
   truth-matching columns `t_g1_*`, `t_g2_*`, `t_shower_pdg_*` must NOT be used in `signal()` / `truth()`
   (they may be used in `category()` or for studies).  Define signals in terms
   of *final-state topology and kinematics* (what a detector sees), never in terms of the GENIE
   interaction type (`mc_intType`, `mc_resID`) or GENIE knob weights (`g_*`).  If you need a quantity
   that is a function of these columns (e.g. TKI from the leading proton and the muon) compute it in
   `truth()`.  If you genuinely need another generator-level column, add it to `genpred.load_events`
   with the same definition and units as the DF column.
3. **The measurement must be sound**: run_measurement prints a closure ratio (unfolding the MC
   prediction as data must recover the truth) that must be ~1, an efficiency, and the fractional
   uncertainty.  Bins with total fractional uncertainty >= 100% are ignored by the corpus metrics.
4. **ME FHC only** (numu, 1.057e21 POT).  The flux normalisation is the standard 0-100 GeV integral.
5. **Compare only with the model under study.**  The simulation's own truth-level prediction (GENIE /
   MnvTune) is an extraction tool, not a benchmark: `run_measurement.py` strips it from the result
   (it goes to the human-only store), and the corpus carries only data and the model.  Do not compute,
   print or quote it (e.g. from the accumulators), nor other generators' predictions from
   `reference/` files or papers, as a target for the model.

## pi0 / photon columns (DF v2)

The MasterAnaDev tuples carry MINERvA's CCpi0 reconstruction (AngleScan shower finding + two photon
candidates; arXiv:1606.07127 Sec. III, Altinok/McGivern CC pi0 analyses).  DF names (full list and
source branches: `makedf/config/reco.cols`, section "v2 additions"):

| group | columns | notes |
|---|---|---|
| showers | `n_showers` (AngleScan shower candidates), `n_good_blobs`, `shower_{E,nclus,ndig,good,score}_{0..3}`, `pi0_prefilter`, `pi0_ml_evt`, `vtx_min_oa` | `pi0_prefilter` (tuple pi0 prefilter flag) is 1 or -1 |
| pi0 | `pi0_mgg` (gamma-gamma invariant mass), `pi0_{E,T,P,px,py,pz}`, `pi0_theta`, `pi0_phi`, `pi0_thetaX`, `pi0_thetaY`, `pi0_oa`, `pi0_cos_oa` | MeV; theta/phi in rad w.r.t. the **detector** z axis (rotate p for the beam frame); `pi0_thetaX/Y` and the opening angle `pi0_oa` in **degrees**; -9999 when no candidate |
| photons | `g{1,2}_{E,px,py,pz,theta,phi}`, `g{1,2}_dir_{0..2}` (unit vector), `g{1,2}_vtx_{0..2}` (conversion point, mm), `g{1,2}_dist_vtx` (conversion distance, mm), `g{1,2}_dEdx` (dE/dx at shower start), `g{1,2}_time` (ns), `g{1,2}_maxZ`, `g{1,2}_evis_{trk,ecal,scalX,scalUV}`, `g{1,2}_E_{trk,ecal,scalX,scalUV}` | g1 is the more energetic photon in ~90% of events; pi0 p = g1 + g2 |
| other | `had_em_frac`, `p_em_frac`, `e_unattached`, `e_dispersed`, `n_muon_iso_blobs`, `muon_iso_blobs_E`, `vtx_iso_blobs_E_out`, `n_prongs` | |
| truth matching (MC) | `t_g{1,2}_{pdg,parent_pdg,parent_E,E,px,py,pz,frac,evis}`, `t_shower_pdg_{0..3}` | Geant4 particle (and parent) behind each photon candidate |
| truth (MC) | `tfs_pi0_{T,px,py,pz}_{0,1}` (signal-safe), `tfs_gamma_{E,px,py,pz}_0`; `t_pi0_{E,px,py,pz,theta}_{0,1}`, `t_gamma_*`, `t_N_gamma_pi0` (simulation only) | `tfs_*` momenta in the detector frame like the other `tfs_` slots |

A CC pi0 selection in the spirit of the published MINERvA analyses: the standard muon/vertex
preselection; `n_showers == 2` (exactly two shower candidates; `>= 2` for CC Npi0);
`g1_dist_vtx > 150` and `g2_dist_vtx > 150` mm (photons convert >= 15 cm from the vertex, rejecting
charged-pion and proton activity); `75 < pi0_mgg < 195` MeV (the mass peak is at ~135 MeV in data
and MC; use the sidebands for the background constraint); optionally no charged-pion candidate / no
Michel (`improved_nmichel == 0`) for CC1pi0 and a reconstructed W cut.  On 3 MC files this simple
selection has ~76% purity for events with >= 1 true pi0.  `g?_dEdx` (photon dE/dx at the start of
the shower) separates photons from protons; `pi0_prefilter` is the tuple's own pi0 prefilter.  The
photon energies are calorimetric; the published analyses calibrated the scale by fitting the mgg
peak separately in data and MC, so compare the data and MC peak positions before relying on it.

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
