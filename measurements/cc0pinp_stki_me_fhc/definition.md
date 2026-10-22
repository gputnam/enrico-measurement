# numu CC0pi Np single-transverse kinematic imbalance on CH, ME FHC (`cc0pinp_stki_me_fhc_*`)

Eight cross sections from one signal definition and one selection (one corpus entry each).
Units are cm^2/nucleon per unit of the variable (per GeV/c per deg for the double differential).

| id | observable | units | bins |
|---|---|---|---|
| cc0pinp_stki_me_fhc_dpt | dpT = \|pT_mu + pT_p\| | GeV/c | 0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.475, 0.55, 0.65, 0.80, 1.0 |
| cc0pinp_stki_me_fhc_dalphat | dalphaT = arccos(-pT_mu . dpT / (pT_mu dpT)) | deg | 0-180 in 15 deg steps |
| cc0pinp_stki_me_fhc_pn | pn = sqrt(dpT^2 + dpL^2), the reconstructed initial-neutron momentum | GeV/c | 0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.75, 1.0 |
| cc0pinp_stki_me_fhc_dphit | dphiT = arccos(-pT_mu . pT_p / (pT_mu pT_p)) | deg | 0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 130, 180 |
| cc0pinp_stki_me_fhc_pp | leading-proton momentum p_p | GeV/c | 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 0.90, 1.0, 1.2 |
| cc0pinp_stki_me_fhc_thetap | leading-proton polar angle theta_p (beam frame) | deg | 0, 10, 20, 25, 30, ..., 65, 70 |
| cc0pinp_stki_me_fhc_nprot | number of protons in the proton phase space: bin 1 = exactly 1, bin 2 = >= 2 | count | 0.5, 1.5, 2.5 |
| cc0pinp_stki_me_fhc_dpt_dalphat | dpT (fastest axis) in bins of dalphaT | GeV/c x deg | dpT 0, 0.1, 0.2, 0.3, 0.45, 1.0 x dalphaT 0, 45, 90, 135, 180 |

The transverse plane is the plane perpendicular to the neutrino beam direction. dpL comes from energy conservation,
as in Lu et al. 2018 (Eqs. 7-9):

- dpL = R/2 - (m_A'^2 + dpT^2)/(2R), with R = m_A + pL_mu + pL_p - E_mu - E_p.
- m_A = 6 m_n + 6 m_p - 92.16 MeV (the 12C mass).
- m_A' = m_A - m_n + b, with a fixed removal (excitation) energy **b = 27.13 MeV**.

## Signal (generator final state after FSI, status-1 particles)

- numu CC with the true vertex in the tracker FV (5980 < z < 8422 mm, apothem 850 mm).
- mu- with 1.5 < p_mu < 20 GeV/c and theta_mu < 20 deg, which is the corpus muon phase space. The LE analysis used
  p_mu < 10 GeV/c.
- No mesons (pi+-, pi0, K, other mesons), no heavy baryons, no leptons other than the muon, and no photon above
  10 MeV. This is the ccqelike_3d "QE-like" final state.
- At least one proton with 0.45 < p_p < 1.2 GeV/c and theta_p < 70 deg, as in Lu et al. 2018.
- The STKI proton is the most energetic proton satisfying the phase space. It is searched among the two leading
  final-state protons (`tfs_p_*_{0,1}`). The case where only the third or a later proton satisfies the phase space
  requires >= 3 protons above 450 MeV/c and is negligible.
- Neutrons and additional protons of any momentum are allowed.
- There is no cut on interaction mode, true W or Q2.
- nprot counts the protons in the same proton phase space, using the two leading slots (1 or >= 2).
- Truth: muon from `mc_lep`, proton from `tfs_p_*`, both rotated to the beam frame.

## Selection (reco)

Muon and vertex:
- MINOS-matched mu- (helicity 1) and reco vertex in the FV.
- theta_mu < 20 deg, 1.5 < p_mu < 20 GeV/c.
- At most 1 dead discriminator.

Proton tracks:
- At least 2 tracks, and every non-muon track is a proton candidate (`multiplicity == 2 + n_secp`).
- Every secondary candidate passes the Q2_QE-dependent CCQENu proton-score1 cut (>0.2/0.1/0.0 for Q2_QE
  <0.2/0.2-0.6/>0.6 GeV^2).
- The primary candidate (MasterAnaDev_proton_*, always the most energetic candidate) has score1 > 0.5. This is an
  ESC-like quality cut. It suppresses protons that reinteract or are misreconstructed, whose range-based momentum is
  biased low. The low tail p_reco/p_true < 0.9 drops from 23% to 15%. Lu et al. instead required elastically
  scattered contained protons, which is not available in the DF.
- The primary candidate is in the phase space: 0.45 < p_p < 1.2 GeV/c and theta_p < 70 deg. Its momentum comes from
  the dE/dx-range fit, and its direction is rotated to the beam frame.

Vetoes:
- No Michel electron (`n_improved_michel == 0`).
- Unattached recoil (`nonvtx_iso_blobs_E + dispr_id_E`) < 150 MeV.

nprot (reco): 1 + [any secondary candidate with 0.45 < secp_P < 1.2 GeV/c and secp_theta < 70 deg] (4 momentum slots).

Performance (full MC):
- Efficiency 0.24 overall. It is 0.24 at dpT < 0.3 GeV/c, falls to 0.14 at 0.8-1.0 GeV/c, is 0.20-0.24 across
  dalphaT, and 0.16-0.27 across pn.
- Purity 0.71. The backgrounds are 1 charged pion 16%, 1 pi0 6%, multi-pion 2%, and other (QE-like outside the
  proton/muon phase space, NC, wrong sign) 5%.
- 179.5k selected data events.

## Backgrounds

Template normalisations are fitted per universe by weighted least squares to the 2D (recoil 0-1.5 GeV + overflow,
Michel 0/>=1) distribution of events passing all other cuts. The free parameters are signal, the 1 charged pion
background and the pi0/multi-pion background. "other" is fixed to the MC. CV scales are 0.99 (signal), 0.94 (1
charged pion, which includes absorbed/missed pions) and 0.96 (pi0/multi-pion). One global set is used for all
observables, so the background shape in each variable comes from the MnvTune.

## Unfolding and systematics

- MnvTune v4.4.1 (with the FSI elastic weight), D'Agostini 6 iterations (4 for pp, thetap, nprot), per-universe
  migration.
- Standard universes: flux (nu-e constrained), GENIE incl. nucleon FSI (FrElas/FrInel/FrAbs/MFP_N), RPA, the
  low-recoil 2p2h tune, muon energy MINOS/MINERvA, MINOS efficiency, beam angle, GEANT proton/neutron/pion, target
  mass, and the recoil response.
- Proton universes. Energy scale (lateral, kinetic-energy shifts from the tuple): Bethe-Bloch, MEU, mass, Birks.
  Proton score1: Bethe-Bloch, MEU, mass, Birks.
- Typical total uncertainty 8-13%. Statistics are 1-2.5%. The largest sources are the proton-score universes (~4%
  each), flux 3.6%, 2p2h tune (up to 12% at dpT < 0.1 GeV/c), muon energy scale (up to 5% at small dalphaT) and the
  proton energy scale (up to 4% at small dalphaT).

## Closure

- MC-as-data with CV background: the closure is exactly 1. The closure printed by run_measurement (1.00-1.04)
  includes the data-fitted background scales.
- Reweighted-MC stress test, unfolding with the CV response:
  - nucleon FSI knobs and GEANT proton: <= 2%.
  - 2p2h-tune variant 3, which shifts the QE/2p2h mix by +25% at low dpT: up to +20% in the two lowest dpT bins
    (response composition dependence, because QE protons reinteract more). The Low_Recoil_2p2h_Tune band
    (12%/10% there) covers this.

## Motivation

- Model v5 has a global RFG (k_F = 221 MeV/c), RS17 2p2h with an np-dominated pair final state, and an Oset-type
  nucleon cascade.
- Its QE sector is constrained in the corpus only through muon kinematics, SumTp and Eavail. There, LFG(+RPA), a
  spectral function, RPA alone and a stiffer axial form factor are degenerate.
- STKI of the leading proton separates three things:
  - the initial-state nucleon momentum (the dpT and pn peak below ~0.25 GeV/c, and the Fermi edge);
  - intranuclear momentum transfer from 2p2h and pion absorption (the dpT/pn tail);
  - proton FSI (the dalphaT rise towards 180 deg, the dphiT tail).
- No proton angle and no muon-proton correlation is in the corpus.

## Model sensitivity

- Ground-state model: RFG vs LFG vs SF, through the pn shape near k_F and the high-momentum tail.
- Removal energy.
- Nucleon FSI strength and elastic/inelastic mix: dalphaT slope, theta_p, p_p.
- 2p2h strength and pair final state: dpT/pn tail, nprot >= 2.
- Pion absorption: tail.

## Caveats

- The proton momentum resolution (dpT resolution ~50 MeV/c core, with a high-side tail) is worse than in the LE ESC
  analysis. 50 MeV/c bins therefore have sizeable bin-to-bin anticorrelations: use the full covariance.
- The recoil cut makes the efficiency depend on the neutron/extra-proton content (2p2h): see the 2p2h band.
- Secondary-proton momenta are available for 4 slots and only the polar angle is stored, so nprot is coarse
  (purity of the >= 2 bin is 48%).

## References

- X.-G. Lu et al. (MINERvA), PRL 121, 022504 (2018), arXiv:1805.05486: signal, pn, b = 27.13 MeV.
- T. Cai et al. (MINERvA), PRD 101, 092001 (2020), arXiv:1910.08658: pn and dpL, binning ideas.
- X.-G. Lu et al., PRC 94, 015503 (2016), arXiv:1512.05748: STKI definitions.
- D. Ruterbories et al. (MINERvA), PRD 99, 012004 (2019), arXiv:1811.02774: CCQENu proton score and recoil cuts.
