# numu CC single-pi0 production on CH, ME FHC (`ccpi0_me_fhc_*`)

Six single-differential cross sections, one per corpus entry, all with the same signal definition and selection.
Units are cm^2/nucleon per unit of the variable.

| id | observable | units | bins |
|---|---|---|---|
| ccpi0_me_fhc_wexp | W_exp (hadronic invariant mass, struck nucleon at rest) | GeV | 1.0, 1.2, 1.35, 1.5, 1.65, 1.8 |
| ccpi0_me_fhc_q2 | Q2 = 2 Enu (E_mu - p_mu cos theta_mu) - m_mu^2 | GeV^2 | 0, 0.1, 0.2, 0.3, 0.45, 0.65, 0.9, 1.3, 2.0, 4.0 |
| ccpi0_me_fhc_tpi0 | pi0 kinetic energy T_pi0 | GeV | 0, 0.08, 0.16, 0.25, 0.35, 0.5, 0.7, 1.0, 1.5 |
| ccpi0_me_fhc_thetapi0 | pi0 polar angle w.r.t. the neutrino beam | deg | 0, 15, 30, 45, 60, 75, 90, 110, 135, 180 |
| ccpi0_me_fhc_pmu | muon momentum | GeV/c | 1.5, 2.5, 3.5, 4.5, 6, 8, 11, 20 |
| ccpi0_me_fhc_thetamu | muon polar angle w.r.t. the beam | deg | 0, 2, 4, 6, 8, 10, 12, 14, 17, 20 |

## Signal (generator final state, post-FSI, status-1 particles)
- numu CC with the true vertex in the tracker fiducial volume (5980 < z < 8422 mm, hexagon apothem 850 mm).
- Muon phase space: theta_mu < 20 deg and 1.5 < p_mu < 20 GeV/c. This is the same as ccinc/ccpip.
- **Exactly one pi0**, no pi+ or pi-, no K+-/K0, no other mesons, and no other charged or neutral leptons
  (`tfs_npi0 == 1`, `tfs_npip == tfs_npim == 0`, `tfs_nkpm + tfs_nk0 == 0`, `tfs_nomeson == 0`, `tfs_nlep == 0`).
  Any number of nucleons and photons is allowed.
- **W_exp < 1.8 GeV**, with W_exp^2 = M_N^2 + 2 M_N q0 - Q2, where q0 = Enu - E_mu, Q2 is the true four-momentum
  transfer and M_N = 938.9 MeV. This is the convention of ccpip_me_fhc and Altinok et al. It is not a cut on the
  generator's true W or interaction mode.
- The pi0 is the undecayed final-state pi0 of the generator record (`tfs_pi0_*_0`). The decay photons come from the
  detector simulation and are not part of the definition.
- Truth observables:
  - T_pi0 = `tfs_pi0_T_0`.
  - theta_pi0 is the angle of the `tfs_pi0` momentum after rotation to the beam frame.
  - Q2 = `mc_Q2`, and W_exp is as above.
  - The muon variables are taken from `mc_lep`, in the beam frame.
- Signal events with W_exp < 1.0 GeV (~1%, Fermi motion) or Q2 > 4 GeV^2, T_pi0 > 1.5 GeV (<1%) are in the signal
  but fall outside the reported bins of that variable.

## Selection (reco)
Muon and vertex:
- MINOS-matched mu- (q/p < 0).
- Reco vertex in the FV.
- At most 1 dead discriminator.
- theta_mu < 20 deg and 1.5 < p_mu < 20 GeV/c.

pi0 (MasterAnaDev CCpi0 / AngleScan reconstruction, McGivern et al. Sec. III, Altinok et al. Sec. III):
- exactly two shower candidates (`n_showers == 2`);
- both photons convert > 15 cm from the vertex (`g{1,2}_dist_vtx > 150` mm);
- at most one hadron track (`had_number <= 1`, i.e. a proton candidate);
- 75 < m_gg < 195 MeV;
- no Michel electron (`improved_nmichel == 0`; this vetoes charged pions);
- reconstructed W_exp < 1.8 GeV.

Reconstructed kinematics:
- q0_reco = recoil_E / 1.16. recoil_E is the calorimetric CC-inclusive recoil; the factor 1.16 is the median
  recoil_E / q0_true of selected signal, because the inclusive spline over-corrects EM energy.
- Enu = E_mu + q0, Q2 = 2 Enu (E_mu - p_mu cos theta_mu) - m_mu^2, and W_exp as above.
- Resolution (selected signal): W_exp core 0.14 GeV; Q2 IQR ±0.05 GeV^2; T_pi0 ~30% (median reco/true 0.94); theta_pi0
  ~8 deg.
- The pi0 momentum is p = p_g1 + p_g2, rotated to the beam frame, and T_pi0 = E_g1 + E_g2 - m_pi0.

**EM energy scale.** The m_gg peak was fit with a Gaussian plus a quadratic over 40-260 MeV, using events before the
m_gg/W/Michel cuts. The peak is at 129.7 ± 0.6 MeV in data and 124.7 ± 0.7 MeV in MC. The fit-range dependence is
±1%.
- Data photon energies are divided by 1.038. This follows the separate data/MC m_gg peak calibration of Altinok et al.
- The pi0 part of the data q0 is corrected by the same factor.
- A ±2% EM-scale lateral universe acts on m_gg, the pi0 kinematics and the pi0 part of q0.

## Backgrounds and fit
The background categories are:
- `pi0_meson`: ≥ 1 pi0 plus charged or other mesons;
- `pi0_other`: pi0 events outside the signal, i.e. W_exp > 1.8, ≥ 2 pi0, outside the muon phase space/FV, NC, or
  wrong sign;
- `charged`: no pi0, ≥ 1 charged or other meson;
- `other`: no mesons.

The first three are normalised by a weighted least-squares fit, repeated in every universe, to the m_gg distribution
(bins 0, 30, 60, 75, 195, 200, 250, 300, 400 MeV) in four regions: {Michel veto passed / failed} × {W_exp < 1.8 /
≥ 1.8 GeV}. The signal region itself is excluded. Signal and `other` are fixed to the MC.

Fit result (CV):
- pi0_meson 0.66, pi0_other 1.05, charged 1.30;
- the scale spread across universes is 0.03/0.09/0.07;
- chi2 = 1716 for 31 points: the m_gg shapes are not perfectly modelled.

Two alternative fits are carried as one-universe systematic bands:
- `BkgFit_noWside`: drops the W_exp ≥ 1.8 regions, giving 0.52/1.92/1.59;
- `BkgFit_noMichel`: drops the Michel regions with pi0_other fixed, giving 0.63/1/1.41.

Selected sample: 43,859 data events; post-fit MC purity 40%.

## Extraction
- D'Agostini unfolding with 4 iterations. Efficiency correction, standard flux (0-100 GeV, nu-e constrained), and CH
  nucleons in the FV.
- MnvTune v4.3.1 (deuterium pion tune, low-Q2 pion suppression, COH) with the FSI elastic/absorption fixes. It is
  used only for efficiency, migration and background templates.
- Systematics:
  - standard universes (flux, GENIE including Rvn/p1pi and Rvn/p2pi, FSI, MINOS, muon energy/angle, GEANT hadrons,
    target mass);
  - recoil response universes (proton, meson, EM, other, neutron, crosstalk) applied to q0;
  - the EM-scale universe;
  - the two background-fit bands.
- Overall efficiency is 5.8%. It is 5-6% versus W_exp, 4.4% → 7.7% from low to high T_pi0, and lowest (3.4-4%) at
  theta_pi0 60-110 deg.
- Closure (MC as data, unscaled backgrounds) is 1.000 in every bin. run_measurement.py prints a "closure" that applies
  the data-fitted background scales to MC pseudo-data, so it differs from 1 by construction.

## Motivation
This is the first neutral-pion measurement in the corpus.

- **Isospin.** CC1pi0 has mixed isospin. Delta(1232) production (I=3/2) gives p pi+ : p pi0 : n pi+ = 9:2:1. The I=1/2
  resonances N(1440), N(1520) and N(1535), and I=1/2 non-resonant strength, enhance pi0 relative to pi+.
- **No coherent pi0 in CC.** The low-Q2 / forward deficit seen in CC1pi+ cannot be coherent here.
- **Charge exchange.** pi+ → pi0 charge exchange and absorption in carbon feed and deplete this channel.

The corpus has only pi+ with W_exp < 1.4 GeV (ccpip_me_fhc_*). This measurement adds two things:
- the W_exp distribution across 1.0-1.8 GeV, separating the Delta peak from the second and third resonance regions;
- a pi0/pi+ comparison in the same muon phase space.

## Model sensitivity
- N-Delta axial coupling (C5A, MA_Delta).
- Higher resonances and non-resonant/DIS pion production at W 1.4-1.8 GeV, including the RES/DIS transition.
- Pion FSI: charge exchange and absorption, which shape the T_pi0 spectrum.
- Pi0 angular distribution (Delta decay versus DIS fragmentation).

## References
- O. Altinok et al. (MINERvA), Phys. Rev. D 96, 072003 (2017), arXiv:1708.03723: numu CC(pi0), W_exp < 1.8 GeV,
  m_gg window, EM scale from the m_gg peak, sideband fit.
- C.L. McGivern et al. (MINERvA), Phys. Rev. D 94, 052005 (2016), arXiv:1606.07127: CC pi+-/pi0 with W_exp < 1.8,
  photon conversion > 15 cm, 75 < m_gg < 195 MeV.
- T. Le et al. (MINERvA), Phys. Lett. B 749, 130 (2015), arXiv:1503.02107: anti-numu CC pi0.
- Granados, Messerly et al. (MINERvA), arXiv:2605.24224: ME CC1pi+ (corpus ccpip_me_fhc_*).
