# ccinc_eavail_q3_me_fhc: low-recoil d2sigma/dEavail dq3 (numu CC inclusive, ME FHC)

Follows MINERvA's ME low-recoil analysis, M.V. Ascencio et al., Phys. Rev. D 106, 032001 (2022),
arXiv:2110.13372. The LE precedent is P.A. Rodrigues et al., PRL 116, 071802 (2016), arXiv:1511.05944.
Data: ME FHC, 1.057e21 POT.

**Signal (topological, any hadronic final state):**
- numu CC (mc_incoming = 14, mc_current = 1) with the true vertex in the tracker FV (z 5980-8422 mm, apothem 850 mm);
- theta_mu < 20 deg and 1.5 < p_mu < 20 GeV/c (the published phase space).

There is no cut on W, interaction mode or hadron content. Events with q3 > 1.2 GeV or Eavail > 1.2 GeV are
signal events in the overflow; they are unfolded but not reported.

**Observables (truth, computed from generator final-state particles):**
- q3 [GeV] = |q| = sqrt(Q2 + q0^2), with q0 = Enu - E_mu (ev.q0q3_true, the MAT TruthFunctions convention).
- Eavail [GeV] = sum of T over protons + sum of T over pi+- + total E of every other final-state particle
  except neutrons (pi0, gamma, K, e, ...). Strange/charmed baryons contribute E - m_p. Neutrons contribute 0.
  This is `tfs_eavail` (makedf Derived.cxx). genpred applies the same definition to the status-1 particles.
  pi0s written undecayed by the generator therefore count with their total energy, the same as their decay photons.
  This is Eq. (1) of Ascencio et al.: the charged-pion mass is not included.

**Binning:**
- q3 [GeV]: 0, 0.2, 0.3, 0.4, 0.6, 0.9, 1.2 (the six panels of Ascencio et al.).
- Eavail [GeV]: 0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80, 1.0, 1.2.

That gives 102 bins, with Eavail varying fastest. The bins above the kinematic limit (Eavail > q3) are empty and are
not used: 85 bins have < 100% uncertainty. Cross section units: cm^2/nucleon/GeV^2.

**Reconstruction and selection:**
- Selection: MINOS-matched mu- (q/p < 0), reco vertex in the FV, theta_mu < 20 deg, 1.5 < p_mu < 20 GeV/c,
  <= 1 dead discriminator.
- q0_reco = calibrated calorimetric recoil energy of all non-muon clusters (MasterAnaDev_recoil_E, the CC-inclusive
  spline-corrected estimator; median reco/true q0 is about 1 from 0.05 to 2 GeV). Then Enu = E_mu + q0,
  Q2 = 2 Enu (E_mu - p_mu cos theta) - m_mu^2, and q3_reco = sqrt(Q2 + q0^2).
- Eavail_reco = 1.17 x (blob_recoil_E_tracker + blob_recoil_E_ecal): tracker + ECAL only, as in MAT LowRecoilFunctions.

The q3 resolution is about 20% (68% half-width) and the Eavail resolution is about 20/30/50/90/175 MeV for
q3 < 0.2/0.3-0.4/0.4-0.6/0.6-0.9/0.9-1.2 GeV. Both match Fig. 1 of the paper.

**Performance (full statistics):**
- 3.94M selected data events (all q3).
- Background (NC, wrong-sign, nue, out of phase space) is 3.6%, subtracted with MC.
- Efficiency is 0.66-0.88: 0.85 at q3 < 0.4 and 0.67 at q3 0.9-1.2.
- Closure is 1.000.
- Statistical uncertainty is 0.2-2%. Total uncertainty is typically 5-18% (median 11.6%).

**Unfolding and systematics:** D'Agostini with 4 iterations (the paper used 2). The response comes from MnvTune v1 plus the
GENIE FSI elastic and absorption fixes (= MnvTune-v1.2 of the paper). Universes: flux (100), GENIE knobs, RPA, low-recoil
2p2h tune, muon energy/angle, MINOS efficiency, Geant hadron reweights, target mass, and lateral recoil-response
universes that shift Eavail_reco and q0_reco together. The response universes cover proton, meson, EM, other,
low/mid/high-energy neutron and crosstalk; each uses the allNonMuonClusters response errors, scaled by 1.17 for
Eavail and by 1.51 (recoil_E / visible) for q0.

Dominant systematics:
- Lowest Eavail bin (< 20 MeV): crosstalk response (10-14%).
- Elsewhere: flux (3.5%), proton/meson response (1.5-8%), the 2p2h-tune and nucleon-FSI unfolding-prior bands
  (2-8%, largest in the dip), and muon energy scale (2%).
- Neutron response is <= 0.5% (median), up to 6% at high q3 and low Eavail.

**Motivation:** the inclusive (pT, p||) deficit (24%) and the QE-like deficit (41-46%) cannot be attributed from muon
kinematics alone. At fixed q3 the Eavail spectrum separates four regions:
- the QE peak: ground state, Pauli blocking and RPA at low q3;
- the dip between QE and Delta: 2p2h/MEC, whose energy goes partly into invisible neutrons and binding;
- the Delta peak: pion absorption moves Delta events down by ~m_pi;
- the high-W/DIS tail.

The lowest Eavail bin collects events with no visible hadronic energy (neutron-only final states, proton FSI and charge
exchange). No corpus entry constrains energy transfer.

**Model sensitivity:** the QE ground state (Fermi motion, Pauli blocking: an RFG proton has T > k_F^2/2M = 26 MeV),
2p2h, nucleon FSI (energy loss, p->n charge exchange), pion absorption, and the size of the RES/SIS/DIS strength at
q3 < 1.2 GeV.

**References:** Ascencio et al., arXiv:2110.13372; Rodrigues et al., arXiv:1511.05944; Ruterbories et al.
(ME inclusive, same muon selection), arXiv:2106.16210.
