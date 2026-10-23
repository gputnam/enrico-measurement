# cccohpip_me_fhc: numu CC coherent-like pi+ production on CH (MINERvA ME FHC)

Corpus ids: `cccohpip_me_fhc_{epi,thetapi,enu,t,q2,pmu,thetamu}` (one unfolded single-differential
cross section each, full covariance; the seven share the same selection and signal).

**Reference analysis:** M.A. Ramirez et al. (MINERvA), "Neutrino-induced coherent pi+ production in C,
CH, Fe and Pb at <E_nu> ~ 6 GeV", PRL 131, 051801 (2023), arXiv:2210.01285 (CH tracker sample only).
LE precedents: A. Higuera et al., PRL 113, 261802 (2014), arXiv:1409.3835; A. Mislivec et al.,
PRD 97, 032014 (2018), arXiv:1711.01178.

## Signal (final-state topology, computable on generator status-1 particles)
numu CC interaction, true vertex in the tracker FV (z 5980-8422 mm, apothem 850 mm), with
- exactly one pi+ and no other meson (no pi-, pi0, K+-, K0, other mesons);
- **no nucleons of any energy** (tfs_np = tfs_nn = 0), no other baryons, no extra leptons, no photon with
  E > 10 MeV (the nuclear remnant / fragments are allowed: in nuflow v7 the 12C recoil is not a status-1
  particle, in GENIE it is a nucleus);
- 2 < E_nu < 20 GeV, theta_mu < 20 deg, theta_pi < 70 deg (zero tracking efficiency above);
- |t| < 0.125 GeV^2 with |t| = |sum_{mu,pi} pT|^2 + (sum_{mu,pi} (E - pL))^2 (beam frame, nucleus recoil
  neglected; arXiv:2210.01285 Eq. 1), from the true muon and pion four-momenta.

This is the "mu- + exactly one pi+ + no other status-1 hadrons, |t| < cut" definition. It is the
coherent final state without a reference to the interaction type. In the MnvTune GENIE sample it is
99.7% coherent. In nuflow v7, 96.6% of the signal cross section is coherent (process 100, which
genpred maps to mc_intType 4), 2.3% is RES and 1.1% DIS/SIS. Those are incoherent events that leave
no final-state nucleon. 95.6% of the v7 coherent cross section passes the signal phase space.

## Observables (true = reco definitions; units)
| id suffix | observable | bins |
|---|---|---|
| epi | E_pi, total pion energy [GeV] | 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 6 |
| thetapi | theta_pi w.r.t. the beam [deg] | 0, 5, 10, 15, 20, 25, 30, 70 |
| enu | E_nu [GeV] (reco: E_mu + E_pi) | 2, 3, 4, 5, 6, 7, 8, 10, 14, 20 |
| t | \|t\| [GeV^2] | 0, 0.02, 0.04, 0.07, 0.125 |
| q2 | Q^2 = -(p_nu - p_mu)^2 [GeV^2] (reco with E_nu = E_mu + E_pi) | 0, .025, .05, .075, .1, .15, .2, .3, .4, .6, 1.0 |
| pmu | p_mu [GeV/c] | 0, 2, 3, 4, 5, 6, 7, 8, 10, 20 |
| thetamu | theta_mu [deg] | 0, 2, 4, 6, 8, 12, 20 |

All are dsigma/dX in cm^2 / nucleon / unit, normalised to the standard 0-100 GeV flux integral
(AUTHORING rule 4). The `enu` entry is therefore the flux-integrated dsigma/dE_nu, not sigma(E_nu).
sigma(E_nu) = (dsigma/dE_nu * Phi_tot) / (dPhi/dE_nu) in each bin. The model/data ratio per bin is
the same either way. E_pi and theta_pi bins follow the published CH tables. theta_pi is merged
above 30 deg, where the signal vanishes.

## Selection
MINOS-matched mu- (q/p < 0), reco vertex in the FV, <= 1 dead discriminator, theta_mu < 20 deg;
exactly two tracks at the vertex (`multiplicity == 2`); the hadron track is pion-like
(dE/dx LLR `had_LLR_0 > 0`); vertex energy `rec_vtx150` (visible non-muon energy within 150 mm of the
vertex) < 50 MeV (vetoes nucleons and extra hadrons; the published analysis used a cylinder of 200 mm
radius and 7 planes, which includes the muon, and a +-1 sigma window); E_pi = calorimetric recoil of all
non-muon clusters (`recoil_E`, median E_pi^reco/E_pi^true = 0.97, about 25% resolution);
pion direction from the track (`pi_theta_0`, `pi_phi_0`); 2 < E_nu = E_mu + E_pi < 20 GeV;
reco |t| < 0.125 GeV^2.

## Backgrounds
Incoherent backgrounds are split into three categories, each normalised by a fit:
- CC numu with no pion;
- CC numu with pions and true W < 1.4 GeV;
- CC numu with pions and W >= 1.4 GeV.

The fit uses the high-|t| sideband (0.2 < |t| < 0.7 GeV^2, otherwise fully selected), binned in |t|
(0.2/0.3/0.45/0.7) x E_pi (0/0.5/1/2/20 GeV). It is redone in every systematic universe. The fitted
normalisations do not depend on the generator's coherent/incoherent split. The remaining categories
come from the MC: NC, wrong-sign, nue, and coherent events outside the signal phase space. Signal in
the sideband is about 3% and is fixed to the MC.

## Extraction
MnvTune v4.4.1 + low-Q^2 pion suppression (the published tune) is used for the efficiency, the
migration and the background templates. The MINERvA coherent reweight is not applied. D'Agostini
unfolding with 4 iterations. Systematics: the standard universes (flux, GENIE, RPA, 2p2h, muon
energy, MINOS efficiency, beam angle, Geant4 hadron, target mass), the recoil-response lateral
universes (which shift E_pi), and a pion-angle lateral universe (MasterAnaDev pion_theta bias
offsets).

## Motivation and model sensitivity
nuflow v7 added CC coherent pi+: Berger-Sehgal PCAC, pi-12C elastic data below T_pi = 0.87 GeV, and
the Rein-Sehgal hadronic model with PDG piN input above 1.7 GeV. It closed the low-Q^2 deficit of
CC1pi+. Two deficits remain: low pT / high pz in ccinc_2d, and forward pions in ccpip_thetapi
(T_pi < 350 MeV only). Both could come either from coherent production being too small at
E_nu > 5 GeV / E_pi > 1 GeV (hypothesis A) or from deficient low-Q^2 incoherent DIS (hypothesis B).
This measurement isolates coherent production directly. Its E_pi, theta_pi and E_nu dependence
tests the high-energy (Rein-Sehgal absorptive) extrapolation. |t| tests the pi-12C slope, and Q^2
and the muon kinematics test the PCAC form factor. No other corpus entry has a nucleon-free,
low-|t| selection.
