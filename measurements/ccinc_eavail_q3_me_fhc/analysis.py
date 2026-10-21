"""MINERvA ME FHC numu CC-inclusive low-recoil double-differential cross section d2sigma/(dEavail dq3).

Follows M.V. Ascencio et al. (MINERvA), Phys. Rev. D 106, 032001 (2022), arXiv:2110.13372
(LE precedent: P.A. Rodrigues et al., PRL 116, 071802 (2016), arXiv:1511.05944).

Signal: numu CC, true vertex in the tracker FV, theta_mu < 20 deg, 1.5 < p_mu < 20 GeV/c (any hadrons).
True observables: q3 = sqrt(Q2 + q0^2), q0 = Enu - Emu (MAT TruthFunctions convention, ev.q0q3_true);
Eavail = sum T(p) + sum T(pi+-) + sum E(any other final-state particle except neutrons; pi0/gamma/K
total energy, strange baryons minus a nucleon mass) = makedf tfs_eavail (Derived.cxx), the same
definition genpred applies to generator final-state (status-1) particles.
Reco: MINOS-matched mu- in FV, theta < 20 deg, 1.5 < p < 20 GeV/c, <= 1 dead discriminator.
  q0_reco  = calorimetric (spline-calibrated) recoil energy of all non-muon clusters (MasterAnaDev_recoil_E,
             the CC-inclusive hadronic energy estimator), Enu = Emu + q0, Q2 = 2 Enu (Emu - p cos th) - m_mu^2,
             q3_reco = sqrt(Q2 + q0^2).
  Eavail_reco = 1.17 x (blob_recoil_E_tracker + blob_recoil_E_ecal)  (MAT LowRecoilFunctions GetEAvailable:
             tracker + ECAL only, which gives the better resolution at low hadronic energy).
Tune: MnvTune v1 + FSI elastic and absorption fixes (= MnvTune-v1.2 of the paper).  Recoil-response
lateral universes (proton, meson, EM, other, neutron, crosstalk) shift both Eavail and q0.
"""
import numpy as np

from xsec.binning import Binning
from xsec.kinematics import M_MU, recoil_shift
from xsec.measurement import Measurement
from xsec.universes import response_universes, standard_universes

ZMIN, ZMAX, APOTHEM = 5980.0, 8422.0, 850.0
Q3 = [0.0, 0.2, 0.3, 0.4, 0.6, 0.9, 1.2]  # GeV (Ascencio et al. panels)
EAVAIL = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80,
          1.0, 1.2]  # GeV
EAVAIL_SCALE = 1.17  # MAT LowRecoilFunctions Eavailable_scale
Q0_CALIB = 1.51      # recoil_E / visible (all non-muon clusters, ID+OD): converts response errors to q0 units
RESPONSE = ("proton", "meson", "em", "other", "low_neutron", "mid_neutron", "high_neutron", "xtalk")


class LowRecoilEavailQ3(Measurement):
    name = "ccinc_eavail_q3_me_fhc"
    tune_name = "v1"
    tune_options = dict(fsi_elastic=True, fsi_absorption=True)
    n_iter = 4
    fv = (ZMIN, ZMAX, APOTHEM)
    samples = ("sel",)
    signal_sample = "sel"
    bkg_categories = ("NC", "wrong_sign", "nue", "out_of_phase_space")
    reco_binning = Binning([EAVAIL, Q3], ["eavail", "q3"])
    true_binning = Binning([EAVAIL, Q3], ["eavail", "q3"])
    extra_reco_columns = ["recoil_E", "blob_recoil_trk", "blob_recoil_ecal"]
    extra_mc_reco_columns = [f"rerr_all_{p}" for p in RESPONSE]
    extra_truth_columns = ["tfs_eavail", "w_fsi_absorption"]

    def universes(self, tune):
        return standard_universes(tune) + response_universes(RESPONSE)

    def preselect(self, ev):
        return ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1)

    def reco(self, ev, shift):
        mu = ev.muon(shift)
        sel = (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (mu["theta"] < np.radians(20.0))
               & (ev["minos_match"] == 1) & (ev["dead"] <= 1) & (ev["minos_qp"] < 0)
               & (mu["p"] > 1500.0) & (mu["p"] < 20000.0))
        dvis = recoil_shift(ev, shift, "rerr_all_")
        vis = np.maximum(ev["blob_recoil_trk"], 0.0) + np.maximum(ev["blob_recoil_ecal"], 0.0)
        eavail = np.maximum(EAVAIL_SCALE * (vis + dvis), 0.0)
        q0 = np.maximum(ev["recoil_E"] + Q0_CALIB * dvis, 0.0)
        enu = mu["E"] + q0
        q2 = 2.0 * enu * (mu["E"] - mu["p"] * np.cos(mu["theta"])) - M_MU * M_MU
        q3 = np.sqrt(np.maximum(q2 + q0 * q0, 0.0))
        return (eavail / 1000.0, q3 / 1000.0), {"sel": sel}

    def truth(self, ev):
        _, q3 = ev.q0q3_true()
        return np.maximum(ev["tfs_eavail"], 0.0) / 1000.0, q3 / 1000.0

    def signal(self, ev):
        t = ev.true_muon()
        return ((ev["mc_incoming"] == 14) & (ev["mc_current"] == 1) & ev.true_vertex_in_fv(ZMIN, ZMAX, APOTHEM)
                & (t["theta"] < np.radians(20.0)) & (t["p"] > 1500.0) & (t["p"] < 20000.0))

    def category(self, ev):
        cat = np.full(len(ev), 3)
        cat[ev["mc_current"] != 1] = 0
        cc = ev["mc_current"] == 1
        cat[cc & (ev["mc_incoming"] == -14)] = 1
        cat[cc & (np.abs(ev["mc_incoming"]) == 12)] = 2
        return cat
