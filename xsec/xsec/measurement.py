"""Base class for a measurement (see runner.py for the interface used by the event loop)."""
import numpy as np

from .universes import standard_universes

# columns every MINERvA reco-based measurement uses
MUON_COLUMNS = ["mu_lepE_0", "mu_lepE_1", "mu_lepE_2", "minos_p", "minos_used_curv", "minos_qp", "mu_thetaX",
                "mu_thetaY", "minos_match", "dead", "vtx_0", "vtx_1", "vtx_2", "numi_pot", "batch_structure",
                "vtx_batch", "has_vtx"]
TRUTH_COLUMNS = ["mc_incoming", "mc_current", "mc_intType", "mc_targetZ", "mc_targetNucleon", "mc_Enu", "mc_Q2",
                 "mc_W", "mc_vtx_0", "mc_vtx_1", "mc_vtx_2", "mc_lep_0", "mc_lep_1", "mc_lep_2", "mc_lep_3",
                 "g_Rvn1pi_m1", "g_Rvn1pi_p1", "g_Rvp1pi_m1", "g_Rvp1pi_p1", "g_MaRES_m1", "g_MaRES_p1"]
GENIE_KNOBS = ["AGKYxF1pi", "AhtBY", "BhtBY", "CCQEPauliSupViaKF", "CV1uBY", "CV2uBY", "EtaNCEL", "FrAbs_N",
               "FrAbs_pi", "FrCEx_N", "FrCEx_pi", "FrElas_N", "FrElas_pi", "FrInel_N", "FrInel_pi", "FrPiProd_N",
               "FrPiProd_pi", "MFP_N", "MFP_pi", "MaCCQE", "MaCCQEshape", "MaNCEL", "MaRES", "MvRES", "NormCCQE",
               "NormCCRES", "NormDISCC", "NormNCRES", "RDecBR1gamma", "Rvn1pi", "Rvn2pi", "Rvp1pi", "Rvp2pi",
               "Theta_Delta2Npi", "VecFFCCQEshape"]
GENIE_COLUMNS = [f"g_{k}_{s}" for k in GENIE_KNOBS for s in ("m1", "p1")]


class Measurement:
    name = "measurement"
    tune_name = "v1"
    tune_options = {}
    minos_p_offset = 0.0
    n_iter = 10
    samples = ("sel",)
    signal_sample = "sel"
    bkg_categories = ("bkg",)
    reco_binning = None
    true_binning = None
    extra_reco_columns = []
    extra_truth_columns = []

    @property
    def reco_columns(self):
        return sorted(set(MUON_COLUMNS + self.extra_reco_columns))

    @property
    def truth_columns(self):
        return sorted(set(TRUTH_COLUMNS + GENIE_COLUMNS + self.extra_truth_columns))

    def universes(self, tune):
        return standard_universes(tune)

    def preselect(self, ev):
        return np.ones(len(ev), dtype=bool)

    def category(self, ev):
        return np.zeros(len(ev), dtype=int)
