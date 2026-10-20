"""MINERvA ME FHC CC-inclusive double-differential cross section in (pT, p||).

Reference: D. Ruterbories, A. Filkins et al., Phys. Rev. D 104, 092007 (2021), arXiv:2106.16210
(the MINERvA-101 tutorial analysis).  MnvTune v1, 10 D'Agostini iterations.

Reco: tracker FV (z in [5980, 8422] mm, apothem 850 mm), theta_mu < 20 deg, MINOS-matched,
      <= 1 dead discriminator upstream, mu- (MINOS q/p < 0).
Signal: numu CC, true vertex in FV, theta_mu <= 20 deg, p|| >= 1.5 GeV/c.
"""
import numpy as np

from xsec.binning import Binning
from xsec.measurement import Measurement

PT = [0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5]
PZ = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60]
ZMIN, ZMAX, APOTHEM = 5980.0, 8422.0, 850.0


class CCInclusive2D(Measurement):
    name = "ccinc_2d_me_fhc"
    tune_name = "v1"
    n_iter = 10
    bkg_categories = ("NC", "wrong_sign", "nue", "out_of_phase_space")
    reco_binning = Binning([PT, PZ], ["pt", "pz"])
    true_binning = Binning([PT, PZ], ["pt", "pz"])

    def __init__(self, mat_signal_convention=False):
        # MINERvA-101 counts selected numu-CC events outside the truth phase space as signal in
        # the efficiency numerator / migration.  Default here: they are background.
        self.mat_signal_convention = mat_signal_convention

    def preselect(self, ev):
        return ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1)

    def reco(self, ev, shift):
        mu = ev.muon(shift)
        sel = (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (mu["theta"] < np.radians(20.0))
               & (ev["minos_match"] == 1) & (ev["dead"] <= 1) & (ev["minos_qp"] < 0))
        return (mu["pt"] / 1000.0, mu["pz"] / 1000.0), {"sel": sel}

    def truth(self, ev):
        t = ev.true_muon()
        return t["pt"] / 1000.0, t["pz"] / 1000.0

    def numu_cc(self, ev):
        return (ev["mc_incoming"] == 14) & (ev["mc_current"] == 1)

    def phase_space(self, ev):
        t = ev.true_muon()
        return ev.true_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (t["theta"] <= np.radians(20.0)) & (t["pz"] >= 1500.0)

    def signal(self, ev):
        if self.mat_signal_convention and not ev.is_truth:
            return self.numu_cc(ev)
        return self.numu_cc(ev) & self.phase_space(ev)

    def category(self, ev):
        cat = np.full(len(ev), 3)
        cat[ev["mc_current"] != 1] = 0
        cc = ev["mc_current"] == 1
        cat[cc & (ev["mc_incoming"] == -14)] = 1
        cat[cc & (np.abs(ev["mc_incoming"]) == 12)] = 2
        return cat
