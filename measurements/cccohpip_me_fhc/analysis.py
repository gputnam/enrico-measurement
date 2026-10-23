"""MINERvA ME FHC numu CC coherent-like pi+ production on the CH tracker.

Follows M.A. Ramirez et al. (MINERvA), PRL 131, 051801 (2023), arXiv:2210.01285 (ME, CH sample;
LE precedents: A. Higuera et al., PRL 113, 261802 (2014), arXiv:1409.3835; A. Mislivec et al.,
PRD 97, 032014 (2018), arXiv:1711.01178).

Signal (final-state topology, generator-computable): numu CC, true vertex in the tracker FV,
mu- + exactly one pi+ and nothing else hadronic: no other mesons (pi-, pi0, K, other), no nucleons of
any energy (tfs_np == tfs_nn == 0), no other baryons, no extra leptons, no photon > 10 MeV (nuclear
remnant / fragments allowed); 2 < E_nu < 20 GeV, theta_mu < 20 deg, theta_pi < 70 deg and
|t| < 0.125 GeV^2 with
    |t| = (pT_mu + pT_pi)^2 + (E_mu - pz_mu + E_pi - pz_pi)^2      (beam frame, nucleus recoil neglected)
computed from the true muon and pion four-momenta.  In nuflow v7 the 12C recoil is not a status-1
particle; in GENIE it is a nuclear fragment.  Neither enters the definition.

Reco: MINOS-matched mu- (q/p < 0), reco vertex in the tracker FV, <= 1 dead discriminator,
theta_mu < 20 deg, exactly two tracks at the vertex (multiplicity == 2), hadron track pion-like
(dE/dx LLR > 0), vertex energy E_vtx (non-muon visible energy within 150 mm of the vertex) < 50 MeV,
E_pi = calorimetric recoil energy of all non-muon clusters (MasterAnaDev_recoil_E, assumes
E_nu = E_mu + E_pi), pion direction from the track (pi_theta_0/pi_phi_0, beam frame),
2 < E_nu = E_mu + E_pi < 20 GeV, reco |t| < 0.125 GeV^2.
Background: incoherent processes (CC with no pion; CC with pions and W < 1.4 GeV; W >= 1.4 GeV) are
normalised by a fit to the high-|t| sideband (0.2 < |t| < 0.7 GeV^2, otherwise selected), binned in
|t| x E_pi; NC, wrong-sign, nue and out-of-phase-space coherent are taken from the MC.
"""
import numpy as np

from xsec.binning import Binning
from xsec.kinematics import M_MU, Shift, recoil_shift, rotate_to_beam
from xsec.measurement import Measurement
from xsec.universes import Universe, response_universes, standard_universes

ZMIN, ZMAX, APOTHEM = 5980.0, 8422.0, 850.0
M_PI = 139.57018
T_CUT = 0.125e6            # MeV^2
T_SIDE = (0.2e6, 0.7e6)    # MeV^2
EVTX_CUT = 50.0            # MeV, rec_vtx150
Q0_CALIB = 1.51            # recoil_E / visible energy, converts response errors to recoil_E units
RESPONSE = ("proton", "meson", "em", "other", "low_neutron", "mid_neutron", "high_neutron", "xtalk")
SIDE_T = [0.2, 0.3, 0.45, 0.7]          # GeV^2
SIDE_E = [0.0, 0.5, 1.0, 2.0, 20.0]     # GeV

VARIABLES = {
    "epi": dict(edges=[0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 6.0], n_iter=4),          # GeV
    "thetapi": dict(edges=[0, 5, 10, 15, 20, 25, 30, 70], n_iter=4),              # deg
    "enu": dict(edges=[2, 3, 4, 5, 6, 7, 8, 10, 14, 20], n_iter=4),                       # GeV
    "t": dict(edges=[0, 0.02, 0.04, 0.07, 0.125], n_iter=4),                              # GeV^2
    "q2": dict(edges=[0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 1.0], n_iter=4),  # GeV^2
    "pmu": dict(edges=[0, 2, 3, 4, 5, 6, 7, 8, 10, 20], n_iter=4),                    # GeV/c
    "thetamu": dict(edges=[0, 2, 4, 6, 8, 12, 20], n_iter=4),                             # deg
}


def abs_t(pxm, pym, pzm, Em, pxp, pyp, pzp, Ep):
    return (pxm + pxp) ** 2 + (pym + pyp) ** 2 + (Em - pzm + Ep - pzp) ** 2


class CCCohPip(Measurement):
    tune_name = "v4.4.1"
    tune_options = dict(lowq2pi=True)
    fv = (ZMIN, ZMAX, APOTHEM)
    samples = ("sel",)
    signal_sample = "sel"
    bkg_categories = ("nopi", "wlow", "whigh", "other")
    extra_reco_columns = ["multiplicity", "recoil_E", "rec_vtx150", "had_LLR_0", "pi_theta_0", "pi_phi_0",
                          "pi_theta_up_0", "pi_theta_dn_0"]
    extra_mc_reco_columns = [f"rerr_all_{p}" for p in RESPONSE]
    extra_truth_columns = ["t_is_fiducial", "tfs_np", "tfs_nn", "tfs_npip", "tfs_npim", "tfs_npi0", "tfs_nkpm",
                           "tfs_nk0", "tfs_nomeson", "tfs_nhbaryon", "tfs_nlep", "tfs_maxgammaE", "tfs_pi_pdg_0",
                           "tfs_pi_T_0", "tfs_pi_px_0", "tfs_pi_py_0", "tfs_pi_pz_0", "mc_targetA"]

    def __init__(self, variable="epi"):
        self.variable = variable
        v = VARIABLES[variable]
        self.name = f"cccohpip_me_fhc_{variable}"
        self.n_iter = v["n_iter"]
        self.reco_binning = Binning([v["edges"]], [variable])
        self.true_binning = Binning([v["edges"]], [variable])

    def universes(self, tune):
        U = standard_universes(tune) + response_universes(RESPONSE)
        for j, s in enumerate((-1, 1)):
            U.append(Universe("Pion_Angle", j, "lat", shift=Shift(f"PiTh{s:+d}", extra=(("pitheta", s),))))
        return U

    # ------------------------------------------------------------ reco
    def preselect(self, ev):
        return (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1) & (ev["multiplicity"] == 2)
                & (ev["had_LLR_0"] > 0))

    def _reco_all(self, ev, shift):
        mu = ev.muon(shift)
        com = (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1) & (ev["minos_qp"] < 0)
               & (ev["dead"] <= 1) & (mu["theta"] < np.radians(20.0)) & (ev["multiplicity"] == 2)
               & (ev["had_LLR_0"] > 0))
        dvis = recoil_shift(ev, shift, "rerr_all_")
        epi = np.maximum(ev["recoil_E"] + Q0_CALIB * dvis, M_PI)
        th = ev["pi_theta_0"].astype(np.float64)
        s = shift.get("pitheta", 0.0)
        if s and ev.is_mc:
            # MasterAnaDev_pion_theta_bias{Up,Down} are angle offsets (|delta| <= 6 mrad), not angles
            delta = ev["pi_theta_up_0"] if s > 0 else ev["pi_theta_dn_0"]
            th = th + np.where(np.isfinite(delta) & (np.abs(delta) < 0.1), delta, 0.0)
        ph = ev["pi_phi_0"].astype(np.float64)
        ppi = np.sqrt(np.maximum(epi * epi - M_PI * M_PI, 0.0))
        mphi = np.arctan2(np.tan(mu["thetaY"]), np.tan(mu["thetaX"]))
        t = abs_t(mu["pt"] * np.cos(mphi), mu["pt"] * np.sin(mphi), mu["pz"], mu["E"],
                  ppi * np.sin(th) * np.cos(ph), ppi * np.sin(th) * np.sin(ph), ppi * np.cos(th), epi)
        enu = mu["E"] + epi
        q2 = np.maximum(2.0 * enu * (mu["E"] - mu["p"] * np.cos(mu["theta"])) - M_MU * M_MU, 0.0)
        com = com & (enu > 2000.0) & (enu < 20000.0) & (th >= 0) & (th < np.pi)
        evtx = ev["rec_vtx150"]
        sel = com & (evtx < EVTX_CUT) & (t < T_CUT)
        side = com & (evtx < EVTX_CUT) & (t > T_SIDE[0]) & (t < T_SIDE[1])
        vals = dict(epi=epi / 1e3, thetapi=np.degrees(th), enu=enu / 1e3, t=t / 1e6, q2=q2 / 1e6,
                    pmu=mu["p"] / 1e3, thetamu=np.degrees(mu["theta"]))
        return vals, sel, side, t, epi

    def reco(self, ev, shift):
        vals, sel, _, _, _ = self._reco_all(ev, shift)
        return (vals[self.variable],), {"sel": sel}

    def extra(self, ev, shift):
        _, _, side, t, epi = self._reco_all(ev, shift)
        return {"tside": (Binning([SIDE_T, SIDE_E], ["t", "epi"], flow=False), (t / 1e6, epi / 1e3), side)}

    # ------------------------------------------------------------ truth
    def _true(self, ev):
        mu = ev.true_muon()
        px, py, pz = rotate_to_beam(ev["tfs_pi_px_0"], ev["tfs_pi_py_0"], ev["tfs_pi_pz_0"])
        p = np.sqrt(px * px + py * py + pz * pz)
        ok = (ev["tfs_pi_pdg_0"] == 211) & (p > 0)
        epi = np.where(ok, ev["tfs_pi_T_0"] + M_PI, np.nan)
        thpi = np.where(ok, np.degrees(np.arccos(np.clip(np.divide(pz, p, out=np.ones_like(p), where=p > 0), -1, 1))),
                        np.nan)
        mphi = np.arctan2(*rotate_to_beam(ev["mc_lep_0"], ev["mc_lep_1"], ev["mc_lep_2"])[1::-1])
        t = abs_t(mu["pt"] * np.cos(mphi), mu["pt"] * np.sin(mphi), mu["pz"], mu["E"], px, py, pz, epi)
        return mu, epi, thpi, t

    def truth(self, ev):
        mu, epi, thpi, t = self._true(ev)
        vals = dict(epi=epi / 1e3, thetapi=thpi, enu=ev["mc_Enu"] / 1e3, t=t / 1e6, q2=ev["mc_Q2"] / 1e6,
                    pmu=mu["p"] / 1e3, thetamu=np.degrees(mu["theta"]))
        return (vals[self.variable],)

    def _topology(self, ev):
        return ((ev["mc_incoming"] == 14) & (ev["mc_current"] == 1) & (ev["t_is_fiducial"] == 1)
                & ev.true_vertex_in_fv(ZMIN, ZMAX, APOTHEM)
                & (ev["tfs_npip"] == 1) & (ev["tfs_npim"] == 0) & (ev["tfs_npi0"] == 0)
                & (ev["tfs_nkpm"] + ev["tfs_nk0"] + ev["tfs_nomeson"] == 0)
                & (ev["tfs_np"] == 0) & (ev["tfs_nn"] == 0) & (ev["tfs_nhbaryon"] == 0)
                & (ev["tfs_nlep"] == 0) & (ev["tfs_maxgammaE"] <= 10.0))

    def signal(self, ev):
        mu, epi, thpi, t = self._true(ev)
        enu = ev["mc_Enu"]
        with np.errstate(invalid="ignore"):
            return (self._topology(ev) & (enu > 2000.0) & (enu < 20000.0) & (mu["theta"] < np.radians(20.0))
                    & (thpi < 70.0) & (t < T_CUT))

    def category(self, ev):
        """0 nopi: CC numu without pions; 1 wlow / 2 whigh: CC numu with pions (incoherent) and true
        W < / >= 1.4 GeV; 3 other: NC, wrong sign, nue, coherent outside the signal phase space."""
        cc = (ev["mc_incoming"] == 14) & (ev["mc_current"] == 1)
        npi = ev["tfs_npip"] + ev["tfs_npim"] + ev["tfs_npi0"]
        coh = ev["mc_intType"] == 4
        cat = np.full(len(ev), 3)
        cat[cc & ~coh & (npi == 0)] = 0
        cat[cc & ~coh & (npi > 0) & (ev["mc_W"] < 1400.0)] = 1
        cat[cc & ~coh & (npi > 0) & (ev["mc_W"] >= 1400.0)] = 2
        return cat


FIT_CATS = ("nopi", "wlow", "whigh")


def background_scales(acc, meas):
    """Fit the nopi / wlow / whigh normalisations (per universe) to the reco |t| x E_pi distribution in
    the high-|t| sideband (0.2 < |t| < 0.7 GeV^2).  Signal and 'other' are fixed to the MC."""
    d = acc["x/tside/data"]
    fixed = acc["x/tside/sig"] + acc["x/tside/bkg/other"]
    T = np.stack([acc[f"x/tside/bkg/{c}"] for c in FIT_CATS], 0)  # (C, bins, U)
    U = fixed.shape[1]
    sig2 = np.maximum(d, 1.0)
    scales = np.ones((len(FIT_CATS), U))
    for u in range(U):
        A = T[:, :, u].T / np.sqrt(sig2)[:, None]
        b = (d - fixed[:, u]) / np.sqrt(sig2)
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        scales[:, u] = np.clip(x, 0.0, None)
    R = acc["data/sel"].shape[0]
    return {c: np.tile(scales[i], (R, 1)) for i, c in enumerate(FIT_CATS)}, scales
