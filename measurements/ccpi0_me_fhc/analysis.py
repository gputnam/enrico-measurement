"""MINERvA ME FHC numu CC single-pi0 production on CH: dsigma/dX for X = W_exp, Q2, T_pi0, theta_pi0,
p_mu, theta_mu (one corpus entry per variable).

Signal definition follows MINERvA's numu CC(pi0) analysis (O. Altinok et al., PRD 96, 072003 (2017),
arXiv:1708.03723; anti-numu precedent: C.L. McGivern et al., PRD 94, 052005 (2016), arXiv:1606.07127),
in the muon phase space of the ME corpus (ccpip_me_fhc / ccinc):
  numu CC, true vertex in the tracker FV (z 5980-8422 mm, apothem 850 mm),
  theta_mu < 20 deg, 1.5 < p_mu < 20 GeV/c,
  exactly one final-state (post-FSI, status 1) pi0, no pi+-, no K, no other mesons, no other
  leptons; any number of nucleons / photons,
  W_exp < 1.8 GeV with W_exp^2 = M_N^2 + 2 M_N q0 - Q2, q0 = Enu - E_mu (struck nucleon at rest;
  the same convention as ccpip_me_fhc).
The pi0 is the undecayed final-state pi0 of the generator record (tfs_pi0_*_0); its decay photons
are produced by the detector simulation and are not part of the signal definition.

Reco (MasterAnaDev CCpi0 reconstruction, AngleScan shower finding):
  MINOS-matched mu-, reco vertex in FV, <= 1 dead discriminator, theta_mu < 20 deg, 1.5 < p_mu < 20 GeV/c;
  exactly two shower candidates, both converting > 15 cm from the vertex (McGivern et al.);
  at most one hadron track (the proton candidate); 75 < m_gg < 195 MeV; no Michel electron
  (improved_nmichel == 0, charged-pion veto); W_exp(reco) < 1.8 GeV.
  q0_reco = recoil_E / 1.16 (calorimetric CC-inclusive recoil, rescaled so that median q0_reco/q0_true = 1
  for selected CC1pi0 signal; the CC-inclusive spline over-corrects EM showers), Enu = E_mu + q0,
  Q2 = 2 Enu (E_mu - p_mu cos theta_mu) - m_mu^2, W_exp as above.
  pi0: p = g1 + g2 (detector frame, rotated to the beam frame), T = E_gg - m_pi0.
EM energy scale: the m_gg peak (Gaussian + quadratic fit, 40-260 MeV, events before the m_gg cut) is
  at 129.7 MeV in data and 124.7 MeV in MC (full data, 1/5 of MC): data photon energies are divided by
  EM_DATA_SCALE = 1.038 (as in Altinok et al., who fit the peak separately in data and MC), and an
  EM-energy-scale lateral universe of +-2% is applied to the MC photon energies (m_gg, pi0 kinematics,
  and the pi0 part of q0).
Backgrounds: categories pi0_meson (>= 1 pi0 + charged/other mesons), pi0_other (pi0 events outside the
  signal: W_exp > 1.8, >= 2 pi0, outside the muon phase space/FV, NC), charged (no pi0, >= 1 charged
  meson), other (no meson).  The first three are normalised by a fit to the m_gg sidebands
  (< 60 MeV, 200-400 MeV) and to the Michel-tagged and W_exp >= 1.8 GeV samples (fit_regions below).
"""
import numpy as np

from xsec.binning import Binning
from xsec.kinematics import M_MU, Shift, recoil_shift, rotate_to_beam
from xsec.measurement import Measurement
from xsec.universes import Universe, response_universes, standard_universes

ZMIN, ZMAX, APOTHEM = 5980.0, 8422.0, 850.0
M_N = (938.272 + 939.565) / 2.0
M_PI0 = 134.977
W_MAX = 1800.0
MGG_LO, MGG_HI = 75.0, 195.0
DIST_MIN = 150.0          # mm, photon conversion distance from the vertex
Q0_SCALE = 1.16           # recoil_E / q0_true for selected CC1pi0 signal (MC)
Q0_CALIB = 1.51           # recoil_E / visible energy (converts response errors to q0 units, as ccinc_eavail_q3)
EM_DATA_SCALE = 1.038     # data/MC m_gg peak ratio
EM_SCALE_ERR = 0.02
RESPONSE = ("proton", "meson", "em", "other", "low_neutron", "mid_neutron", "high_neutron", "xtalk")

VARIABLES = {
    "wexp": dict(edges=[1.0, 1.2, 1.35, 1.5, 1.65, 1.8], n_iter=4),                       # GeV
    "q2": dict(edges=[0.0, 0.1, 0.2, 0.3, 0.45, 0.65, 0.9, 1.3, 2.0, 4.0], n_iter=4),       # GeV^2
    "tpi0": dict(edges=[0.0, 0.08, 0.16, 0.25, 0.35, 0.5, 0.7, 1.0, 1.5], n_iter=4),        # GeV
    "thetapi0": dict(edges=[0, 15, 30, 45, 60, 75, 90, 110, 135, 180], n_iter=4),          # deg
    "pmu": dict(edges=[1.5, 2.5, 3.5, 4.5, 6, 8, 11, 20], n_iter=4),                       # GeV/c
    "thetamu": dict(edges=[0, 2, 4, 6, 8, 10, 12, 14, 17, 20], n_iter=4),                  # deg
}
# background-fit histogram: m_gg x region (0 signal-like, 1 Michel, 2 W>=1.8, 3 Michel & W>=1.8)
FIT_MGG = [0.0, 30.0, 60.0, MGG_LO, MGG_HI, 200.0, 250.0, 300.0, 400.0]
FIT_REGION = [-0.5, 0.5, 1.5, 2.5, 3.5]
FIT_CATS = ("pi0_meson", "pi0_other", "charged")
# alternative background fits (systematic bands): without the W_exp >= 1.8 GeV rows (regions 2, 3);
# without the Michel-tagged rows (regions 1, 3) and with pi0_other fixed to the MC (it is then unconstrained)
BKGFIT_BANDS = ("BkgFit_noWside", "BkgFit_noMichel")


def _wexp(q2, q0):
    w2 = M_N * M_N - q2 + 2.0 * M_N * q0
    return np.where(w2 > 0, np.sqrt(np.maximum(w2, 0.0)), 0.0)


def _theta_deg(px, py, pz):
    p = np.sqrt(px * px + py * py + pz * pz)
    return np.degrees(np.arccos(np.clip(np.divide(pz, p, out=np.ones_like(p), where=p > 0), -1, 1)))


def em_universes():
    return [Universe("EM_energy_scale", j, "lat", shift=Shift(f"EM{s:+d}", extra=(("em_scale", s),)))
            for j, s in enumerate((-1, 1))]


class CCPi0(Measurement):
    tune_name = "v4.3.1"
    tune_options = dict(fsi_absorption=True)
    fv = (ZMIN, ZMAX, APOTHEM)
    samples = ("sel",)
    signal_sample = "sel"
    bkg_categories = ("pi0_meson", "pi0_other", "charged", "other")
    extra_reco_columns = ["recoil_E", "improved_nmichel", "n_showers", "had_number", "pi0_mgg", "pi0_E", "pi0_px",
                          "pi0_py", "pi0_pz", "g1_dist_vtx", "g2_dist_vtx"]
    extra_mc_reco_columns = [f"rerr_all_{p}" for p in RESPONSE]
    extra_truth_columns = ["t_is_fiducial", "tfs_npip", "tfs_npim", "tfs_npi0", "tfs_nkpm", "tfs_nk0", "tfs_nomeson",
                           "tfs_nlep", "tfs_pi0_T_0", "tfs_pi0_px_0", "tfs_pi0_py_0", "tfs_pi0_pz_0", "w_fsi_absorption",
                           "tfs_pi_T_0", "tfs_pi_px_0", "tfs_pi_py_0", "tfs_pi_pz_0"]  # tfs_pi_*: MnvTune COH weight

    def __init__(self, variable="wexp"):
        self.variable = variable
        v = VARIABLES[variable]
        self.name = f"ccpi0_me_fhc_{variable}"
        self.n_iter = v["n_iter"]
        self.reco_binning = Binning([v["edges"]], [variable])
        self.true_binning = Binning([v["edges"]], [variable])

    def universes(self, tune):
        # the two BkgFit_* universes are copies of the CV whose background scales come from alternative
        # sideband fits (see background_scales); they must stay the last two universes
        return (standard_universes(tune) + response_universes(RESPONSE) + em_universes()
                + [Universe(b, 0, "vert") for b in BKGFIT_BANDS])

    # ------------------------------------------------------------ reco
    def preselect(self, ev):
        return (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1) & (ev["n_showers"] == 2)
                & (ev["g1_dist_vtx"] > DIST_MIN) & (ev["g2_dist_vtx"] > DIST_MIN) & (ev["had_number"] <= 1))

    def _reco_all(self, ev, shift):
        mu = ev.muon(shift)
        com = (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1) & (ev["minos_qp"] < 0)
               & (ev["dead"] <= 1) & (mu["theta"] < np.radians(20.0)) & (mu["p"] > 1500.0) & (mu["p"] < 20000.0)
               & (ev["n_showers"] == 2) & (ev["g1_dist_vtx"] > DIST_MIN) & (ev["g2_dist_vtx"] > DIST_MIN)
               & (ev["had_number"] <= 1))
        s_em = (1.0 + EM_SCALE_ERR * shift.get("em_scale")) if ev.is_mc else 1.0 / EM_DATA_SCALE
        mgg = ev["pi0_mgg"] * s_em
        epi = np.maximum(ev["pi0_E"], 0.0) * s_em
        px, py, pz = rotate_to_beam(ev["pi0_px"] * s_em, ev["pi0_py"] * s_em, ev["pi0_pz"] * s_em)
        dvis = recoil_shift(ev, shift, "rerr_all_")
        q0 = np.maximum((ev["recoil_E"] + Q0_CALIB * dvis) / Q0_SCALE + (s_em - 1.0) * epi, 0.0)
        enu = mu["E"] + q0
        q2 = np.maximum(2.0 * enu * (mu["E"] - mu["p"] * np.cos(mu["theta"])) - M_MU * M_MU, 0.0)
        w = _wexp(q2, q0)
        nomich = ev["improved_nmichel"] == 0
        win = (mgg > MGG_LO) & (mgg < MGG_HI)
        sel = com & win & nomich & (w < W_MAX)
        vals = dict(wexp=w / 1e3, q2=q2 / 1e6, tpi0=np.maximum(epi - M_PI0, 0.0) / 1e3, thetapi0=_theta_deg(px, py, pz),
                    pmu=mu["p"] / 1e3, thetamu=np.degrees(mu["theta"]))
        region = (~nomich).astype(float) + 2.0 * (w >= W_MAX)
        fitmask = com & (mgg > 0.0) & (mgg < FIT_MGG[-1])
        return vals, sel, mgg, region, fitmask

    def reco(self, ev, shift):
        vals, sel, _, _, _ = self._reco_all(ev, shift)
        return (vals[self.variable],), {"sel": sel}

    def extra(self, ev, shift):
        _, _, mgg, region, fitmask = self._reco_all(ev, shift)
        return {"fit": (Binning([FIT_MGG, FIT_REGION], ["mgg", "region"], flow=True), (mgg, region), fitmask)}

    # ------------------------------------------------------------ truth
    def wexp_true(self, ev):
        return _wexp(ev["mc_Q2"], ev["mc_Enu"] - ev["mc_lep_3"])

    def truth(self, ev):
        t = ev.true_muon()
        px, py, pz = rotate_to_beam(ev["tfs_pi0_px_0"], ev["tfs_pi0_py_0"], ev["tfs_pi0_pz_0"])
        has = ev["tfs_npi0"] >= 1
        vals = dict(wexp=self.wexp_true(ev) / 1e3, q2=ev["mc_Q2"] / 1e6,
                    tpi0=np.where(has, ev["tfs_pi0_T_0"] / 1e3, np.nan),
                    thetapi0=np.where(has, _theta_deg(px, py, pz), np.nan),
                    pmu=t["p"] / 1e3, thetamu=np.degrees(t["theta"]))
        return (vals[self.variable],)

    def _topology(self, ev):
        return ((ev["tfs_npi0"] == 1) & (ev["tfs_npip"] == 0) & (ev["tfs_npim"] == 0)
                & (ev["tfs_nkpm"] + ev["tfs_nk0"] == 0) & (ev["tfs_nomeson"] == 0) & (ev["tfs_nlep"] == 0))

    def signal(self, ev):
        t = ev.true_muon()
        return ((ev["mc_incoming"] == 14) & (ev["mc_current"] == 1) & ev.true_vertex_in_fv(ZMIN, ZMAX, APOTHEM)
                & (t["theta"] < np.radians(20.0)) & (t["p"] > 1500.0) & (t["p"] < 20000.0)
                & self._topology(ev) & (self.wexp_true(ev) < W_MAX))

    def category(self, ev):
        npi0 = ev["tfs_npi0"]
        other_mesons = ev["tfs_npip"] + ev["tfs_npim"] + ev["tfs_nkpm"] + ev["tfs_nk0"] + ev["tfs_nomeson"]
        return np.select([(npi0 >= 1) & (other_mesons > 0), npi0 >= 1, other_mesons > 0], [0, 1, 2], 3)


def _fit(acc, use, cats, u):
    d = acc["x/fit/data"][use]
    fixed = acc["x/fit/sig"][use, u] + acc["x/fit/bkg/other"][use, u]
    for c in FIT_CATS:
        if c not in cats:
            fixed = fixed + acc[f"x/fit/bkg/{c}"][use, u]
    T = np.stack([acc[f"x/fit/bkg/{c}"][use, u] for c in cats], 1)
    wt = 1.0 / np.sqrt(np.maximum(d, 1.0))
    x, *_ = np.linalg.lstsq(T * wt[:, None], (d - fixed) * wt, rcond=None)
    pred = fixed + T @ x
    return np.maximum(x, 0.0), float((((d - pred) * wt) ** 2).sum()), int(use.sum())


def fit_regions():
    """(nominal, noWside, noMichel) masks over the global bins of the m_gg x region fit histogram."""
    B = Binning([FIT_MGG, FIT_REGION], ["mgg", "region"], flow=True)
    nx = len(FIT_MGG) + 1
    ix = np.arange(B.size) % nx
    iy = np.arange(B.size) // nx  # 1 + region for in-range bins
    nominal = B.inner_mask() & ~((ix == FIT_MGG.index(MGG_LO) + 1) & (iy == 1))
    return nominal, nominal & ((iy == 1) | (iy == 2)), nominal & ((iy == 1) | (iy == 3))


def background_scales(acc, meas):
    """Normalisations of the pi0_meson, pi0_other and charged backgrounds from a weighted least-squares fit
    (per universe) to the m_gg x region histogram, excluding the signal region (75 < m_gg < 195 MeV, no
    Michel, W_exp < 1.8 GeV) and the flow bins.  Signal and the 'other' category are fixed to the MC.
    The last two universes (BkgFit_*) use the alternative fits instead (systematic bands)."""
    nominal, no_w, no_michel = fit_regions()
    U = acc["x/fit/sig"].shape[1]
    coef = np.ones((len(FIT_CATS), U))
    for u in range(U):
        if u == U - 2:
            coef[:, u] = _fit(acc, no_w, FIT_CATS, u)[0]
        elif u == U - 1:
            x = _fit(acc, no_michel, ("pi0_meson", "charged"), u)[0]
            coef[:, u] = (x[0], 1.0, x[1])
        else:
            coef[:, u] = _fit(acc, nominal, FIT_CATS, u)[0]
    R = acc["data/sel"].shape[0]
    return {c: np.tile(coef[i], (R, 1)) for i, c in enumerate(FIT_CATS)}, coef
