"""MINERvA ME FHC numu CC1pi+ differential cross sections (tracked + untracked pions).

Reference: Granados, Messerly et al., arXiv:2605.24224 (code: github.com/MinervaExpt/CC-CH-pip-ana).
Variables: tpi (mixed tracked/untracked), thetapi (mixed, Tpi > 20 MeV signal), pmu, pzmu, ptmu,
thetamu, q2.

Signal (SignalDefinition::OnePi / OnePiThetaPi): numu CC, truth_is_fiducial, true vertex
5991.29 < z < 8400.91 mm, apothem 850 mm, theta_mu < 0.3491 rad, 1.5 < p_mu < 20 GeV/c,
W_exp(true) < 1.4 GeV, exactly one pi+ (Tpi_min < Tpi < 350 MeV), no other mesons, no photon
> 10 MeV, no e+-, any baryons.
Reco: vertex (same volume), MINOS-matched mu-, p_mu, theta_mu, < 2 isolated prongs, and either
  tracked : exactly one hadron with a quality endpoint Michel, LLR > 0 and node cuts, W_exp < 1.4
  untracked: good vertex/cluster Michel (LowRecoilPion), exactly one Michel candidate,
             Tpi_min < Tpi(range) < 350 MeV, pT_mu < 1.8 GeV/c, W_exp(Eavail) <= 1.4
The overlap uses the tracked reconstruction.  Background: W_exp sideband (reco W >= 1.5 GeV,
otherwise selected) fit of the true-W 1.4-1.8 and > 1.8 GeV background normalisations.

Not reproducible from public code: the MnvTune v4.7.1 (Tpi, Q2) charged-pion weight and the
"untracked pion" weight (GetChargedPionTuneWeight / GetUntrackedPionWeight are not public),
so the model is the v4.3.1 composition used by this analysis.
"""
import numpy as np

from xsec.binning import Binning
from xsec.kinematics import M_MU
from xsec.measurement import Measurement
from xsec.universes import standard_universes

M_PI_CC = 139.569  # CCNuPionIncConsts::CHARGED_PION_MASS
M_PI_MNV = 139.57018  # MinervaUnits::M_pion
M_HAD = (938.3 + 940.6) / 2  # CalcWexp nucleon mass
ZMIN, ZMAX, APOTHEM = 5991.29, 8400.91, 850.0
NSLOT = 5
W_MAX, W_SIDEBAND = 1400.0, 1500.0
WFIT_BINS = [0., 1000., 1100., 1200., 1300., 1400., 1500., 1750., 2000., 2500., 3000.]

# published binning (reference tables) and units; value -> table unit conversion of our MeV/rad
VARIABLES = {
    "tpi": dict(edges=[0, 20, 35, 50, 65, 80, 100, 125, 165, 200, 350], n_iter=10, table="tpi", cov="mixtpi"),
    "thetapi": dict(edges=[0, 15, 30, 45, 60, 76, 108, 122, 136, 150, 165, 180], n_iter=10, table="thetapi",
                    cov=None),
    "pmu": dict(edges=[1, 2, 3, 4, 5.5, 7.5, 10, 13, 20], n_iter=4, table="pmu", cov="pmu"),
    "pzmu": dict(edges=[0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20], n_iter=4, table="pzmu", cov="pzmu"),
    "ptmu": dict(edges=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5], n_iter=10, table="ptmu", cov="ptmu"),
    "thetamu": dict(edges=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20], n_iter=4, table="thetamu",
                    cov="thetamu_deg"),
    "q2": dict(edges=[0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.3, 2.0, 3.0], n_iter=10, table="q2",
               cov="q2"),
}


def calc_q2(enu, emu, theta):
    q2 = 2.0 * enu * (emu - np.sqrt(np.maximum(emu * emu - 105.659 ** 2, 0.0)) * np.cos(theta)) - 105.659 ** 2
    return np.maximum(q2, 0.0)


def calc_wexp(q2, ehad):
    w2 = M_HAD * M_HAD - q2 + 2.0 * M_HAD * ehad
    return np.where(w2 > 0, np.sqrt(np.maximum(w2, 0.0)), 0.0)


def fix_angle(a):
    a = np.abs(a)
    return np.where(a > np.pi, 2 * np.pi - a, a)


class CCPip(Measurement):
    tune_name = "v4.7.1"
    tune_options = dict(lowq2pi_channel="MENU1PI", coh_use_epi=True, diffractive_scale=1.4368, fsi_elastic=True,
                        fsi_absorption=True, geant_neutron_cv=True, geant_neutron_cv_reco_only=True)
    samples = ("sel", "wside")
    signal_sample = "sel"
    bkg_categories = ("wlow", "wmid", "whigh")
    extra_reco_columns = (["n_nonvtx_iso_blobs_all", "had_number", "had_recoil_def", "had_recoil_2trk",
                           "lrp_good", "lrp_ncand", "lrp_tpi", "lrp_theta", "lrp_eavail", "lrp_trk_hadmask", "nu_helicity"]
                          + [f"{c}_{i}" for c in ("had_isExiting", "had_isForked", "had_isODMatch", "had_isSideECAL",
                                                  "had_isTracker", "had_LLR", "had_pi_E_rcorr", "pi_E", "pi_Px",
                                                  "pi_Py", "pi_Pz", "pi_theta", "pi_Q0", "pi_Q1", "pi_Q2", "pi_Q3",
                                                  "pi_Q4", "pi_Q5") for i in range(NSLOT)])
    extra_truth_columns = ["t_is_fiducial", "tfs_npip", "tfs_npim", "tfs_npi0", "tfs_nkpm", "tfs_nk0", "tfs_nomeson",
                           "tfs_maxgammaE", "tfs_nlep", "tfs_pi_pdg_0", "tfs_pi_T_0", "tfs_pi_px_0", "tfs_pi_py_0",
                           "tfs_pi_pz_0", "tfs_pi_pdg_1", "tfs_pi_T_1", "w_fsi_absorption"]

    def __init__(self, variable="tpi"):
        self.variable = variable
        v = VARIABLES[variable]
        self.name = f"ccpip_me_fhc_{variable}"
        self.n_iter = v["n_iter"]
        self.tpi_min = 20.0 if variable == "thetapi" else 0.0
        self.reco_binning = Binning([v["edges"]], [variable])
        self.true_binning = Binning([v["edges"]], [variable])

    def universes(self, tune):
        return standard_universes(tune)

    # ------------------------------------------------------------ reco
    def preselect(self, ev):
        return (ev["minos_match"] == 1) & ev.reco_vertex_in_fv(ZMIN - 1, ZMAX + 1, APOTHEM)

    def _common(self, ev, mu):
        z = ev["vtx_2"]
        return ((z > ZMIN) & (z < ZMAX) & ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1)
                & (ev["minos_qp"] < 0) & (mu["p"] > 1500.0) & (mu["p"] < 20000.0) & (mu["theta"] < 0.3491)
                & (ev["n_nonvtx_iso_blobs_all"] < 2))

    def _tracked(self, ev, mu):
        """Tracked pion candidate (index or -1), tracked W_exp, Ehad."""
        n = len(ev)
        ns = sum(ev.has(f"had_isForked_{i}") for i in range(NSLOT))  # DFs may store fewer slots
        nhad = ev["had_number"]
        mask = ev["lrp_trk_hadmask"].astype(np.int64)
        cand = np.zeros((n, ns), dtype=bool)
        anyq = np.zeros(n, dtype=bool)
        for i in range(ns):
            q = ((ev[f"had_isForked_{i}"] == 0) & (ev[f"had_isExiting_{i}"] == 0) & (ev[f"had_isSideECAL_{i}"] == 0)
                 & (ev[f"had_isODMatch_{i}"] == 0) & (ev[f"had_isTracker_{i}"] == 1) & (i < nhad))
            anyq |= q
            e01 = ev[f"pi_Q0_{i}"] + ev[f"pi_Q1_{i}"]
            node = ((6 < e01) & (e01 < 32) & (2 < ev[f"pi_Q2_{i}"]) & (ev[f"pi_Q2_{i}"] < 22) & (0 < ev[f"pi_Q3_{i}"])
                    & (ev[f"pi_Q3_{i}"] < 19) & (0 < ev[f"pi_Q4_{i}"]) & (ev[f"pi_Q4_{i}"] < 31) & (0 < ev[f"pi_Q5_{i}"])
                    & (ev[f"pi_Q5_{i}"] < 60))
            cand[:, i] = ((mask >> i) & 1).astype(bool) & (ev[f"had_LLR_{i}"] > 0) & node & (i < nhad)
        ncand = cand.sum(1)
        passes = anyq & (ncand == 1)
        idx = np.where(passes, cand.argmax(1), -1)
        rows = np.arange(n)
        pick = lambda col: np.where(idx >= 0, np.stack([ev[f"{col}_{i}"] for i in range(ns)], 1)[rows, np.maximum(idx, 0)], np.nan)
        epi = pick("pi_E")
        ecal_pi = pick("had_pi_E_rcorr")
        ecal_nopi = ev["had_recoil_def"] - np.where(idx >= 0, ecal_pi, 0.0)
        # |t| of the pion candidate (GetCalRecoilEnergyNoPi_Corrected)
        phi = np.arctan2(np.tan(mu["thetaY"]), np.tan(mu["thetaX"]))
        pxmu, pymu, pzmu = mu["pt"] * np.cos(phi), mu["pt"] * np.sin(phi), mu["pz"]
        t = ((epi + mu["E"] - pzmu - pick("pi_Pz")) ** 2 + (pick("pi_Px") + pxmu) ** 2 + (pick("pi_Py") + pymu) ** 2)
        # original quirk: bins in MeV multiplied again by 1e3 -> +60 MeV for 0 < ecal_nopi < 1000 when |t| >= 125e3
        corr = np.where((idx >= 0) & ~(t < 125e3) & (ecal_nopi > 0) & (ecal_nopi < 1000), ecal_nopi + 60.0, ecal_nopi)
        ecal = np.where(ecal_nopi > 1000, ev["had_recoil_2trk"], corr)
        etrk = np.where((idx >= 0) & ~(ecal_nopi > 1000), epi, 0.0)
        ehad = ecal + etrk
        w = calc_wexp(calc_q2(mu["E"] + ehad, mu["E"], mu["theta"]), ehad)
        tpi = epi - M_PI_CC
        theta = np.degrees(fix_angle(pick("pi_theta")))
        return passes, w, ehad, tpi, theta

    def _untracked(self, ev, mu):
        ehad = ev["lrp_eavail"] + ev["lrp_ncand"] * M_PI_MNV
        w = calc_wexp(calc_q2(mu["E"] + ehad, mu["E"], mu["theta"]), ehad)
        tpi = ev["lrp_tpi"]
        passes = ((ev["lrp_good"] == 1) & (ev["lrp_ncand"] == 1) & (tpi > self.tpi_min) & (tpi < 350.0)
                  & (mu["pt"] < 1800.0) & (w > 0.0))
        return passes, w, ehad, tpi, np.degrees(ev["lrp_theta"])

    def _reco_all(self, ev, shift):
        mu = ev.muon(shift)
        com = self._common(ev, mu)
        tp, tw, teh, ttpi, tth = self._tracked(ev, mu)
        up, uw, ueh, utpi, uth = self._untracked(ev, mu)
        tr_except_w, un_except_w = com & tp, com & up
        tr_sel = tr_except_w & (tw < W_MAX)
        un_sel = un_except_w & ~(uw > W_MAX)
        sel = tr_sel | un_sel
        tr_side = tr_except_w & (tw >= W_SIDEBAND)
        un_side = un_except_w & (uw > W_MAX) & (uw >= W_SIDEBAND)
        side = (tr_side | un_side) & ~sel
        # GetMixed*: tracked cuts > untracked cuts > tracked sideband > untracked sideband > tracked except W
        use_tr = np.select([tr_sel, un_sel, tr_side, un_side, tr_except_w], [True, False, True, False, True], False)
        ehad = np.where(use_tr, teh, ueh)
        wmix = np.where(use_tr, tw, uw)
        tpi = np.where(use_tr, ttpi, utpi)
        thpi = np.where(use_tr, tth, uth)
        q2 = calc_q2(mu["E"] + ehad, mu["E"], mu["theta"])
        vals = dict(tpi=tpi, thetapi=thpi, pmu=mu["p"] / 1e3, pzmu=mu["pz"] / 1e3, ptmu=mu["pt"] / 1e3,
                    thetamu=np.degrees(mu["theta"]), q2=q2 / 1e6)
        return vals, sel, side, tr_except_w | un_except_w, wmix

    def reco(self, ev, shift):
        vals, sel, side, _, _ = self._reco_all(ev, shift)
        return (vals[self.variable],), {"sel": sel, "wside": side}

    def extra(self, ev, shift):
        _, _, _, except_w, wmix = self._reco_all(ev, shift)
        return {"wfit": (Binning([WFIT_BINS], ["w"], flow=True), (wmix,), except_w)}

    # ------------------------------------------------------------ truth
    def _true_pi(self, ev):
        is_pip = ev["tfs_pi_pdg_0"] == 211
        from xsec.kinematics import rotate_to_beam
        px, py, pz = rotate_to_beam(ev["tfs_pi_px_0"], ev["tfs_pi_py_0"], ev["tfs_pi_pz_0"])
        p = np.sqrt(px * px + py * py + pz * pz)
        th = np.degrees(np.arccos(np.clip(np.divide(pz, p, out=np.ones_like(p), where=p > 0), -1, 1)))
        return np.where(is_pip, ev["tfs_pi_T_0"], np.nan), np.where(is_pip, th, np.nan)

    def wexp_true(self, ev):
        ehad = ev["mc_Enu"] - ev["mc_lep_3"]
        return calc_wexp(ev["mc_Q2"], ehad)

    def truth(self, ev):
        t = ev.true_muon()
        tpi, thpi = self._true_pi(ev)
        vals = dict(tpi=tpi, thetapi=thpi, pmu=t["p"] / 1e3, pzmu=t["pz"] / 1e3, ptmu=t["pt"] / 1e3,
                    thetamu=np.degrees(t["theta"]), q2=ev["mc_Q2"] / 1e6)
        return (vals[self.variable],)

    def signal(self, ev):
        t = ev.true_muon()
        tpi, _ = self._true_pi(ev)
        z = ev["mc_vtx_2"]
        w = self.wexp_true(ev)
        topo = ((ev["tfs_npip"] == 1) & (ev["tfs_npim"] == 0) & (ev["tfs_npi0"] == 0) & (ev["tfs_nkpm"] + ev["tfs_nk0"] == 0)
                & (ev["tfs_nomeson"] == 0) & (ev["tfs_maxgammaE"] <= 10.0) & (ev["tfs_nlep"] == 0))
        return ((ev["mc_incoming"] == 14) & (ev["mc_current"] == 1) & (ev["t_is_fiducial"] == 1) & (z > ZMIN) & (z < ZMAX)
                & ev.true_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (t["theta"] < 0.3491) & (w > 0.0) & (w < W_MAX)
                & topo & (tpi > self.tpi_min) & (tpi < 350.0) & (t["p"] > 1500.0) & (t["p"] < 20000.0))

    def category(self, ev):
        w = self.wexp_true(ev)
        return np.where(w > 1800.0, 2, np.where(w > 1400.0, 1, 0))


def background_scales(acc, meas, wmin=W_SIDEBAND):
    """Fit the wmid and whigh normalisations to the reco W_exp distribution of events passing all cuts
    except W, over W >= 1.5 GeV (TFractionFitter in the original; weighted least squares here)."""
    B = Binning([WFIT_BINS], ["w"], flow=True)
    lo = np.searchsorted(B.edges[0], wmin)  # flow index of the first fit bin
    sel = np.zeros(B.size, dtype=bool)
    sel[lo + 1:] = True  # in-range bins above wmin and the overflow
    d = acc["x/wfit/data"][sel]
    fixed = acc["x/wfit/sig"][sel] + acc["x/wfit/bkg/wlow"][sel]
    tm, th = acc["x/wfit/bkg/wmid"][sel], acc["x/wfit/bkg/whigh"][sel]
    U = fixed.shape[1]
    sig2 = np.maximum(d, 1.0)[:, None]
    scales = np.ones((2, U))
    for u in range(U):
        A = np.stack([tm[:, u], th[:, u]], 1) / np.sqrt(sig2)
        b = (d - fixed[:, u]) / np.sqrt(sig2[:, 0])
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        scales[:, u] = np.maximum(x, 0.0)
    R = acc["data/sel"].shape[0]
    return {"wmid": np.tile(scales[0], (R, 1)), "whigh": np.tile(scales[1], (R, 1))}, scales
