"""MINERvA ME FHC numu CC0pi Np (>=1 proton) single-transverse kinematic imbalance (STKI).

References: X.-G. Lu et al. (MINERvA), PRL 121, 022504 (2018), arXiv:1805.05486 (LE signal definition,
pn with b = 27.13 MeV); T. Cai et al. (MINERvA), PRD 101, 092001 (2020), arXiv:1910.08658;
X.-G. Lu et al., PRC 94, 015503 (2016), arXiv:1512.05748 (STKI definitions).

Signal (generator final state after FSI): numu CC, true vertex in the tracker FV; mu- with
1.5 < p_mu < 20 GeV/c and theta_mu < 20 deg; no mesons (pi+-, pi0, K, other), no heavy baryons, no
leptons other than the muon, no photon above 10 MeV (the ccqelike_3d definition); at least one proton
with 0.45 < p_p < 1.2 GeV/c and theta_p < 70 deg.  The STKI proton is the most energetic proton in
that phase space (searched among the two leading final-state protons, the generator-summary slots).
Neutrons and additional protons are allowed.

Reco: MINOS-matched mu- (helicity), reco vertex in FV, theta_mu < 20 deg, 1.5 < p_mu < 20 GeV/c,
<= 1 dead discriminator; >= 2 tracks, all non-muon tracks are proton candidates (multiplicity ==
2 + n_secp) with a Q2_QE-dependent proton score1 cut (CCQENu bands) and score1 > 0.5 for the primary
candidate (ESC-like quality, cf. Lu et al.); the primary proton candidate
(MasterAnaDev_proton_*, momentum from the dE/dx-range fit) in 0.45 < p_p < 1.2 GeV/c, theta_p < 70 deg;
no Michel electron; unattached recoil (non-vertex isolated blobs + dispersed energy) < 150 MeV.
Background: pion-production template scales from a fit of the (recoil, Michel) distribution of
events passing all other cuts (signal, 1 charged pion, >= 1 pi0 or multi-pion floated).
"""
import numpy as np

from xsec.binning import Binning
from xsec.kinematics import M_MU, M_N, M_P, Shift, q2_qe, recoil_shift, rotate_to_beam
from xsec.measurement import Measurement
from xsec.universes import Universe, response_universes, standard_universes

ZMIN, ZMAX, APOTHEM = 5980.0, 8422.0, 850.0
PMU_MIN, PMU_MAX, THMU_MAX = 1500.0, 20000.0, np.radians(20.0)
PP_MIN, PP_MAX, THP_MAX = 450.0, 1200.0, 70.0  # MeV/c, deg
RECOIL_MAX = 150.0  # MeV
PRIMARY_SCORE1_MIN = 0.5  # clean (Bragg-peak) primary proton: suppresses reinteracting/misreconstructed protons
NSECP_SLOTS = 6
NSECP_P_SLOTS = 4
PSCORE_Q2 = [0.2, 0.6]
PSCORE_MIN = [0.2, 0.1, 0.0]
# pn (Lu et al. 2018): m_A = 12C mass, m_A' = m_A - m_n + b, b = 27.13 MeV
M_A = 6 * M_N + 6 * M_P - 92.16
M_AP = M_A - M_N + 27.13
RFIT_RECOIL = [0, 50, 100, 150, 200, 300, 400, 600, 800, 1000, 1500]  # MeV
PROTON_VARIANTS = (("Proton_Energy_BetheBloch", ("bb_dn", "bb_up")), ("Proton_Energy_MEU", ("meu_dn", "meu_up")),
                   ("Proton_Energy_Mass", ("mass_dn", "mass_up")), ("Proton_Energy_Birks", ("birks",)))

DPT = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.475, 0.55, 0.65, 0.8, 1.0]
VARIABLES = {
    "dpt": dict(edges=[DPT], names=["dpt"], n_iter=6),
    "dalphat": dict(edges=[[0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]], names=["dalphat"], n_iter=6),
    "pn": dict(edges=[[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75, 1.0]], names=["pn"], n_iter=6),
    "dphit": dict(edges=[[0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 130, 180]], names=["dphit"], n_iter=6),
    "pp": dict(edges=[[0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.2]], names=["pp"], n_iter=4),
    "thetap": dict(edges=[[0, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]], names=["thetap"], n_iter=4),
    "nprot": dict(edges=[[0.5, 1.5, 2.5]], names=["nprot"], n_iter=4),
    "dpt_dalphat": dict(edges=[[0, 0.1, 0.2, 0.3, 0.45, 1.0], [0, 45, 90, 135, 180]], names=["dpt", "dalphat"],
                        n_iter=6),
}


def proton_score_min(q2qe_gev2):
    return np.asarray(PSCORE_MIN)[np.searchsorted(PSCORE_Q2, q2qe_gev2, side="right")]


def stki(mpx, mpy, mpz, emu, ppx, ppy, ppz):
    """STKI from beam-frame muon and proton momenta [MeV]; returns dict in GeV/c and deg."""
    with np.errstate(invalid="ignore", divide="ignore"):
        dpx, dpy = mpx + ppx, mpy + ppy
        dpt = np.hypot(dpx, dpy)
        ptm, ptp = np.hypot(mpx, mpy), np.hypot(ppx, ppy)
        dalphat = np.degrees(np.arccos(np.clip(-(mpx * dpx + mpy * dpy) / (ptm * dpt), -1, 1)))
        dphit = np.degrees(np.arccos(np.clip(-(mpx * ppx + mpy * ppy) / (ptm * ptp), -1, 1)))
        pp = np.sqrt(ppx * ppx + ppy * ppy + ppz * ppz)
        ep = np.sqrt(pp * pp + M_P * M_P)
        R = M_A + mpz + ppz - emu - ep
        dpl = 0.5 * R - (M_AP * M_AP + dpt * dpt) / (2 * R)
        pn = np.hypot(dpt, dpl)
        thp = np.degrees(np.arccos(np.clip(ppz / pp, -1, 1)))
    return dict(dpt=dpt / 1e3, dalphat=dalphat, dphit=dphit, pn=pn / 1e3, pp=pp / 1e3, thetap=thp)




class STKI(Measurement):
    tune_name = "v4.4.1"
    tune_options = dict(fsi_elastic=True)
    fv = (ZMIN, ZMAX, APOTHEM)
    samples = ("sel",)
    signal_sample = "sel"
    bkg_categories = ("1chpi", "1pi0", "multipi", "other")
    extra_reco_columns = (["nu_helicity", "multiplicity", "n_improved_michel", "nonvtx_iso_blobs_E", "dispr_id_E",
                           "n_secp", "p_score1", "p_P", "p_Px", "p_Py", "p_Pz", "p_T"]
                          + [f"p_E_{v}" for v in ("bb_up", "bb_dn", "meu_up", "meu_dn", "mass_up", "mass_dn", "birks")]
                          + [f"p_s1_{v}" for v in ("bb_up", "bb_dn", "meu_up", "meu_dn", "mass_up", "mass_dn", "birks")]
                          + [f"secp_score1_{i}" for i in range(NSECP_SLOTS)]
                          + [f"secp_{c}_{i}" for c in ("P", "theta") for i in range(NSECP_P_SLOTS)])
    extra_mc_reco_columns = [f"rerr_dispr_{p}" for p in ("proton", "meson", "em", "other")]
    extra_truth_columns = (["tfs_np", "tfs_npip", "tfs_npim", "tfs_npi0", "tfs_nkpm", "tfs_nk0", "tfs_nomeson",
                            "tfs_nhbaryon", "tfs_maxgammaE", "tfs_nlep", "w_fsi_elastic"]
                           + [f"tfs_p_{c}_{k}" for c in ("T", "px", "py", "pz") for k in (0, 1)])

    def __init__(self, variable="dpt"):
        self.variable = variable
        v = VARIABLES[variable]
        self.name = f"cc0pinp_stki_me_fhc_{variable}"
        self.n_iter = v["n_iter"]
        self.reco_binning = Binning(v["edges"], v["names"])
        self.true_binning = Binning(v["edges"], v["names"])

    def universes(self, tune):
        U = standard_universes(tune) + response_universes()
        for band, variants in PROTON_VARIANTS:
            for j, var in enumerate(variants):
                U.append(Universe(band, j, "lat", shift=Shift(f"p_{var}", proton=var)))
        for j, var in enumerate(("bb_dn", "bb_up")):
            U.append(Universe("Proton_Score_BetheBloch", j, "lat", shift=Shift(f"ps_{var}", extra=(("pscore", var),))))
        for j, var in enumerate(("meu_dn", "meu_up")):
            U.append(Universe("Proton_Score_MEU", j, "lat", shift=Shift(f"ps_{var}", extra=(("pscore", var),))))
        for j, var in enumerate(("mass_dn", "mass_up")):
            U.append(Universe("Proton_Score_Mass", j, "lat", shift=Shift(f"ps_{var}", extra=(("pscore", var),))))
        U.append(Universe("Proton_Score_Birks", 0, "lat", shift=Shift("ps_birks", extra=(("pscore", "birks"),))))
        return U

    # ------------------------------------------------------------ reco
    def preselect(self, ev):
        return (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1) & (ev["has_vtx"] == 1)
                & (ev["multiplicity"] >= 2) & (ev["p_P"] > 0))

    def _reco_all(self, ev, shift):
        mu = ev.muon(shift)
        phi = np.arctan2(np.tan(mu["thetaY"]), np.tan(mu["thetaX"]))
        mpx, mpy, mpz = mu["pt"] * np.cos(phi), mu["pt"] * np.sin(phi), mu["pz"]
        # proton (detector frame -> beam frame), energy-scale variants shift the kinetic energy
        T = ev["p_T"].astype(np.float64)
        P = ev["p_P"].astype(np.float64)
        if shift.proton and ev.is_mc:
            d = ev[f"p_E_{shift.proton}"].astype(np.float64)
            d = np.where(np.abs(d) < 1000.0, d, 0.0)
            Tn = np.maximum(T + d, 1.0)
            Pn = np.sqrt(Tn * Tn + 2 * M_P * Tn)
            fac = np.divide(Pn, P, out=np.ones_like(P), where=P > 0)
        else:
            Pn, fac = P, 1.0
        ppx, ppy, ppz = rotate_to_beam(ev["p_Px"] * fac, ev["p_Py"] * fac, ev["p_Pz"] * fac)
        k = stki(mpx, mpy, mpz, mu["E"], ppx, ppy, ppz)
        # proton identification
        q2 = q2_qe(mu["E"], mu["p"], mu["theta"]) / 1e6
        q2 = np.where(np.isfinite(q2), q2, 0.0)
        smin = proton_score_min(np.maximum(q2, 0.0))
        s1 = ev["p_score1"].astype(np.float64)
        pv = shift.get("pscore", "")
        if pv and ev.is_mc:
            ds = ev[f"p_s1_{pv}"].astype(np.float64)
            s1 = s1 + np.where(np.abs(ds) < 1.0, ds, 0.0)
        nsec = ev["n_secp"]
        sec_ok = np.ones(len(ev), dtype=bool)
        for i in range(NSECP_SLOTS):
            sec_ok &= ~((i < nsec) & (ev[f"secp_score1_{i}"] < smin))
        pid = (s1 >= np.maximum(smin, PRIMARY_SCORE1_MIN)) & sec_ok & (ev["multiplicity"] == 2 + nsec)
        base = (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["has_vtx"] == 1) & (ev["minos_match"] == 1)
                & (ev["nu_helicity"] == 1) & (mu["theta"] < THMU_MAX) & (mu["p"] > PMU_MIN) & (mu["p"] < PMU_MAX)
                & (ev["dead"] <= 1) & (ev["multiplicity"] >= 2) & (P > 0) & pid
                & (k["pp"] * 1e3 > PP_MIN) & (k["pp"] * 1e3 < PP_MAX) & (k["thetap"] < THP_MAX))
        recoil = ev["nonvtx_iso_blobs_E"] + ev["dispr_id_E"] + recoil_shift(ev, shift, "rerr_dispr_")
        michel = ev["n_improved_michel"] > 0
        # proton multiplicity in the STKI phase space: primary + secondary candidates (4 momentum slots)
        nsec_ps = np.zeros(len(ev))
        for i in range(NSECP_P_SLOTS):
            p2, th2 = ev[f"secp_P_{i}"], np.degrees(ev[f"secp_theta_{i}"])
            nsec_ps += ((i < nsec) & (p2 > PP_MIN) & (p2 < PP_MAX) & (th2 < THP_MAX)).astype(float)
        k["nprot"] = 1.0 + np.minimum(nsec_ps, 1.0)
        return k, base, recoil, michel

    def reco(self, ev, shift):
        k, base, recoil, michel = self._reco_all(ev, shift)
        sel = base & ~michel & (recoil < RECOIL_MAX)
        return tuple(k[n] for n in self.reco_binning.names), {"sel": sel}

    def extra(self, ev, shift):
        _, base, recoil, michel = self._reco_all(ev, shift)
        B = Binning([RFIT_RECOIL, [-0.5, 0.5, 1.5]], ["recoil", "michel"])
        return {"rfit": (B, (recoil, michel.astype(float)), base)}

    # ------------------------------------------------------------ truth
    def _true(self, ev):
        lx, ly, lz = rotate_to_beam(ev["mc_lep_0"].astype(np.float64), ev["mc_lep_1"].astype(np.float64),
                                    ev["mc_lep_2"].astype(np.float64))
        emu = ev["mc_lep_3"].astype(np.float64)
        out, ok, nps = None, np.zeros(len(ev), dtype=bool), np.zeros(len(ev))
        for slot in (0, 1):
            px, py, pz = rotate_to_beam(ev[f"tfs_p_px_{slot}"].astype(np.float64),
                                        ev[f"tfs_p_py_{slot}"].astype(np.float64),
                                        ev[f"tfs_p_pz_{slot}"].astype(np.float64))
            k = stki(lx, ly, lz, emu, px, py, pz)
            inps = (np.nan_to_num(k["pp"]) * 1e3 > PP_MIN) & (np.nan_to_num(k["pp"]) * 1e3 < PP_MAX) & (
                np.nan_to_num(k["thetap"], nan=180.0) < THP_MAX)
            use = inps & ~ok
            if out is None:
                out = {n: np.where(use, v, np.nan) for n, v in k.items()}
            else:
                for n, v in k.items():
                    out[n] = np.where(use, v, out[n])
            ok |= inps
            nps += inps
        out["nprot"] = np.where(ok, np.minimum(nps, 2.0), np.nan)
        return out, ok

    def truth(self, ev):
        k, _ = self._true(ev)
        return tuple(k[n] for n in self.true_binning.names)

    def meson_free(self, ev):
        return ((ev["tfs_npip"] + ev["tfs_npim"] + ev["tfs_npi0"] + ev["tfs_nkpm"] + ev["tfs_nk0"]
                 + ev["tfs_nomeson"] + ev["tfs_nhbaryon"] + ev["tfs_nlep"]) == 0) & (ev["tfs_maxgammaE"] < 10.0)

    def signal(self, ev):
        t = ev.true_muon()
        _, okp = self._true(ev)
        return ((ev["mc_incoming"] == 14) & (ev["mc_current"] == 1) & self.meson_free(ev)
                & ev.true_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (t["theta"] < THMU_MAX) & (t["p"] > PMU_MIN)
                & (t["p"] < PMU_MAX) & okp)

    def category(self, ev):
        npic = ev["tfs_npip"] + ev["tfs_npim"]
        npi0 = ev["tfs_npi0"]
        cc = (ev["mc_incoming"] == 14) & (ev["mc_current"] == 1)
        cat = np.full(len(ev), 3)
        cat[cc & (npic == 1) & (npi0 == 0)] = 0
        cat[cc & (npic == 0) & (npi0 == 1)] = 1
        cat[cc & (npic + npi0 >= 2)] = 2
        return cat


FIT_GROUPS = (("signal",), ("1chpi",), ("1pi0", "multipi"))


def background_scales(acc, meas):
    """Weighted least-squares fit (per universe) of the (recoil, Michel) distribution of events that pass
    every cut except the Michel veto and the recoil cut: free normalisations for the signal, the
    1-charged-pion background and the pi0/multi-pion background; 'other' fixed to the MC."""
    B = Binning([RFIT_RECOIL, [-0.5, 0.5, 1.5]], ["recoil", "michel"])
    sx = B.shape[0]
    g = np.arange(B.size)
    ix, iy = g % sx, g // sx
    use = (ix >= 1) & (iy >= 1) & (iy <= 2)  # recoil in range or overflow, both Michel bins
    d = acc["x/rfit/data"][use]
    fixed = acc["x/rfit/bkg/other"][use]
    temps = []
    for grp in FIT_GROUPS:
        t = 0.0
        for c in grp:
            t = t + (acc["x/rfit/sig"] if c == "signal" else acc[f"x/rfit/bkg/{c}"])[use]
        temps.append(t)
    U = fixed.shape[1]
    w = 1.0 / np.sqrt(np.maximum(d, 1.0))
    coef = np.ones((len(FIT_GROUPS), U))
    for u in range(U):
        A = np.stack([t[:, u] for t in temps], 1) * w[:, None]
        b = (d - fixed[:, u]) * w
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        coef[:, u] = np.maximum(x, 0.0)
    R = acc["data/sel"].shape[0]
    out = {}
    for j, grp in enumerate(FIT_GROUPS):
        for c in grp:
            if c != "signal":
                out[c] = np.tile(coef[j], (R, 1))
    return out, coef
