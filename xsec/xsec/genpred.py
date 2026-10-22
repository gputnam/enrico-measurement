"""Generator (nuflow) predictions in the signal definition and binning of a measurement.

A nuflow ROOT file (`events` + `run` trees, see model/README.md) is converted into a pandas
dataframe with the same *truth* columns the makedf DFs carry (mc_*, tfs_*, t_*), in the same
units (MeV, mm, detector frame).  A measurement's own `signal(ev)` and `truth(ev)` are then
evaluated on it, so the generator prediction uses exactly the signal definition and binning of the
measurement.

Normalisation (unweighted or weighted events, per nucleon of the CH target):
    dsigma/dX [cm^2/nucleon/unit] = flux_avg_xsec * 1e-38 / nucleons_per_target
                                    * sum_{signal in bin} w / sum_all w / bin_volume
with flux_avg_xsec the generator's flux-averaged total cross section per target unit over its
flux, which must be the measurement flux (export_flux.py) over 0-100 GeV.

Columns provided (anything else a measurement's truth/signal needs must be added here):
    mc_incoming, mc_current, mc_intType, mc_targetZ, mc_targetA, mc_targetNucleon, mc_Enu, mc_Q2, mc_W,
    mc_x, mc_y, mc_nu_0..3, mc_vtx_0..3 (a fixed point inside every tracker FV), mc_lep_0..3,
    mc_resID (-1), mc_charm (0), t_is_fiducial (1), t_N_pip/pim/pi0,
    tfs_* (TruthSummary of makedf/src/Derived.cxx, incl. tfs_pi0_* and tfs_gamma_*_0),
    t_theta_mu, t_p_mu, t_pt_mu, t_pz_mu, t_q0, t_q3.
Not available (MINERvA-simulation only, must not be used in signal()/truth()): t_pi0_* (tuple order),
t_N_gamma_pi0, t_gamma_* (pi0 decay photons); use tfs_pi0_* instead.
Weight columns (w_fsi_*, geant_*, cv_wgt_tuple) are 1.  GENIE knob columns (g_*) are absent: a
signal definition or truth variable must not use them.
"""
import numpy as np
import pandas as pd
import uproot

from .kinematics import BEAM_ANGLE, Event

M_P_DF, M_N_DF, M_PI_DF = 938.272, 939.565, 139.570  # MeV (makedf Derived.cxx kMp, kMn, kMpi)
M_PI0_DF = 134.977  # MeV (makedf Derived.cxx kMpi0)
K_SLOTS = 2  # leading protons / pions / pi0s stored per event (makedf TruthSummary::K)
VERTEX = (0.0, 0.0, 7000.0)  # mm: inside the tracker FV of every measurement

# NuHepMC signal process id -> GENIE mc_intType (1 QE, 2 RES, 3 DIS, 4 COH, 5 AMNUGAMMA, 7 NUEEL, 8 MEC/2p2h)
# NuHepMC ranges: 100 coherent, 200 CCQE, 250 NC elastic, 300 CC MEC, 350 NC MEC, 400 CC RES, 450 NC RES,
# 500 CC SIS/non-resonant pion, 550 NC SIS, 600 CC DIS, 650 NC DIS.
INTTYPE_BY_RANGE = [(100, 200, 4), (200, 300, 1), (300, 400, 8), (400, 500, 2), (500, 600, 3), (600, 700, 3)]


def int_type(process_id):
    out = np.zeros(len(process_id), dtype=int)
    for lo, hi, t in INTTYPE_BY_RANGE:
        out[(process_id >= lo) & (process_id < hi)] = t
    return out


def _from_beam(px, py, pz):
    """Inverse of kinematics.rotate_to_beam: beam frame (nu along +z) -> MINERvA detector frame."""
    c, s = np.cos(BEAM_ANGLE), np.sin(BEAM_ANGLE)
    return px, c * py + s * pz, -s * py + c * pz


def read_run(path):
    with uproot.open(path) as f:
        r = f["run"].arrays(library="np")
    return {k: v[0] for k, v in r.items()}


def load_events(path, max_events=None):
    """Read a nuflow ROOT file into a truth dataframe (makedf column conventions)."""
    with uproot.open(path) as f:
        t = f["events"]
        a = t.arrays(["weight", "channel", "nu_pdg", "target_pdg", "Enu", "Q2", "q0", "q3", "W", "hit_pdg",
                      "pdg", "status", "E", "px", "py", "pz"], entry_stop=max_events, library="ak")
    import awkward as ak
    n = len(a)
    fs = a["status"] == 1
    pdg = a["pdg"][fs]
    E = a["E"][fs] * 1e3
    px, py, pz = a["px"][fs] * 1e3, a["py"][fs] * 1e3, a["pz"][fs] * 1e3
    apdg = abs(pdg)
    nu = ak.to_numpy(a["nu_pdg"]).astype(int)

    # primary lepton: leading charged lepton of the neutrino flavour (CC) else the outgoing neutrino (NC)
    clep = abs(nu) - 1
    is_cl = apdg == ak.Array(clep)[:, np.newaxis]
    has_cl = ak.to_numpy(ak.any(is_cl, axis=1))
    is_nl = apdg == ak.Array(abs(nu))[:, np.newaxis]
    lep_mask = ak.where(ak.Array(has_cl)[:, np.newaxis], is_cl, is_nl)
    lep_E = ak.where(lep_mask, E, -1.0)
    li = ak.argmax(lep_E, axis=1, keepdims=True)
    lep = dict(E=ak.to_numpy(ak.fill_none(ak.firsts(E[li]), 0.0)),
               px=ak.to_numpy(ak.fill_none(ak.firsts(px[li]), 0.0)),
               py=ak.to_numpy(ak.fill_none(ak.firsts(py[li]), 0.0)),
               pz=ak.to_numpy(ak.fill_none(ak.firsts(pz[li]), 0.0)))
    # remove the primary lepton from the hadronic final state
    idx = ak.local_index(pdg, axis=1)
    keep = idx != ak.fill_none(ak.firsts(li), -1)
    hp, hE, hx, hy, hz = pdg[keep], E[keep], px[keep], py[keep], pz[keep]
    ha = abs(hp)

    def count(m):
        return ak.to_numpy(ak.sum(m, axis=1)).astype(int)

    def ssum(v, m):
        return ak.to_numpy(ak.sum(ak.where(m, v, 0.0), axis=1)).astype(float)

    is_p, is_n = hp == 2212, hp == 2112
    is_pic, is_pi0 = ha == 211, hp == 111
    is_kpm, is_k0 = ha == 321, (ha == 311) | (ha == 130) | (ha == 310)
    is_g = hp == 22
    is_lep = (ha >= 11) & (ha <= 16)
    is_frag = ha > 1000000000
    is_om = (ha >= 100) & (ha < 1000) & ~is_pic & ~is_pi0 & ~is_kpm & ~is_k0
    is_hb = (ha >= 1000) & (ha < 10000) & ~is_p & ~is_n
    is_other = ~(is_p | is_n | is_pic | is_pi0 | is_kpm | is_k0 | is_g | is_lep | is_frag | is_om | is_hb)
    Tp = np.maximum(hE - M_P_DF, 0.0)
    Tn = np.maximum(hE - M_N_DF, 0.0)
    Tpi = hE - M_PI_DF
    eavail = (ssum(Tp, is_p) + ssum(Tpi, is_pic) + ssum(hE, is_pi0 | is_kpm | is_k0 | is_g | is_om | is_other)
              + ssum(hE - M_P_DF, is_hb))

    d = {}
    d["mc_incoming"] = nu
    d["mc_current"] = np.where(has_cl, 1, 2)
    d["mc_intType"] = int_type(ak.to_numpy(a["channel"]).astype(int))
    tgt = ak.to_numpy(a["target_pdg"]).astype(np.int64)
    d["mc_targetZ"] = (tgt // 10000) % 1000
    d["mc_targetA"] = (tgt // 10) % 1000
    d["mc_targetNucleon"] = ak.to_numpy(a["hit_pdg"]).astype(int)
    d["mc_Enu"] = ak.to_numpy(a["Enu"]) * 1e3
    d["mc_Q2"] = ak.to_numpy(a["Q2"]) * 1e6
    d["mc_W"] = ak.to_numpy(a["W"]) * 1e3
    q0 = ak.to_numpy(a["q0"]) * 1e3
    d["mc_x"] = np.divide(d["mc_Q2"], 2 * M_N_DF * q0, out=np.zeros(n), where=q0 > 0)
    d["mc_y"] = q0 / d["mc_Enu"]
    nx, ny, nz = _from_beam(np.zeros(n), np.zeros(n), d["mc_Enu"])
    d["mc_nu_0"], d["mc_nu_1"], d["mc_nu_2"], d["mc_nu_3"] = nx, ny, nz, d["mc_Enu"]
    d["mc_resID"] = np.full(n, -1)
    d["mc_charm"] = np.zeros(n, dtype=int)
    for i, v in enumerate(VERTEX):
        d[f"mc_vtx_{i}"] = np.full(n, v)
    lx, ly, lz = _from_beam(lep["px"], lep["py"], lep["pz"])
    d["mc_lep_0"], d["mc_lep_1"], d["mc_lep_2"], d["mc_lep_3"] = lx, ly, lz, lep["E"]
    d["t_is_fiducial"] = np.ones(n, dtype=int)
    d["mc_vtx_3"] = np.zeros(n)
    d["tfs_np"], d["tfs_nn"] = count(is_p), count(is_n)
    d["tfs_npip"], d["tfs_npim"], d["tfs_npi0"] = count(hp == 211), count(hp == -211), count(is_pi0)
    d["tfs_nkpm"], d["tfs_nk0"], d["tfs_nomeson"] = count(is_kpm), count(is_k0), count(is_om)
    d["tfs_nhbaryon"], d["tfs_ngamma"], d["tfs_nlep"] = count(is_hb), count(is_g), count(is_lep)
    d["tfs_np_T120"] = count(is_p & (Tp > 120.0))
    d["tfs_nnucfrag"] = count(is_frag)
    d["tfs_sumTp"], d["tfs_sumTn"] = ssum(Tp, is_p), ssum(Tn, is_n)
    d["tfs_maxgammaE"] = ak.to_numpy(ak.fill_none(ak.max(ak.where(is_g, hE, 0.0), axis=1), 0.0)).astype(float)
    d["tfs_eavail"] = eavail
    nuc = is_p | is_n
    nE = ak.where(nuc, hE, 0.0)
    ni = ak.argmax(nE, axis=1, keepdims=True)
    lead_E = ak.to_numpy(ak.fill_none(ak.firsts(nE[ni]), 0.0))
    d["tfs_leadnuc_E"] = lead_E
    d["tfs_leadnuc_pdg"] = np.where(lead_E > 0, ak.to_numpy(ak.fill_none(ak.firsts(hp[ni]), 0)), 0)
    # leading-K protons / charged pions by kinetic energy (detector frame momenta)
    for name, m, T in (("p", is_p, Tp), ("pi", is_pic, Tpi)):
        Ts = ak.where(m, T, -np.inf)
        order = ak.argsort(Ts, axis=1, ascending=False)
        for k in range(K_SLOTS):
            sel = order[:, k:k + 1]
            valid = ak.to_numpy(ak.fill_none(ak.firsts(ak.where(m, True, False)[sel]), False))
            g = lambda v: np.where(valid, ak.to_numpy(ak.fill_none(ak.firsts(v[sel]), np.nan)), np.nan)
            bx, by, bz = _from_beam(g(hx), g(hy), g(hz))
            d[f"tfs_{name}_T_{k}"] = g(T)
            d[f"tfs_{name}_px_{k}"], d[f"tfs_{name}_py_{k}"], d[f"tfs_{name}_pz_{k}"] = bx, by, bz
            if name == "pi":
                d[f"tfs_pi_pdg_{k}"] = np.where(valid, ak.to_numpy(ak.fill_none(ak.firsts(hp[sel]), 0)), 0)
    # leading-K pi0s by kinetic energy (makedf TruthSummary tfs_pi0_*), detector frame
    Tpi0 = hE - M_PI0_DF
    order = ak.argsort(ak.where(is_pi0, Tpi0, -np.inf), axis=1, ascending=False)
    for k in range(K_SLOTS):
        sel = order[:, k:k + 1]
        valid = ak.to_numpy(ak.fill_none(ak.firsts(is_pi0[sel]), False))
        g = lambda v: np.where(valid, ak.to_numpy(ak.fill_none(ak.firsts(v[sel]), np.nan)), np.nan)
        bx, by, bz = _from_beam(g(hx), g(hy), g(hz))
        d[f"tfs_pi0_T_{k}"] = g(Tpi0)
        d[f"tfs_pi0_px_{k}"], d[f"tfs_pi0_py_{k}"], d[f"tfs_pi0_pz_{k}"] = bx, by, bz
    # leading final-state photon by energy (makedf TruthSummary tfs_gamma_*_0; pi0s are undecayed)
    gi = ak.argmax(ak.where(is_g, hE, -np.inf), axis=1, keepdims=True)
    gvalid = ak.to_numpy(ak.any(is_g, axis=1))
    gg = lambda v: np.where(gvalid, ak.to_numpy(ak.fill_none(ak.firsts(v[gi]), np.nan)), np.nan)
    d["tfs_gamma_E_0"] = gg(hE)
    d["tfs_gamma_px_0"], d["tfs_gamma_py_0"], d["tfs_gamma_pz_0"] = _from_beam(gg(hx), gg(hy), gg(hz))
    # true muon summary in the beam frame
    p = np.sqrt(lep["px"] ** 2 + lep["py"] ** 2 + lep["pz"] ** 2)
    d["t_p_mu"] = p
    d["t_theta_mu"] = np.arccos(np.clip(np.divide(lep["pz"], p, out=np.ones_like(p), where=p > 0), -1, 1))
    d["t_pt_mu"] = np.hypot(lep["px"], lep["py"])
    d["t_pz_mu"] = lep["pz"]
    d["t_q0"] = d["mc_Enu"] - lep["E"]
    d["t_q3"] = np.sqrt(np.maximum(d["mc_Q2"] + d["t_q0"] ** 2, 0.0))
    d["t_N_pip"], d["t_N_pim"], d["t_N_pi0"] = d["tfs_npip"], d["tfs_npim"], d["tfs_npi0"]
    for w in ("w_fsi_elastic", "w_fsi_absorption", "geant_neutron_cv", "cv_wgt_tuple") + tuple(
            f"geant_{p}_{s}" for p in ("p", "n", "pi") for s in ("dn", "up")):
        d[w] = np.ones(n)
    d["weight"] = ak.to_numpy(a["weight"]).astype(float)
    return pd.DataFrame(d)


class GeneratorSample:
    """nuflow events (truth dataframe) plus run-level normalisation."""

    def __init__(self, path, max_events=None):
        self.path = path
        self.run = read_run(path)
        self.df = load_events(path, max_events)
        self.ev = Event(self.df, True, None, is_truth=True)
        self.weight = self.df["weight"].to_numpy()
        nucleons = float(self.run.get("nucleons_per_target", 13.0))
        # dsigma per event weight unit [cm^2 / nucleon]
        self.norm = float(self.run["flux_avg_xsec"]) * 1e-38 / nucleons / self.weight.sum()

    def predict(self, meas):
        """Differential cross section in meas.true_binning (inner bins, cm^2/nucleon/unit)."""
        TB = meas.true_binning
        sig = np.asarray(meas.signal(self.ev), dtype=bool)
        idx = TB.index(*meas.truth(self.ev))
        m = sig & (idx >= 0)
        h = np.bincount(idx[m], weights=self.weight[m], minlength=TB.size)
        return TB.inner(h) * self.norm / TB.widths()

    def signal_fraction(self, meas):
        return float(self.weight[np.asarray(meas.signal(self.ev), dtype=bool)].sum() / self.weight.sum())
