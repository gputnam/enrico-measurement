"""MINERvA ME FHC QE-like triple-differential cross section d3sigma/(dpT dp|| dSumTp).

Reference: D. Ruterbories et al., Phys. Rev. Lett. 129, 021803 (2022), arXiv:2203.08022.
MnvTune v4.4.1, 10 D'Agostini iterations, background fit to 4 samples per (pT, SumTp) bin.

Signal: numu CC with only muon + nucleons (+ photons < 10 MeV) in the final state, true vertex
in the tracker FV, theta_mu < 20 deg, p|| > 1.5 GeV/c.  SumTp = summed true proton KE.
Reco: MINOS-matched mu- in FV, theta < 20 deg, <= 1 dead discriminator; every non-muon
track proton-like by dE/dx (Q2QE-dependent proton-score cut, CCQENu default).
Samples (Michel x isolated blobs): A (0, <=1) signal, B (>=1, <=1), C (0, >=2), D (>=1, >=2).
Isolated blobs: start z > 4750 mm and < 10 MeV per hit (denser clusters are neutron-induced
protons; the PRL wording is inverted, but this choice reproduces the published sample sizes
A 1.3M / B 0.23M / C 0.22M / D 57k to ~5%).
Recoil (non-vertex: nonvtx_iso_blobs_energy + dispr_id_energy) < 500 MeV, as in PRD 99 012004.
Reco SumTp proxy: calorimetric available energy as in MAT GetEAvailable,
1.17 x (blob_recoil_E_tracker + blob_recoil_E_ecal), which includes tracked protons and vertex energy.
"""
import numpy as np

from xsec.binning import Binning
from xsec.kinematics import CV, q2_qe, recoil_shift
from xsec.measurement import Measurement
from xsec.universes import response_universes, standard_universes

SUMTP = [0, 0.02, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.40, 0.60, 0.799]
PT = [0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 2.5]
PZ = [1.5, 3.5, 4.5, 7.0, 8.0, 10.0, 20.0]
ZMIN, ZMAX, APOTHEM = 5980.0, 8422.0, 850.0
BLOB_ZMIN = 4750.0
BLOB_MAX_E_PER_HIT = 10.0  # MeV
RECOIL_MAX = 0.5  # GeV
EAVAIL_SCALE = 1.17  # MAT LowRecoilFunctions Eavailable_scale
NBLOB_SLOTS = 5
NSECP_SLOTS = 6
# CCQENu default proton score bands: Q2QE < 0.2 -> score > 0.2; [0.2, 0.6) -> 0.1; >= 0.6 -> 0.0
PSCORE_Q2 = [0.2, 0.6]
PSCORE_MIN = [0.2, 0.1, 0.0]


def proton_score_min(q2qe_gev2):
    idx = np.searchsorted(PSCORE_Q2, q2qe_gev2, side="right")
    return np.asarray(PSCORE_MIN)[idx]


class QELike3D(Measurement):
    name = "ccqelike_3d_me_fhc"
    tune_name = "v4.4.1"
    tune_options = dict(fsi_elastic=True)
    n_iter = 10
    samples = ("A", "B", "C", "D")
    signal_sample = "A"
    bkg_categories = ("1chpi", "1pi0", "multipi", "other")
    fit_categories = ("signal", "1chpi", "1pi0", "multipi")
    reco_binning = Binning([SUMTP, PT, PZ], ["sumtp", "pt", "pz"])
    true_binning = Binning([SUMTP, PT, PZ], ["sumtp", "pt", "pz"])
    extra_reco_columns = (["nu_helicity", "multiplicity", "n_improved_michel", "n_isoblob_prong",
                           "nonvtx_iso_blobs_E", "dispr_id_E", "p_score", "n_secp",
                           "blob_recoil_trk", "blob_recoil_ecal"]
                          + [f"isoblob_z_{i}" for i in range(NBLOB_SLOTS)]
                          + [f"isoblob_E_{i}" for i in range(NBLOB_SLOTS)]
                          + [f"isoblob_nhits_{i}" for i in range(NBLOB_SLOTS)]
                          + [f"secp_score_{i}" for i in range(NSECP_SLOTS)])
    extra_mc_reco_columns = ([f"rerr_dispr_{p}" for p in ("proton", "meson", "em", "other")]
                             + [f"rerr_all_{p}" for p in ("proton", "meson", "em", "other")])
    extra_truth_columns = ["tfs_np", "tfs_npip", "tfs_npim", "tfs_npi0", "tfs_nkpm", "tfs_nk0", "tfs_nomeson",
                           "tfs_nhbaryon", "tfs_maxgammaE", "tfs_nlep", "tfs_sumTp", "w_fsi_elastic"]

    def __init__(self, fit_backgrounds=True):
        self.fit_backgrounds = fit_backgrounds

    def universes(self, tune):
        return standard_universes(tune) + response_universes()

    # ------------------------------------------------------------ reco
    def preselect(self, ev):
        return ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1) & (ev["has_vtx"] == 1)

    def n_blobs(self, ev):
        n = ev["n_isoblob_prong"]
        cnt = np.zeros(len(ev), dtype=int)
        for i in range(NBLOB_SLOTS):
            z = ev[f"isoblob_z_{i}"]
            eph = ev[f"isoblob_E_{i}"] / np.maximum(ev[f"isoblob_nhits_{i}"], 1)
            cnt += ((i < n) & (z > BLOB_ZMIN) & (eph < BLOB_MAX_E_PER_HIT)).astype(int)
        return cnt + np.maximum(n - NBLOB_SLOTS, 0)  # beyond the stored slots: assume in tracker

    def protons_ok(self, ev, mu):
        q2 = q2_qe(mu["E"], mu["p"], mu["theta"]) / 1e6
        q2 = np.where(np.isfinite(q2), q2, 0.0)
        smin = proton_score_min(np.maximum(q2, 0.0))
        multi = ev["multiplicity"] >= 2
        prim = (ev["p_score"] >= 0) & (ev["p_score"] >= smin)
        nsec = ev["n_secp"]
        sec_ok = np.ones(len(ev), dtype=bool)
        for i in range(NSECP_SLOTS):
            sc = ev[f"secp_score_{i}"]
            sec_ok &= ~((i < nsec) & (sc < smin))
        return ~multi | (prim & sec_ok)

    def reco(self, ev, shift):
        mu = ev.muon(shift)
        base = (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["has_vtx"] == 1) & (ev["minos_match"] == 1)
                & (ev["nu_helicity"] == 1) & (mu["theta"] < np.radians(20.0)) & (ev["dead"] <= 1)
                & self.protons_ok(ev, mu))
        michel = ev["n_improved_michel"] > 0
        blobs = self.n_blobs(ev) >= 2
        recoil = (ev["nonvtx_iso_blobs_E"] + ev["dispr_id_E"] + recoil_shift(ev, shift, "rerr_dispr_")) / 1000.0
        base = base & (recoil < RECOIL_MAX)
        samples = {"A": base & ~michel & ~blobs, "B": base & michel & ~blobs,
                   "C": base & ~michel & blobs, "D": base & michel & blobs}
        return (self.eavail(ev, shift), mu["pt"] / 1000.0, mu["pz"] / 1000.0), samples

    def eavail(self, ev, shift):
        e = np.maximum(ev["blob_recoil_trk"], 0.0) + np.maximum(ev["blob_recoil_ecal"], 0.0)
        return EAVAIL_SCALE * (e + recoil_shift(ev, shift, "rerr_all_")) / 1000.0

    # ------------------------------------------------------------ truth
    def truth(self, ev):
        t = ev.true_muon()
        # clamp: v1 DFs carry -O(keV) SumTp for protons at rest (mass rounding), which would fall in
        # the underflow instead of the first SumTp bin
        return np.maximum(ev["tfs_sumTp"], 0.0) / 1000.0, t["pt"] / 1000.0, t["pz"] / 1000.0

    def qelike_fs(self, ev):
        return ((ev["tfs_npip"] + ev["tfs_npim"] + ev["tfs_npi0"] + ev["tfs_nkpm"] + ev["tfs_nk0"]
                 + ev["tfs_nomeson"] + ev["tfs_nhbaryon"] + ev["tfs_nlep"]) == 0) & (ev["tfs_maxgammaE"] < 10.0)

    def signal(self, ev):
        t = ev.true_muon()
        return ((ev["mc_incoming"] == 14) & (ev["mc_current"] == 1) & self.qelike_fs(ev)
                & ev.true_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (t["theta"] < np.radians(20.0)) & (t["pz"] > 1500.0))

    def category(self, ev):
        npic = ev["tfs_npip"] + ev["tfs_npim"]
        npi0 = ev["tfs_npi0"]
        cc = (ev["mc_incoming"] == 14) & (ev["mc_current"] == 1)
        cat = np.full(len(ev), 3)
        cat[cc & (npic == 1) & (npi0 == 0)] = 0
        cat[cc & (npic == 0) & (npi0 == 1)] = 1
        cat[cc & (npic + npi0 >= 2)] = 2
        return cat


def fit_cells(binning):
    """Cell id per reco global bin = (SumTp, pT) bin including flow; p|| summed over in-range bins."""
    sx, sy, sz = binning.shape
    g = np.arange(binning.size)
    ix, iy, iz = g % sx, (g // sx) % sy, g // (sx * sy)
    cell = ix + sx * iy
    keep = (iz >= 1) & (iz <= sz - 2)
    return cell, keep


def background_scales(acc, meas, pz_points=True):
    """pz_points: in-range p|| bins of a (SumTp, pT) cell are separate fit points (PRL); False
    integrates over p|| (one point per sample and cell, LE/ME paper)."""
    from xsec.sideband import fit_scales
    cell, keep = fit_cells(meas.reco_binning)
    sx, sy, _ = meas.reco_binning.shape
    point = np.arange(meas.reco_binning.size) // (sx * sy) if pz_points else None
    scales, coef = fit_scales(acc, meas, cell, meas.fit_categories, keep=keep, point_of_bin=point)
    return {c: v for c, v in scales.items() if c != "signal"}, coef
