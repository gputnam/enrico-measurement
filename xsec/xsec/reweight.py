"""Model and detector weights (ports of the MINERvA MAT reweighters, numpy-vectorized).

All inputs are numpy arrays; energies in GeV unless noted.
"""
import functools

import numpy as np

from .paths import aux, npz

NONRES_PI_WEIGHT = 0.43
NONRES_PI_SHIFT = 0.04


# ---------------------------------------------------------------- GENIE non-resonant pion
def is_nonres_pi(rvn1pi_m1, rvp1pi_m1):
    return (rvn1pi_m1 < 1.0) | (rvp1pi_m1 < 1.0)


def nonres_pi_weight(rvn1pi_m1, rvp1pi_m1):
    return np.where(is_nonres_pi(rvn1pi_m1, rvp1pi_m1), NONRES_PI_WEIGHT, 1.0)


def is_cc_res(intType, current):
    return (intType == 2) & (current == 1)


# deuterium pion fit (MnvTune v4): CC RES normalization 1.15, MaRES 1.12 -> 0.94
D2_MARES, GENIE_MARES, GENIE_MARES_1SIG = 0.94, 1.12, 0.22
D2_RESNORM = 1.15


def genie_par_reweight(knob_m1, knob_p1, new_val, cv_val, one_sigma):
    """MAT GetGenieParReweight: linear interpolation/extrapolation of the knob weights."""
    shift = (new_val - cv_val) / one_sigma
    return np.where(shift >= 0, 1 + shift * (knob_p1 - 1), 1 - shift * (knob_m1 - 1))


def d2_pion_tune_weight(intType, current, mares_m1, mares_p1):
    """CC resonant pion tune to deuterium data (MAT GENIEReweighter with useDeuteriumPionTune)."""
    w = D2_RESNORM * genie_par_reweight(mares_m1, mares_p1, D2_MARES, GENIE_MARES, GENIE_MARES_1SIG)
    return np.where(is_cc_res(intType, current), w, 1.0)


# ---------------------------------------------------------------- 2p2h low-recoil fit
_2P2H_FILES = {
    0: "fit-mec-2d-noScaleDown-penalty00300-best-fit.txt",       # CV
    1: "fit-mec-2d-nn-only-noScaleDown-penalty00300-best-fit.txt",  # nn+pp
    2: "fit-mec-2d-np-only-noScaleDown-penalty02000-best-fit.txt",  # np
    3: "fit-qe-gaussian-noScaleDown-penalty02000-best-fit.txt",      # QE
}


@functools.lru_cache(maxsize=None)
def _2p2h_params(variation):
    vals = []
    for line in open(aux("MParamFiles", "data", "Reweight", _2P2H_FILES[variation])):
        for tok in line.split():
            try:
                vals.append(float(tok))
            except ValueError:
                pass
    return vals[:6]


def _gauss2d(q0, q3, p):
    norm, mq0, mq3, sq0, sq3, corr = p
    z = ((q0 - mq0) / sq0) ** 2 + ((q3 - mq3) / sq3) ** 2 - 2 * corr * (q0 - mq0) * (q3 - mq3) / (sq0 * sq3)
    return 1 + norm * np.exp(-0.5 * z / (1 - corr * corr))


def lowrecoil_2p2h_weight(q0, q3, intType, targetZ, targetNucleon, variation=0):
    """MnvTune 2p2h enhancement.  q0, q3 in GeV.

    variation 0: CV (all 2p2h);  1: nn/pp pairs only;  2: np pairs only;  3: apply to QE instead.
    """
    q0 = np.asarray(q0, dtype=float)
    q3 = np.asarray(q3, dtype=float)
    w = _gauss2d(q0, q3, _2p2h_params(variation))
    if variation == 3:
        apply = (intType == 1) & (targetZ >= 2)
    else:
        apply = (intType == 8) & (targetZ >= 2)
        if variation in (1, 2):
            pair = targetNucleon - 2000000200  # 0: nn, 1: np, 2: pp
            if variation == 1:
                apply = apply & ((pair == 0) | (pair == 2))
            else:
                apply = apply & (pair == 1)
    return np.where(apply, w, 1.0)


# ---------------------------------------------------------------- Valencia RPA
_RPA_OFFSET_NX = {(14, 82): 35, (14, 26): 39, (14, 6): 27, (14, 8): 25,
                  (-14, 82): 33, (-14, 26): 40, (-14, 6): 30, (-14, 8): 29}
_RPA_FILE = {6: "12C", 8: "16O", 26: "56Fe", 82: "208Pb"}


class RPA:
    """weightRPA (NX processing, RPAMaterials on).  Weights for one (nu mode, target Z)."""

    def __init__(self, nu_mode=14, targetZ=6, use_nx=True, rpa_materials=True):
        tz = targetZ if (rpa_materials and targetZ in _RPA_FILE) else 6
        nu = "anu" if (rpa_materials and nu_mode < 0) else "nu"
        d = npz("reweight", f"outNievesRPAratio-{nu}{_RPA_FILE[tz]}-20GeV-20170202.npz")
        self.rel = d["hrelratio__c"]
        self.nonrel = d["hnonrelratio__c"]
        self.q2rel_x, self.q2rel = d["hQ2relratio__x"], d["hQ2relratio__c"]
        self.q2nonrel = d["hQ2nonrelratio__c"]
        if rpa_materials:
            self.offset = _RPA_OFFSET_NX.get((nu_mode, targetZ), 27) if use_nx else 10
        else:
            self.offset = 27 if use_nx else 10

    def _q2hist(self, q2, rel):
        c = self.q2rel if rel else self.q2nonrel
        b = np.searchsorted(self.q2rel_x, q2, side="right")  # ROOT bin number (1-based; 0 underflow)
        b = np.clip(b, 0, len(c) - 1)
        out = c[b]
        return np.where((q2 < 0) | (q2 > 9.0), 1.0, out)

    def weights(self, q0, q3):
        """Returns (cv, lowQ2_up, lowQ2_dn, highQ2_up, highQ2_dn).  q0, q3 in GeV."""
        q0 = np.asarray(q0, dtype=float)
        q3 = np.asarray(q3, dtype=float)
        Q2 = q3 * q3 - q0 * q0
        lim = 3000
        q3bin = np.where(q3 >= 3.0, lim - 1, (q3 * 1000.0).astype(np.int64))
        q0bin = np.where(q0 >= 3.0, lim - 1, (q0 * 1000.0).astype(np.int64))
        q0bin = np.where(q0 < 0.018, 18 + self.offset, q0bin)
        yb = q0bin - self.offset
        nx, ny = self.rel.shape

        def get(h, xb, yb):
            # TH2::GetBinContent clamps out-of-range bin numbers to under/overflow
            return h[np.clip(xb, 0, nx - 1), np.clip(yb, 0, ny - 1)]

        def lookup(h, rel):
            w = get(h, q3bin, yb)
            w = np.where(w <= 0.001, 1.0, w)
            alt = get(h, q3bin + 150, yb)
            w = np.where((q0 < 0.15) & (w > 0.9), alt, w)
            w = np.where(Q2 >= 9.0, 1.0, np.where(Q2 > 3.0, self._q2hist(Q2, rel), w))
            return np.where((w >= 0.001) & (w <= 2.0), w, 1.0)

        cv = lookup(self.rel, True)
        ext = lookup(self.nonrel, False)
        low_up = np.where(cv < 1, cv + 0.25 * (1 - cv), cv)
        low_dn = np.where(cv < 1, cv - 0.25 * (1 - cv), cv)
        r = 0.6
        hp = cv + r * (ext - cv)
        hp = hp + np.where(Q2 < 0.9, 1.5 * (0.9 - Q2) * (ext - hp), 0.0)
        hp = np.minimum(hp, ext)
        hp = np.maximum(hp, cv + 0.03)
        hm = cv - r * (ext - cv)
        hm = np.minimum(hm, cv - 0.03)
        hm = np.where((Q2 > 1.0) & (hm < 1.0), 1.0, hm)
        return cv, low_up, low_dn, hp, hm


def rpa_weights(q0, q3, intType, targetZ, nu_mode=14, use_nx=True, rpa_materials=True):
    """(5, n) RPA weights; 1 for non-QE or targetZ < 6.  nu_mode is the analysis nu PDG."""
    q0 = np.asarray(q0, dtype=float)
    out = np.ones((5, len(q0)))
    apply = (intType == 1) & (targetZ >= 6)
    for z in np.unique(targetZ[apply]):
        m = apply & (targetZ == z)
        out[:, m] = np.array(_rpa(nu_mode, int(z), use_nx, rpa_materials).weights(q0[m], np.asarray(q3)[m]))
    return out


@functools.lru_cache(maxsize=None)
def _rpa(nu_mode, z, use_nx, rpa_materials):
    return RPA(nu_mode, z, use_nx, rpa_materials)


# ---------------------------------------------------------------- MINOS efficiency (ME)
_MINOS_NU = dict(lo=[0.705549, 0.499846, -0.364875, 0.131063, -0.0228281, 0.00154069],
                 hi=[0.959697, -0.139194, 0.175986, -0.0802815, 0.0162377, -0.00122259], pot=(3.9385, 8.0311))
_MINOS_NUBAR = dict(lo=[1.06348, 0.49754, -1.07352, 0.810622], hi=[1.04643, 0.37246, -0.75515, 0.571880],
                    pot=(6.74048, 8.26951))


def batch_pot(numi_pot, batch_structure, vtx_batch):
    """Per-batch POT (1e12) used for the ME MINOS efficiency correction."""
    bp = np.where(batch_structure == 0, numi_pot / 6.0, 0.0)
    bp = np.where(batch_structure == 1, np.where(vtx_batch < 3, numi_pot / 4.0, numi_pot / 8.0), bp)
    bp = np.where(batch_structure == 2, np.where(vtx_batch < 5, numi_pot / 5.0, numi_pot / 10.0), bp)
    return bp


def minos_efficiency_correction(p_minos_gev, bpot, neutrino=True):
    p = np.clip(p_minos_gev, 1.0, 4.0)
    if neutrino:
        c = _MINOS_NU
        lo = sum(a * p ** i for i, a in enumerate(c["lo"]))
        hi = sum(a * p ** i for i, a in enumerate(c["hi"]))
    else:
        c = _MINOS_NUBAR
        a = c["lo"]
        lo = a[0] - a[1] / p - a[2] / p ** 2 - a[3] / p ** 3
        a = c["hi"]
        hi = a[0] - a[1] / p - a[2] / p ** 2 - a[3] / p ** 3
    plo, phi = c["pot"]
    # solve 1 + A pot + B pot^2 through (plo, lo) and (phi, hi)
    a1, b1, c1 = phi, phi * phi, hi - 1.0
    a2, b2, c2 = plo, plo * plo, lo - 1.0
    det = a1 * b2 - b1 * a2
    A = (c1 * b2 - b1 * c2) / det
    B = (a1 * c2 - c1 * a2) / det
    return 1.0 + A * bpot + B * bpot * bpot


def minos_efficiency_error(theta_deg, neutrino=True):
    t = np.clip(theta_deg, 0.0, 40.0)
    err2 = 0.01 ** 2 + (neutrino * (1.0 - (0.99079 - (-4.24729e-4) * t - 7.20874e-5 * t * t))) ** 2
    return np.sqrt(err2)


# ---------------------------------------------------------------- low-Q2 pion suppression (MnvTune v2 / v4.3)
def lowq2pi_weight(Q2_gev2, intType, current, channel="JOINT", shift=0.0):
    d = npz("reweight", "lowQ2pi_weights.npz")
    gx, gy = d[f"{channel}_weights__gx"], d[f"{channel}_weights__gy"]
    sx, sy = d[f"{channel}_weight_shifts__gx"], d[f"{channel}_weight_shifts__gy"]
    q2 = np.asarray(Q2_gev2, dtype=float)
    w = np.interp(q2, gx, gy) + shift * np.interp(q2, sx, sy)
    apply = is_cc_res(intType, current) & (q2 >= 0) & (q2 < 0.7)
    return np.where(apply, w, 1.0)


# ---------------------------------------------------------------- z-expansion axial form factor (MnvTune vX.4)
def _zexp():
    d = npz("reweight", "Z_Expansion_Reweight_v2126.npz")
    return d["zexpvar__x"], d["zexpvar__c"][1:-1], d["zexpvar__band_zexpCov"]


def zexp_weight(Q2_gev2, intType, universe=None):
    """CV (universe=None) or universe z-expansion weight for CCQE; Q2 clamped at 8 GeV^2."""
    from .flux import th1_interpolate
    x, cv, band = _zexp()
    q2 = np.minimum(np.asarray(Q2_gev2, dtype=float), 8.0)
    h = cv if universe is None else band[universe]
    return np.where(intType == 1, th1_interpolate(x, h, q2), 1.0)


# ---------------------------------------------------------------- coherent pion reweight (MnvTune vX.1-3)
def coh_weight(theta_deg, tpi_gev, intType):
    """weightCoherentPi::get_combined_weight (theta in deg, pion KE in GeV), COH events only."""
    d = npz("reweight", "cohPi_weights_tracker_ME.npz")
    ex, ec = d["h1d_epi_xsec_ratio__x"], d["h1d_epi_xsec_ratio__c"]
    tx, tc = d["h1d_theta_xsec_ratio__x"], d["h1d_theta_xsec_ratio__c"]

    def content(edges, c, x):  # GetBinContent(FindBin(x)) incl. under/overflow
        return c[np.clip(np.searchsorted(edges, x, side="right"), 0, len(c) - 1)]
    th, e = np.asarray(theta_deg, float), np.asarray(tpi_gev, float)
    we, wt = content(ex, ec, e), content(tx, tc, th)
    hi_t, hi_e = th > tx[-1], e > ex[-1]
    w = np.where(hi_t & hi_e, 1.0, np.where(hi_t, we, np.where(hi_e, wt, we * wt)))
    return np.where(intType == 4, w, 1.0)
