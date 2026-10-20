"""Systematic universes (MINERvA conventions).

Each Universe has a band name and is either
  vertical : CV selection, event weight multiplied by ratio(ev)
  lateral  : reconstructed quantities recomputed in `shift` (optionally times a ratio)
  norm     : only the normalization changes (e.g. target mass)
Universe 0 is always the CV.  `cov_weight` is the universe weight in the band covariance
(nu-e constraint weights for the flux band).
"""
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from . import reweight as rw
from .kinematics import CV, Shift
from .targets import TARGET_MASS_ERR_CH

GENIE_STANDARD = ["AGKYxF1pi", "AhtBY", "BhtBY", "CCQEPauliSupViaKF", "CV1uBY", "CV2uBY", "EtaNCEL", "FrAbs_N",
                  "FrAbs_pi", "FrCEx_N", "FrCEx_pi", "FrElas_N", "FrElas_pi", "FrInel_N", "FrPiProd_N", "FrPiProd_pi",
                  "MFP_N", "MFP_pi", "MaNCEL", "NormDISCC", "NormNCRES", "RDecBR1gamma", "Rvn2pi", "Rvp2pi",
                  "Theta_Delta2Npi", "VecFFCCQEshape"]
GENIE_RES = ["MaRES", "MvRES"]


@dataclass
class Universe:
    band: str
    index: int
    kind: str  # "cv", "vert", "lat", "norm"
    shift: Shift = CV
    ratio: Optional[Callable] = None  # (ev, tune) -> array; applied to reco and truth rows
    reco_only: bool = False           # ratio is 1 for truth rows
    cov_weight: float = 1.0
    norm: dict = field(default_factory=dict)  # {"targets": factor}

    @property
    def name(self):
        return f"{self.band}_{self.index}"


def _genie(knob, sign):
    idx = "p1" if sign > 0 else "m1"

    def f(ev, tune):
        w = ev[f"g_{knob}_{idx}"].astype(np.float64)
        if knob in ("MFP_N", "MFP_pi"):
            w = np.where(w > 10.0, 5.0, np.where(w < 0.0, 0.0, w))
        return w
    return f


def _rvx1pi(knob, sign):
    idx = "p1" if sign > 0 else "m1"

    def f(ev, tune):
        w = ev[f"g_{knob}_{idx}"].astype(np.float64)
        if tune.cfg["nonrespi"]:
            nonres = rw.is_nonres_pi(ev["g_Rvn1pi_m1"], ev["g_Rvp1pi_m1"])
            w = np.where(nonres, 1 + sign * rw.NONRES_PI_SHIFT / rw.NONRES_PI_WEIGHT, w)
        return w
    return f


def _d2_mares(sign):
    """GENIE_D2_MaRES: MaRES 0.94 +- 0.05 relative to the tuned CV (CC RES only).

    Note: MAT's GenieMaResUniverse::GetWeightRatioToCV returns the CV GENIE weight (not 1) for
    non-CC-RES events, which looks like a bug; we use ratio 1 for non-RES events.
    """
    def f(ev, tune):
        m1, p1 = ev["g_MaRES_m1"], ev["g_MaRES_p1"]
        num = rw.genie_par_reweight(m1, p1, rw.D2_MARES + sign * 0.05, rw.GENIE_MARES, rw.GENIE_MARES_1SIG)
        den = rw.genie_par_reweight(m1, p1, rw.D2_MARES, rw.GENIE_MARES, rw.GENIE_MARES_1SIG)
        r = np.divide(num, den, out=np.ones(len(ev)), where=den != 0)
        return np.where(rw.is_cc_res(ev["mc_intType"], ev["mc_current"]), r, 1.0)
    return f


def _d2_norm(sign):
    def f(ev, tune):
        return np.where(rw.is_cc_res(ev["mc_intType"], ev["mc_current"]), (rw.D2_RESNORM + sign * 0.07) / rw.D2_RESNORM, 1.0)
    return f


def _ep_mvres(sign):
    """GENIE_EP_MvRES: MvRES 0.84 +- 3% (electroproduction fit), CC RES only."""
    def f(ev, tune):
        r = rw.genie_par_reweight(ev["g_MvRES_m1"], ev["g_MvRES_p1"], 0.84 * (1 + sign * 0.03), 0.84, 0.084)
        return np.where(rw.is_cc_res(ev["mc_intType"], ev["mc_current"]), r, 1.0)
    return f


def _flux(i):
    def f(ev, tune):
        if "flux_univ" not in ev._cache:
            ev._cache["flux_univ"] = tune.flux_universe_ratios(ev)
        return ev._cache["flux_univ"][i]
    return f


def _zexp(i):
    def f(ev, tune):
        q2, it = ev["mc_Q2"] / 1e6, ev["mc_intType"]
        cv = rw.zexp_weight(q2, it)
        return np.divide(rw.zexp_weight(q2, it, i), cv, out=np.ones(len(ev)), where=cv != 0)
    return f


def _rpa(var):
    def f(ev, tune):
        c = tune.components(ev)
        if "rpa_vars" not in c:
            return np.ones(len(ev))
        r = c["rpa_vars"]
        return np.divide(r[var], r[0], out=np.ones(len(ev)), where=r[0] != 0)
    return f


def _2p2h(var):
    def f(ev, tune):
        c = tune.components(ev)
        q0, q3 = ev.q0q3_true()
        w = rw.lowrecoil_2p2h_weight(q0 / 1000.0, q3 / 1000.0, ev["mc_intType"], ev["mc_targetZ"],
                                     ev["mc_targetNucleon"], var)
        return np.divide(w, c["2p2h"], out=np.ones(len(ev)), where=c["2p2h"] != 0)
    return f


def _minos_eff(sign):
    def f(ev, tune):
        mu = ev.muon(CV)
        th = np.degrees(mu["theta"])
        err = rw.minos_efficiency_error(th, tune.neutrino)
        cv = tune.minos_eff(ev, CV)
        return 1 + sign * err / cv
    return f


def _muon_minos_flux(sign):
    def f(ev, tune):
        return tune.muon_flux_corr.weight(ev["mc_Enu"] / 1000.0, 0 if sign < 0 else 1)
    return f


def standard_universes(tune, genie=True, flux=True, rpa=True, twop2h=True, muon=True, minos_eff=True,
                       beam_angle=True, muon_angle_res=False, target_mass=True, n_flux=None):
    """MINERvA-101-like standard set (no response / Geant bands)."""
    U = [Universe("CV", 0, "cv")]
    if flux:
        w = tune.flux_universe_weights("1D")
        for i in range(n_flux or tune.n_flux_universes):
            U.append(Universe("Flux", i, "vert", ratio=_flux(i), cov_weight=float(w[i])))
    if genie:
        knobs = GENIE_STANDARD + (GENIE_RES if not tune.cfg["d2pion"] else []) + ([] if tune.cfg["zexp"] else ["MaCCQE"])
        for k in knobs:
            for j, s in enumerate((-1, 1)):
                U.append(Universe(f"GENIE_{k}", j, "vert", ratio=_genie(k, s)))
        if tune.cfg["d2pion"]:
            for band, fn in (("GENIE_D2_MaRES", _d2_mares), ("GENIE_D2_NormCCRES", _d2_norm), ("GENIE_EP_MvRES", _ep_mvres)):
                for j, s in enumerate((-1, 1)):
                    U.append(Universe(band, j, "vert", ratio=fn(s)))
        if tune.cfg["zexp"]:
            for i in range(100):
                U.append(Universe("GENIE_MaZExpCCQE", i, "vert", ratio=_zexp(i)))
        for k in ("Rvn1pi", "Rvp1pi"):
            for j, s in enumerate((-1, 1)):
                U.append(Universe(f"GENIE_{k}", j, "vert", ratio=_rvx1pi(k, s)))
    if rpa and tune.cfg["rpa"]:
        for band, vars_ in (("RPA_HighQ2", (3, 4)), ("RPA_LowQ2", (1, 2))):
            for j, v in enumerate(vars_):
                U.append(Universe(band, j, "vert", ratio=_rpa(v)))
    if twop2h and tune.cfg["twop2h"]:
        for j, v in enumerate((1, 2, 3)):
            U.append(Universe("Low_Recoil_2p2h_Tune", j, "vert", ratio=_2p2h(v)))
    if muon:
        for j, s in enumerate((-1, 1)):
            U.append(Universe("Muon_Energy_MINERvA", j, "lat", shift=Shift(f"MINERvA{s:+d}", minerva_nsig=s)))
        for j, s in enumerate((-1, 1)):
            U.append(Universe("Muon_Energy_MINOS", j, "lat", shift=Shift(f"MINOS{s:+d}", minos_nsig=s),
                              ratio=_muon_minos_flux(s)))
        for j, s in enumerate((-1, 1)):
            U.append(Universe("Muon_Energy_Resolution", j, "lat", shift=Shift(f"Res{s:+d}", res_nsig=s)))
    if minos_eff and tune.cfg["minos_eff"]:
        for j, s in enumerate((-1, 1)):
            U.append(Universe("MINOS_Reconstruction_Efficiency", j, "vert", ratio=_minos_eff(s), reco_only=True))
    if beam_angle:
        for j, s in enumerate((-1, 1)):
            U.append(Universe("BeamAngleX", j, "lat", shift=Shift(f"BAX{s:+d}", beam_dx=s * 0.001)))
        for j, s in enumerate((-1, 1)):
            U.append(Universe("BeamAngleY", j, "lat", shift=Shift(f"BAY{s:+d}", beam_dy=s * 0.0009)))
    if muon_angle_res:
        for j, s in enumerate((-1, 1)):
            U.append(Universe("MuonAngleXResolution", j, "lat", shift=Shift(f"AXR{s:+d}", angx_res_nsig=s)))
        for j, s in enumerate((-1, 1)):
            U.append(Universe("MuonAngleYResolution", j, "lat", shift=Shift(f"AYR{s:+d}", angy_res_nsig=s)))
    if target_mass:
        for j, s in enumerate((-1, 1)):
            U.append(Universe("Target_Mass_CH", j, "norm", norm={"targets": 1 + s * TARGET_MASS_ERR_CH}))
    return U


def response_universes(particles=("proton", "meson", "em", "other")):
    """Recoil response lateral universes (band names follow MAT: response_<particle>)."""
    U = []
    for p in particles:
        for j, s in enumerate((-1, 1)):
            U.append(Universe(f"response_{p}", j, "lat", shift=Shift(f"resp_{p}{s:+d}", recoil=((p, s),))))
    return U


def bands(universes):
    """Ordered {band: [universe indices]} excluding the CV."""
    out = {}
    for i, u in enumerate(universes):
        if u.kind != "cv":
            out.setdefault(u.band, []).append(i)
    return out
