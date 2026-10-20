"""Reconstructed and true kinematics, including lateral systematic shifts.

`Shift` describes a lateral universe; `Event` wraps one dataframe chunk and returns
quantities in a given shift.  Units: MeV, mm, radians.
"""
from dataclasses import dataclass, field

import numpy as np

M_MU, M_P, M_N, M_PI = 105.6583745, 938.272013, 939.56536, 139.57018
BEAM_ANGLE = -0.05887  # rad, NuMI beam pitch (MinervaUnits::numi_beam_angle_rad)
NUCLEAR_TARGET_Z_END = 5835.0


@dataclass(frozen=True)
class Shift:
    """Lateral systematic state (all zero = central value)."""
    name: str = "CV"
    minerva_nsig: float = 0.0     # Muon_Energy_MINERvA
    minos_nsig: float = 0.0       # Muon_Energy_MINOS
    res_nsig: float = 0.0         # Muon_Energy_Resolution (0.4% of p_true - p_reco)
    beam_dx: float = 0.0          # BeamAngleX offset [rad]
    beam_dy: float = 0.0          # BeamAngleY offset [rad]
    angx_res_nsig: float = 0.0    # MuonAngleXResolution (2% of thetaX_true - thetaX_reco)
    angy_res_nsig: float = 0.0
    recoil: tuple = ()            # ((particle, nsig), ...) for recoil response universes
    proton: str = ""              # proton dE/dx variant ("bb_up", "meu_dn", ...)
    extra: tuple = field(default=())  # free-form (key, value) pairs for measurement-specific shifts

    def get(self, key, default=0.0):
        return dict(self.extra).get(key, default)


CV = Shift()


def theta3d(tx, ty):
    sx, sy = 1.0 / np.cos(tx), 1.0 / np.cos(ty)
    return np.arccos(1.0 / np.sqrt(sx * sx + sy * sy - 1.0))


def rotate_to_beam(px, py, pz):
    c, s = np.cos(BEAM_ANGLE), np.sin(BEAM_ANGLE)
    return px, c * py - s * pz, s * py + c * pz


class Event:
    """Accessor over one dataframe (reco or truth rows) with per-shift caching."""

    def __init__(self, df, is_mc, playlist=None, minos_p_offset=0.0, is_truth=False):
        self.df = df
        self.is_mc = is_mc
        self.is_truth = is_truth
        self.playlist = playlist
        self.minos_p_offset = minos_p_offset
        self._cache = {}

    def __len__(self):
        return len(self.df)

    def __getitem__(self, k):
        return self.df[k].to_numpy()

    def has(self, k):
        return k in self.df.columns

    # ------------------------------------------------------------ muon
    def muon(self, s=CV):
        key = ("muon", s.minerva_nsig, s.minos_nsig, s.res_nsig, s.beam_dx, s.beam_dy, s.angx_res_nsig, s.angy_res_nsig)
        if key in self._cache:
            return self._cache[key]
        d = self.df
        lx, ly, lz = d["mu_lepE_0"].to_numpy(np.float64), d["mu_lepE_1"].to_numpy(np.float64), d["mu_lepE_2"].to_numpy(np.float64)
        p_nom = np.sqrt(lx * lx + ly * ly + lz * lz)
        pminos_nom = d["minos_p"].to_numpy(np.float64)
        pminerva = p_nom - pminos_nom
        pminos = pminos_nom * (1.0 + self.minos_p_offset)
        if s.minerva_nsig:
            in_targets = d["vtx_2"].to_numpy() < NUCLEAR_TARGET_Z_END
            err = np.where(in_targets, np.hypot(40.0, 17.0), np.hypot(30.0, 11.0))
            pminerva = pminerva + s.minerva_nsig * err
        if s.minos_nsig:
            curv = d["minos_used_curv"].to_numpy() > 0
            loc = np.where(pminos > 1000.0, 0.006, 0.025)
            err = np.sqrt((pminos * 0.00984) ** 2 + np.where(curv, (pminos * loc) ** 2, 0.0))
            pminos_shifted = pminos + s.minos_nsig * err
        else:
            pminos_shifted = pminos
        p = pminerva + pminos_shifted
        tx = d["mu_thetaX"].to_numpy(np.float64) + s.beam_dx
        ty = d["mu_thetaY"].to_numpy(np.float64) + s.beam_dy
        if s.res_nsig or s.angx_res_nsig or s.angy_res_nsig:
            t = self.true_muon()
            if s.res_nsig:
                p = p + (t["p"] - p) * s.res_nsig * 0.004
            if s.angx_res_nsig:
                tx = tx + (t["thetaX"] - tx) * s.angx_res_nsig * 0.02
            if s.angy_res_nsig:
                ty = ty + (t["thetaY"] - ty) * s.angy_res_nsig * 0.02
        th = theta3d(tx, ty)
        out = dict(p=p, E=np.sqrt(p * p + M_MU * M_MU), theta=th, pt=p * np.sin(th), pz=p * np.cos(th),
                   thetaX=tx, thetaY=ty, p_minos=pminos_shifted, p_minerva=pminerva)
        self._cache[key] = out
        return out

    def true_muon(self):
        if "true_muon" in self._cache:
            return self._cache["true_muon"]
        d = self.df
        px, py, pz = rotate_to_beam(d["mc_lep_0"].to_numpy(np.float64), d["mc_lep_1"].to_numpy(np.float64),
                                    d["mc_lep_2"].to_numpy(np.float64))
        p = np.sqrt(px * px + py * py + pz * pz)
        th = np.arccos(np.clip(np.divide(pz, p, out=np.ones_like(p), where=p > 0), -1, 1))
        dx = np.sqrt(px * px + pz * pz)
        dy = np.sqrt(py * py + pz * pz)
        thx = np.sign(px) * np.arccos(np.clip(np.divide(pz, dx, out=np.ones_like(dx), where=dx > 0), -1, 1))
        thy = np.sign(py) * np.arccos(np.clip(np.divide(pz, dy, out=np.ones_like(dy), where=dy > 0), -1, 1))
        thx = np.where(px > 0, thx, -np.abs(thx))
        thy = np.where(py > 0, thy, -np.abs(thy))
        out = dict(p=p, theta=th, pt=p * np.sin(th), pz=p * np.cos(th), E=d["mc_lep_3"].to_numpy(np.float64),
                   thetaX=thx, thetaY=thy)
        self._cache["true_muon"] = out
        return out

    # ------------------------------------------------------------ truth helpers
    def q0q3_true(self):
        """MAT convention (TruthFunctions.h): q0 = Enu - Elep, q3 = sqrt(Q2 + q0^2).  MeV."""
        d = self.df
        q0 = d["mc_Enu"].to_numpy(np.float64) - d["mc_lep_3"].to_numpy(np.float64)
        q3 = np.sqrt(np.maximum(d["mc_Q2"].to_numpy(np.float64) + q0 * q0, 0.0))
        return q0, q3

    def true_vertex_in_fv(self, zmin=5980.0, zmax=8422.0, apothem=850.0):
        d = self.df
        return in_fv(d["mc_vtx_0"].to_numpy(), d["mc_vtx_1"].to_numpy(), d["mc_vtx_2"].to_numpy(), zmin, zmax, apothem)

    def reco_vertex_in_fv(self, zmin=5980.0, zmax=8422.0, apothem=850.0):
        d = self.df
        return in_fv(d["vtx_0"].to_numpy(), d["vtx_1"].to_numpy(), d["vtx_2"].to_numpy(), zmin, zmax, apothem)


def recoil_shift(ev, shift, prefix):
    """Lateral recoil-response shift [MeV]: sum_p nsig_p * max(err_p, 0), err from columns <prefix><p>
    (MAT ResponseUniverse, legacy error branches).  Data rows get no shift."""
    if not shift.recoil or not ev.is_mc:
        return 0.0
    out = 0.0
    for particle, nsig in shift.recoil:
        out = out + nsig * np.maximum(ev[prefix + particle], 0.0)
    return out


def in_hexagon(x, y, apothem=850.0):
    x, y = np.abs(x), np.abs(y)
    side = apothem * 2.0 / np.sqrt(3.0)
    return (x * x + y * y < apothem * apothem) | ((x <= apothem) & (y < side - x / np.sqrt(3.0)))


def in_fv(x, y, z, zmin, zmax, apothem):
    return (z >= zmin) & (z <= zmax) & in_hexagon(x, y, apothem)


def enu_qe(Emu, pmu, theta, Eb=34.0, neutrino=True):
    """CCQE neutrino energy (MeV), MAT PhysicsVariables convention."""
    mn, mp = (M_N, M_P) if neutrino else (M_P, M_N)
    num = mp * mp - (mn - Eb) ** 2 - M_MU * M_MU + 2 * (mn - Eb) * Emu
    den = 2 * (mn - Eb - Emu + pmu * np.cos(theta))
    return num / den


def q2_qe(Emu, pmu, theta, Eb=34.0, neutrino=True):
    e = enu_qe(Emu, pmu, theta, Eb, neutrino)
    return 2 * e * (Emu - pmu * np.cos(theta)) - M_MU * M_MU
