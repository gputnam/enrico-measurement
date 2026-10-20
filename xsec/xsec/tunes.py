"""MINERvA tune compositions (central-value event weight).

  v1     : flux (nu-e constrained) x non-res pi (0.43) x 2p2h low-recoil fit x RPA x MINOS eff.
  v4.x   : v1 x deuterium pion tune (CC RES norm 1.15, MaRES 0.94)
  v4.3.x : v4 x low-Q2 pion suppression (JOINT)    [COH reweight not yet applied]
  v4.4.x : v4 with z-expansion axial FF (GENIE_MaZExpCCQE universes replace MaCCQE)
  vX.Y.1 : elastic FSI bug fix (w_fsi_elastic column from makedf)
MINOS efficiency is a reco-level weight (1 for truth rows).
"""
import numpy as np

from . import reweight as rw
from .flux import MuonEnergyFluxCorrection, flux_model

TUNES = {
    "v1": dict(d2pion=False, lowq2pi=False),
    "v4.3.1": dict(d2pion=True, lowq2pi=True, coh=True, fsi_elastic=True),
    "v4.4.1": dict(d2pion=True, lowq2pi=False, zexp=True, fsi_elastic=True),
    # v4.7.1 = v4.3.1 x (Tpi, Q2) 1pi weight, which is not publicly available -> same as v4.3.1 here
    "v4.7.1": dict(d2pion=True, lowq2pi=True, coh=True, fsi_elastic=True),
    "genie": dict(nonrespi=False, twop2h=False, rpa=False, d2pion=False, lowq2pi=False),
}


class Tune:
    def __init__(self, name="v1", nu_pdg=14, nue_constraint=True, n_flux_universes=100, neutrino=True, **over):
        cfg = dict(nonrespi=True, twop2h=True, rpa=True, minos_eff=True, d2pion=False, lowq2pi=False, zexp=False,
                   fsi_elastic=False, geant_neutron_cv=False, coh=False)
        cfg.update(TUNES[name])
        cfg.update(over)
        self.name, self.cfg = name, cfg
        self.nu_pdg = nu_pdg
        self.neutrino = neutrino
        self.nue_constraint = nue_constraint
        self.n_flux_universes = n_flux_universes
        self.muon_flux_corr = MuonEnergyFluxCorrection(antinu=not neutrino)

    # ----------------------------------------------------------- flux
    def flux_for(self, playlist, pdg_sign):
        return flux_model(playlist, 14 * pdg_sign if abs(self.nu_pdg) == 14 else self.nu_pdg * pdg_sign,
                          self.nue_constraint, self.n_flux_universes)

    def _enu_sign(self, ev):
        enu = ev["mc_Enu"] / 1000.0
        sign = np.where(ev["mc_incoming"] > 0, 1, -1)
        return enu, sign

    def flux_cv(self, ev):
        enu, sign = self._enu_sign(ev)
        w = np.ones(len(enu))
        for sg in (1, -1):
            m = sign == sg
            if m.any():
                w[m] = self.flux_for(ev.playlist, sg).cv_weight(enu[m])
        return w

    def flux_universe_ratios(self, ev):
        enu, sign = self._enu_sign(ev)
        out = np.ones((self.n_flux_universes, len(enu)))
        for sg in (1, -1):
            m = sign == sg
            if m.any():
                out[:, m] = self.flux_for(ev.playlist, sg).universe_ratios(enu[m])
        return out

    def flux_universe_weights(self, playlist):
        return self.flux_for(playlist, 1).univ_w

    # ----------------------------------------------------------- model pieces (cached on the event)
    def components(self, ev):
        if "tune" in ev._cache:
            return ev._cache["tune"]
        c = {}
        q0, q3 = ev.q0q3_true()
        q0, q3 = q0 / 1000.0, q3 / 1000.0
        it, tz = ev["mc_intType"], ev["mc_targetZ"]
        c["flux"] = self.flux_cv(ev)
        c["nonrespi"] = rw.nonres_pi_weight(ev["g_Rvn1pi_m1"], ev["g_Rvp1pi_m1"]) if self.cfg["nonrespi"] else np.ones(len(ev))
        c["2p2h"] = rw.lowrecoil_2p2h_weight(q0, q3, it, tz, ev["mc_targetNucleon"], 0) if self.cfg["twop2h"] else np.ones(len(ev))
        if self.cfg["rpa"]:
            r = rw.rpa_weights(q0, q3, it, tz, nu_mode=self.nu_pdg)
            c["rpa"], c["rpa_vars"] = r[0], r
        else:
            c["rpa"] = np.ones(len(ev))
        c["d2pion"] = (rw.d2_pion_tune_weight(it, ev["mc_current"], ev["g_MaRES_m1"], ev["g_MaRES_p1"])
                       if self.cfg["d2pion"] else np.ones(len(ev)))
        c["lowq2pi"] = (rw.lowq2pi_weight(ev["mc_Q2"] / 1e6, it, ev["mc_current"])
                        if self.cfg["lowq2pi"] else np.ones(len(ev)))
        c["zexp"] = rw.zexp_weight(ev["mc_Q2"] / 1e6, it) if self.cfg["zexp"] else np.ones(len(ev))
        c["fsi"] = ev["w_fsi_elastic"] if self.cfg["fsi_elastic"] else np.ones(len(ev))
        if self.cfg["coh"]:
            px, py, pz = ev["tfs_pi_px_0"], ev["tfs_pi_py_0"], ev["tfs_pi_pz_0"]
            from .kinematics import rotate_to_beam
            px, py, pz = rotate_to_beam(px, py, pz)
            th = np.degrees(np.arccos(np.clip(pz / np.sqrt(px * px + py * py + pz * pz), -1, 1)))
            c["coh"] = rw.coh_weight(np.nan_to_num(th, nan=0.0), np.nan_to_num(ev["tfs_pi_T_0"] / 1000.0, nan=0.0), it)
        else:
            c["coh"] = np.ones(len(ev))
        c["geant_ncv"] = ev["geant_neutron_cv"] if self.cfg["geant_neutron_cv"] else np.ones(len(ev))
        ev._cache["tune"] = c
        return c

    def minos_eff(self, ev, shift=None):
        from .kinematics import CV
        if ev.is_truth or not self.cfg["minos_eff"]:
            return np.ones(len(ev))
        mu = ev.muon(shift or CV)
        bpot = rw.batch_pot(ev["numi_pot"], ev["batch_structure"], ev["vtx_batch"])
        return rw.minos_efficiency_correction(mu["p_minos"] / 1000.0, bpot, self.neutrino)

    def cv_weight(self, ev, shift=None):
        """CV weight; for reco rows the MINOS efficiency is evaluated in `shift`."""
        c = self.components(ev)
        w = c["flux"] * c["nonrespi"] * c["2p2h"] * c["rpa"] * c["d2pion"] * c["lowq2pi"] * c["zexp"] * c["fsi"] * c["geant_ncv"] * c["coh"]
        return w * self.minos_eff(ev, shift)
