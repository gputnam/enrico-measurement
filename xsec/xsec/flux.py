"""NuMI flux model following the MINERvA FluxReweighter conventions.

* CV event weight     = constrained_flux(Enu) / generated_flux(Enu)    (TH1::Interpolate)
* Flux universe ratio = universe_flux(Enu) / constrained_flux(Enu)
* nu-e constraint     : CV scaled bin-by-bin by (weighted mean / unweighted mean) of all
                        flux universes; universes keep their shape and carry the constraint
                        weights into the covariance.
Energies in GeV.  Fluxes are per m^2 per POT (per bin), as in the MINERvA files.
"""
import functools

import numpy as np

from .paths import aux, npz

# playlist -> flux file playlist string (MAT FluxReweighter::playlistString)
FLUX_PLAYLIST = {
    **{p: "minervame1D" for p in ["1A", "1B", "1C", "1D", "1E", "1F"]},
    **{p: "minervame1M" for p in ["1G", "1L", "1M"]},
    **{p: "minervame1N" for p in ["1N", "1O", "1P"]},
    "5A": "minervame5A", "6A": "minervame6A", "6B": "minervame6B", "6C": "minervame6C", "6D": "minervame6D",
    "6E": "minervame6E", "6F": "minervame6F", "6G": "minervame6G", "6H": "minervame6H", "6I": "minervame6I",
    "6J": "minervame6J",
    "ME_FHC_AVG": "minervame1D1M1NWeightedAve", "ME_RHC_AVG": "minervame5A6AWeightedAve",
    # LE (Eroica era); files exist for these
    "minerva1": "minerva1", "minerva5": "minerva5", "minerva7": "minerva1", "minerva9": "minerva1",
    "minerva13": "minerva13", "LE_FHC": "minervaLE-FHC",
}
ME_PLAYLISTS = set("1A 1B 1C 1D 1E 1F 1G 1L 1M 1N 1O 1P 5A 6A 6B 6C 6D 6E 6F 6G 6H 6I 6J".split())


def th1_interpolate(edges, vals, x):
    """ROOT TH1::Interpolate for one (nb,) or many (U, nb) histograms at points x.

    Linear between bin centers, constant beyond the first/last bin centers.
    """
    c = 0.5 * (edges[1:] + edges[:-1])
    x = np.asarray(x, dtype=float)
    i1 = np.clip(np.searchsorted(c, x, side="right"), 1, len(c) - 1)
    i0 = i1 - 1
    t = np.clip((x - c[i0]) / (c[i1] - c[i0]), 0.0, 1.0)
    vals = np.asarray(vals)
    return vals[..., i0] * (1 - t) + vals[..., i1] * t


def th1_content(edges, vals, x):
    """GetBinContent(FindBin(x)); 0 outside the axis range."""
    b = np.searchsorted(edges, x, side="right") - 1
    ok = (b >= 0) & (b < len(edges) - 1)
    out = np.asarray(vals)[..., np.clip(b, 0, len(edges) - 2)]
    return np.where(ok, out, 0.0)


def read_constraint(path):
    w = {}
    for line in open(path):
        p = line.split()
        if len(p) == 3 and p[0] == "Flux":
            w[int(p[1])] = float(p[2])
    out = np.zeros(max(w) + 1)
    for k, v in w.items():
        out[k] = v
    return out


class FluxModel:
    """Flux for one neutrino species (sign of `pdg`) and flux playlist."""

    def __init__(self, playlist, pdg=14, nue_constraint=True, n_universes=100):
        self.playlist = playlist
        self.pdg = pdg
        pl = FLUX_PLAYLIST[playlist]
        self.is_me = playlist in ME_PLAYLISTS or playlist.startswith("ME_")
        suffix = "_rearrangedUniverses" if (nue_constraint and self.is_me) else ""
        gen_version = "g4numiv6" if self.is_me else "g4numiv5"
        rw = npz("flux", f"flux-gen2thin-pdg{pdg}-{pl}{suffix}__flux_E_cvweighted.npz")
        gen = npz("flux", f"flux-{gen_version}-pdg{pdg}-{pl}{suffix}__flux_E_unweighted.npz")
        self.edges = rw["edges"]
        self.gen = gen["cv"]
        cv = rw["cv"]
        univ = rw["band_Flux"]
        w = np.ones(len(univ))
        if nue_constraint:
            cfile = "sorted_NuEConstraint_FHC_RHC_IMD.txt" if self.is_me else "nu+e_eroica.txt"
            cw = read_constraint(aux("MParamFiles", "data", "FluxConstraints", cfile))
            n = min(len(cw), len(univ))
            w = np.zeros(len(univ))
            w[:n] = cw[:n]
            used = w != 0  # universes without a constraint weight are not used in the means
            weighted = (w[:, None] * univ).sum(0) / w.sum()
            unweighted = univ[used].mean(0)
            corr = np.divide(weighted, unweighted, out=np.zeros_like(weighted), where=unweighted != 0)
            cv = cv * corr
        if n_universes and n_universes < len(univ):
            univ, w = univ[:n_universes], w[:n_universes]
        self.cv = cv
        self.univ = univ
        self.univ_w = w

    # ---- event weights
    def cv_weight(self, enu):
        enu = np.asarray(enu, dtype=float)
        num = th1_interpolate(self.edges, self.cv, enu)
        den = th1_interpolate(self.edges, self.gen, enu)
        hi = enu > 75.0
        if hi.any():
            num = np.where(hi, th1_content(self.edges, self.cv, enu), num)
            den = np.where(hi, th1_content(self.edges, self.gen, enu), den)
        ok = (num != 0) & (den != 0)
        return np.where(ok, num / np.where(den == 0, 1, den), 1.0)

    def universe_ratios(self, enu):
        """(n_universes, n) flux universe weight relative to the (constrained) CV."""
        enu = np.asarray(enu, dtype=float)
        num = th1_interpolate(self.edges, self.univ, enu)
        den = th1_interpolate(self.edges, self.cv, enu)
        hi = enu > 75.0
        if hi.any():
            num = np.where(hi, th1_content(self.edges, self.univ, enu), num)
            den = np.where(hi, th1_content(self.edges, self.cv, enu), den)
        ok = (num != 0) & (den != 0)
        return np.where(ok, num / np.where(den == 0, 1, den), 1.0)

    # ---- integrals
    def _integral(self, vals, emin, emax):
        # TH1::Integral(FindBin(emin), FindBin(emax), "width"): includes the full end bins
        e = self.edges
        b0 = max(np.searchsorted(e, emin, side="right") - 1, 0)
        b1 = min(np.searchsorted(e, emax, side="right") - 1, len(e) - 2)
        wid = np.diff(e)[b0:b1 + 1]
        return (np.asarray(vals)[..., b0:b1 + 1] * wid).sum(-1)

    def integral(self, emin=0.0, emax=100.0):
        """(cv, universes[n_univ]) integrated flux per m^2 per POT."""
        return self._integral(self.cv, emin, emax), self._integral(self.univ, emin, emax)


@functools.lru_cache(maxsize=None)
def flux_model(playlist, pdg=14, nue_constraint=True, n_universes=100):
    return FluxModel(playlist, pdg, nue_constraint, n_universes)


class MuonEnergyFluxCorrection:
    """Flux correlation of the MINOS muon energy scale (FluxReweighter::GetSysUniFluxWeightCorrection)."""

    def __init__(self, antinu=False):
        name = "MuonBar_Energy_MnvH1D" if antinu else "Muon_Energy_MnvH1D"
        d = npz("flux", f"{name}__{name}.npz")
        self.edges, self.cv, self.univ = d["edges"], d["cv"], d["band_Muon_Energy_MINOS"]

    def weight(self, enu, universe):
        """universe 0 = -1 sigma, 1 = +1 sigma.  1 above 15 GeV."""
        num = th1_interpolate(self.edges, self.univ[universe], enu)
        den = th1_interpolate(self.edges, self.cv, enu)
        ok = (num != 0) & (den != 0) & (np.asarray(enu) <= 15.0)
        return np.where(ok, num / np.where(den == 0, 1, den), 1.0)

    def integral_ratio(self, universe, emin=0.0, emax=100.0):
        e = self.edges
        wid = np.diff(e)
        return (self.univ[universe] * wid).sum() / (self.cv * wid).sum()
