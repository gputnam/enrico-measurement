"""Cross-section extraction from accumulated histograms (MINERvA procedure).

per universe u:
  N_u      = data - sum_c s_c,u * B_c,u           (MC already scaled to data POT)
  U_u      = Unfold(N_u; migration_u, n_iter)   (universe migration as in MnvUnfold; CV if absent)
  eps_u    = effnum_u / effden_u
  sigma_u  = U_u / (eps_u * Phi_u * T_u * POT_data) * 1e4 / binvolume      [cm^2 / nucleon / bin unit]
Statistical covariance: data Poisson errors propagated analytically through the unfolding.
Band covariance: weighted spread of the universes about their weighted mean.
"""
import numpy as np

from .flux import flux_model
from .targets import tracker_nucleons
from .unfold import unfold
from .universes import bands


def flux_integrals(tune, universes, pot_by_playlist, emin=0.0, emax=100.0, muon_correlation=True):
    """Data-POT-weighted flux integral (per m^2 per POT) for every universe."""
    tot = sum(pot_by_playlist.values())
    cv, univ = 0.0, None
    for pl, pot in pot_by_playlist.items():
        c, u = tune.flux_for(pl, 1).integral(emin, emax)
        cv += c * pot / tot
        univ = u * pot / tot if univ is None else univ + u * pot / tot
    out = np.full(len(universes), cv)
    for j, u in enumerate(universes):
        if u.band == "Flux":
            out[j] = univ[u.index]
        elif u.band == "Muon_Energy_MINOS" and muon_correlation:
            out[j] = cv * tune.muon_flux_corr.integral_ratio(u.index)
    return out


def background_total(acc, meas, sample, scales=None):
    """(R, U) sum of background categories, optionally scaled per category (dict cat -> (R,U) or scalar)."""
    B = 0.0
    for c in meas.bkg_categories:
        b = acc[f"bkg/{sample}/{c}"]
        if scales is not None and c in scales:
            b = b * scales[c]
        B = B + b
    return B


def covariance_from_universes(x, universes, cv=None):
    """{band: cov} from x (nbins, U)."""
    out = {}
    for band, idx in bands(universes).items():
        X = x[:, idx]
        w = np.array([universes[i].cov_weight for i in idx], dtype=float)
        if len(idx) == 1:
            ref = x[:, 0] if cv is None else cv
            d = X - ref[:, None]
            out[band] = d @ d.T
            continue
        mean = (X * w).sum(1) / w.sum()
        d = X - mean[:, None]
        out[band] = (d * w) @ d.T / w.sum()
    return out


def extract(acc, meas, info, n_iter=None, bkg_scales=None, zmin=5980.0, zmax=8422.0, apothem=850.0,
            emin=0.0, emax=100.0, closure=False):
    """closure=True: replace data by the CV MC prediction (signal + backgrounds) and also return
    the generator-level cross section from the truth tree ('truth_xsec') for comparison."""
    universes, tune = info["universes"], info["tune"]
    RB, TB = meas.reco_binning, meas.true_binning
    n_iter = meas.n_iter if n_iter is None else n_iter
    s = meas.signal_sample
    pot = sum(info["pot_data"].values())
    data = acc[f"data/{s}"]
    if closure:
        data = acc[f"sig/{s}"][:, 0] + background_total(acc, meas, s)[:, 0]
    B = background_total(acc, meas, s, bkg_scales)
    N = data[:, None] - B
    mig = acc["mig"]
    if "mig_u" in acc:
        Uall = np.column_stack([unfold(N[:, u], acc["mig_u"].dense(mig.shape, u), n_iter) for u in range(N.shape[1])])
    else:
        Uall = unfold(N, mig, n_iter)
    U0, cov_unf = unfold(N[:, 0], mig, n_iter, cov_measured=np.clip(data, 0, None))
    effnum, effden = acc["effnum"], acc["effden"]
    eff = np.divide(effnum, effden, out=np.zeros_like(effnum), where=effden > 0)
    phi = flux_integrals(tune, universes, info["pot_data"], emin, emax)
    T0 = tracker_nucleons(zmin, zmax, apothem)
    T = np.array([T0 * u.norm.get("targets", 1.0) for u in universes])
    inner = TB.inner_mask()
    width = TB.widths()
    norm = (phi * T * pot)[None, :] * 1e-4  # flux per m^2 -> per cm^2
    xs_all = np.divide(Uall, eff * norm, out=np.zeros_like(Uall), where=eff > 0)[inner] / width[:, None]
    scale0 = np.divide(1.0, eff[:, 0] * norm[0, 0], out=np.zeros(len(eff)), where=eff[:, 0] > 0)
    xs0 = (U0 * scale0)[inner] / width
    J = (scale0[:, None] * cov_unf * scale0[None, :])[np.ix_(inner, inner)]
    cov_stat = J / np.outer(width, width)
    covs = covariance_from_universes(xs_all, universes, cv=xs0)
    cov_sys = sum(covs.values()) if covs else np.zeros_like(cov_stat)
    truth_xs = (effden[:, 0] / norm[0, 0])[inner] / width
    return dict(xsec=xs0, truth_xsec=truth_xs, universes=xs_all, cov_stat=cov_stat, cov_bands=covs, cov_sys=cov_sys,
                cov_total=cov_stat + cov_sys, unfolded=U0[inner], eff=eff[inner, 0], flux_integral=phi[0],
                n_targets=T0, pot=pot, bkgsub=N[:, 0], data=data, bkg=B[:, 0])
