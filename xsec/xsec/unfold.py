"""Iterative Bayesian (D'Agostini) unfolding with analytic error propagation.

Response convention: migration[r, t] = weighted count of selected signal events with reco
bin r and true bin t.  Efficiency is handled outside the unfolding (MINERvA convention), so
the response is normalized per true bin over reco bins.  The first-iteration prior is the
truth projection of the migration matrix.  Error propagation follows Adye (RooUnfoldBayes):
the full Jacobian d(unfolded)/d(measured) is propagated through the iterations.
"""
import numpy as np


def unfold(measured, migration, n_iter, prior=None, cov_measured=None, return_jacobian=False):
    """Unfold `measured` (n_reco,) with `migration` (n_reco, n_true).

    Returns unfolded (n_true,) and, if cov_measured is given, its covariance (n_true, n_true).
    measured may also be (n_reco, n_univ) to unfold many universes with the same response
    (no error propagation in that case).
    """
    mig = np.asarray(migration, dtype=float)
    tsum = mig.sum(0)
    M = np.divide(mig, tsum, out=np.zeros_like(mig), where=tsum > 0)  # P(r | t)
    eff = M.sum(0)  # 1 where the true bin is populated, 0 otherwise
    u = (tsum if prior is None else np.asarray(prior, dtype=float)).copy()
    n = np.asarray(measured, dtype=float)
    multi = n.ndim == 2
    if multi:
        u = np.repeat(u[:, None], n.shape[1], axis=1)
    track = (cov_measured is not None or return_jacobian) and not multi
    J = np.zeros((mig.shape[1], mig.shape[0])) if track else None
    inv_eff = np.divide(1.0, eff, out=np.zeros_like(eff), where=eff > 0)
    for _ in range(n_iter):
        f = M @ u  # folded prediction per reco bin
        r = np.divide(n, f, out=np.zeros_like(n), where=f > 0)
        A = M.T @ r  # sum_j M_ji n_j / f_j
        if multi:
            u_new = u * A * inv_eff[:, None]
        else:
            u_new = u * A * inv_eff
        if track:
            invf = np.divide(1.0, f, out=np.zeros_like(f), where=f > 0)
            # d u'_i = (A_i/e_i) du_i + (u_i/e_i) sum_j M_ji (dn_j/f_j - n_j/f_j^2 sum_l M_jl du_l)
            term = M.T * invf[None, :] - (M.T * (n * invf * invf)[None, :]) @ M @ J
            J = (A * inv_eff)[:, None] * J + (u * inv_eff)[:, None] * term
        u = u_new
    if not track:
        return u
    out = [u]
    if cov_measured is not None:
        cov_m = np.diag(cov_measured) if np.ndim(cov_measured) == 1 else np.asarray(cov_measured)
        out.append(J @ cov_m @ J.T)
    if return_jacobian:
        out.append(J)
    return tuple(out)
