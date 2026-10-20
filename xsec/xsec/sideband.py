"""Background-template fits to signal + sideband samples (per universe).

For each fit cell (a group of reco bins) and universe, solve
    D_s - O_s = sum_c f_c T_{s,c}         s in samples, c in fit categories
by weighted least squares (weights 1/max(D_s, 1)) using an SVD that drops singular values
below `rcond` x the largest (MINERvA: condition number cut 1e-3); negative scales are set to 0.
O_s is the sum of unfitted categories.  Returns {category: (R, U) scale} for the reco bins.
"""
import numpy as np


def fit_scales(acc, meas, cell_of_bin, fit_cats, samples=None, rcond=1e-3, keep=None):
    """cell_of_bin: (R,) int cell id per reco global bin (-1 = not fitted, scale 1).
    keep: optional (R,) bool, reco bins that enter the cell sums (default all with cell >= 0)."""
    samples = samples or meas.samples
    R, U = acc[f"sig/{samples[0]}"].shape
    ncell = cell_of_bin.max() + 1
    use = cell_of_bin >= 0 if keep is None else (cell_of_bin >= 0) & keep

    def cellsum(h):  # (R,) or (R,U) -> (ncell,) / (ncell,U)
        idx = cell_of_bin[use]
        if h.ndim == 1:
            return np.bincount(idx, weights=h[use], minlength=ncell)
        out = np.zeros((ncell, h.shape[1]))
        np.add.at(out, idx, h[use])
        return out

    S, C = len(samples), len(fit_cats)
    T = np.zeros((ncell, U, S, C))
    Y = np.zeros((ncell, U, S))
    Wt = np.zeros((ncell, S))
    for i, s in enumerate(samples):
        d = cellsum(acc[f"data/{s}"])
        Wt[:, i] = 1.0 / np.sqrt(np.maximum(d, 1.0))
        other = np.zeros((ncell, U))
        for c in meas.bkg_categories:
            if c not in fit_cats:
                other += cellsum(acc[f"bkg/{s}/{c}"])
        Y[:, :, i] = d[:, None] - other
        for j, c in enumerate(fit_cats):
            key = f"sig/{s}" if c == "signal" else f"bkg/{s}/{c}"
            T[:, :, i, j] = cellsum(acc[key])
    A = T * Wt[:, None, :, None]
    b = Y * Wt[:, None, :]
    u_, sv, vt = np.linalg.svd(A, full_matrices=False)  # (ncell,U,S,C), (ncell,U,C), (ncell,U,C,C)
    keep_sv = sv > rcond * sv[..., :1]
    inv = np.where(keep_sv, 1.0 / np.where(sv > 0, sv, 1.0), 0.0)
    coef = np.einsum("nuij,nui,nuki,nuk->nuj", vt, inv, u_, b)  # V diag(1/s) U^T b
    coef = np.where(np.isfinite(coef), np.maximum(coef, 0.0), 1.0)
    empty = (T.sum(2) <= 0)  # category absent in a cell -> scale 1
    coef = np.where(empty, 1.0, coef)
    out = {}
    for j, c in enumerate(fit_cats):
        sc = np.ones((R, U))
        m = cell_of_bin >= 0
        sc[m] = coef[cell_of_bin[m], :, j]
        out[c] = sc
    return out, coef
