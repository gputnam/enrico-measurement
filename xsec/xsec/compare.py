"""Comparisons of extracted cross sections with published results."""
import numpy as np


def chi2(x, ref, cov, mask=None, rcond=1e-12):
    """chi^2 of (x - ref) with covariance `cov` over bins in mask (pseudo-inverse)."""
    if mask is None:
        mask = np.ones(len(x), dtype=bool)
    d = (x - ref)[mask]
    c = cov[np.ix_(mask, mask)]
    return float(d @ np.linalg.pinv(c, rcond=rcond) @ d), int(mask.sum())


def summary(x, ref, cov_ref, cov_x=None, mask=None, label=""):
    if mask is None:
        mask = (ref != 0) & np.isfinite(x)
    r = x[mask] / ref[mask]
    c2, n = chi2(x, ref, cov_ref, mask)
    out = dict(label=label, nbins=n, chi2_refcov=c2, ratio_mean=float(np.mean(r)), ratio_median=float(np.median(r)),
               ratio_rms=float(np.std(r)), integral_ratio=float(x[mask].sum() / ref[mask].sum()))
    if cov_x is not None:
        out["chi2_sumcov"] = chi2(x, ref, cov_ref + cov_x, mask)[0]
    return out


def fmt(s):
    keys = [k for k in s if k != "label"]
    return f"{s['label']:<28}" + "  ".join(f"{k}={s[k]:.4g}" if isinstance(s[k], float) else f"{k}={s[k]}" for k in keys)
