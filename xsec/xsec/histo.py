"""Universe-aware histogram filling."""
import numpy as np


def fill(nbins, idx, w):
    """Weighted histogram over global bins; idx < 0 dropped.  w: (n,) or (n, K) -> (nbins,) or (nbins, K)."""
    idx = np.asarray(idx)
    ok = idx >= 0
    if w.ndim == 1:
        return np.bincount(idx[ok], weights=w[ok], minlength=nbins).astype(np.float64)
    i, W = idx[ok], w[ok]
    K = W.shape[1]
    flat = (i[:, None] * K + np.arange(K)[None, :]).ravel()
    return np.bincount(flat, weights=W.ravel(), minlength=nbins * K).reshape(nbins, K)


def fill2d(nx, ny, ix, iy, w):
    ok = (ix >= 0) & (iy >= 0)
    return np.bincount(ix[ok] * ny + iy[ok], weights=w[ok], minlength=nx * ny).reshape(nx, ny)


class Accumulator(dict):
    """dict of numpy arrays that adds element-wise."""

    def add(self, key, arr):
        if key in self:
            self[key] = self[key] + arr
        else:
            self[key] = np.array(arr, dtype=np.float64)

    def merge(self, other):
        for k, v in other.items():
            self.add(k, v)
        return self
