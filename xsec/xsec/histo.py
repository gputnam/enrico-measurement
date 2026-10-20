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


class SparseCols:
    """Sparse 2D histogram per universe: sorted cell keys (n,) and values (n, U).

    Used for the per-universe migration matrix, whose dense form (n_reco x n_true x U) would not
    fit in memory for fine 3D binnings while only a small fraction of cells is populated."""

    def __init__(self, keys, vals):
        self.keys = np.asarray(keys, dtype=np.int64)
        self.vals = np.asarray(vals, dtype=np.float64)

    @classmethod
    def fill(cls, key, w, cols, U):
        """key: (n,) cell index (<0 dropped); w: (n, len(cols)) weights for universes `cols`."""
        ok = key >= 0
        uk, inv = np.unique(key[ok], return_inverse=True)
        vals = np.zeros((len(uk), U))
        vals[:, cols] = fill(len(uk), inv, w[ok])
        return cls(uk, vals)

    def __add__(self, other):
        keys = np.union1d(self.keys, other.keys)
        vals = np.zeros((len(keys), self.vals.shape[1]))
        vals[np.searchsorted(keys, self.keys)] += self.vals
        vals[np.searchsorted(keys, other.keys)] += other.vals
        return SparseCols(keys, vals)

    def dense(self, shape, u):
        out = np.zeros(int(np.prod(shape)))
        out[self.keys] = self.vals[:, u]
        return out.reshape(shape)


class Accumulator(dict):
    """dict of numpy arrays (or SparseCols) that adds element-wise."""

    def add(self, key, arr):
        if isinstance(arr, SparseCols):
            self[key] = self[key] + arr if key in self else arr
        elif key in self:
            self[key] = self[key] + arr
        else:
            self[key] = np.array(arr, dtype=np.float64)

    def merge(self, other):
        for k, v in other.items():
            self.add(k, v)
        return self
