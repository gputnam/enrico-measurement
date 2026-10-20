"""N-dimensional binning with a linearized global bin index.

With flow=True every axis gets an underflow (index 0) and overflow (index n+1) bin, as in
ROOT; the global index is ROOT-like: ix + (nx+2)*(iy + (ny+2)*iz).  Results are usually
reported for in-range bins only (`inner()`).
"""
import numpy as np


class Binning:
    def __init__(self, edges, names=None, flow=True):
        self.edges = [np.asarray(e, dtype=float) for e in edges]
        self.names = list(names) if names else [f"x{i}" for i in range(len(self.edges))]
        self.flow = flow
        self.nb = [len(e) - 1 for e in self.edges]
        self.shape = tuple(n + 2 for n in self.nb) if flow else tuple(self.nb)
        self.size = int(np.prod(self.shape))

    @property
    def ndim(self):
        return len(self.edges)

    def axis_index(self, d, x):
        e = self.edges[d]
        i = np.searchsorted(e, x, side="right")  # 0 underflow, 1..n in range, n+1 overflow
        i = np.where(np.isnan(x), -1, i)
        if self.flow:
            return i
        return np.where((i >= 1) & (i <= self.nb[d]), i - 1, -1)

    def index(self, *xs):
        """Global bin index for each event; -1 for events outside (flow=False) or NaN."""
        assert len(xs) == self.ndim
        g = np.zeros(len(np.atleast_1d(xs[0])), dtype=np.int64)
        stride = 1
        bad = np.zeros_like(g, dtype=bool)
        for d, x in enumerate(xs):
            i = self.axis_index(d, np.asarray(x, dtype=float))
            bad |= i < 0
            g += i * stride
            stride *= self.shape[d]
        return np.where(bad, -1, g)

    def inner_mask(self):
        """Boolean mask over global bins that are in range on every axis."""
        if not self.flow:
            return np.ones(self.size, dtype=bool)
        m = np.ones(self.shape[::-1], dtype=bool)  # numpy order: last axis = x
        for d in range(self.ndim):
            sl = [slice(None)] * self.ndim
            ax = self.ndim - 1 - d
            sl[ax] = 0
            m[tuple(sl)] = False
            sl[ax] = self.shape[d] - 1
            m[tuple(sl)] = False
        return m.ravel()

    def inner(self, a, axis=0):
        """Select in-range global bins along `axis`."""
        return np.compress(self.inner_mask(), a, axis=axis)

    def widths(self):
        """Bin volume (product of widths) for in-range global bins, in global order."""
        grids = np.meshgrid(*[np.diff(e) for e in self.edges], indexing="ij")
        return np.prod(grids, axis=0).ravel(order="F")  # x fastest

    def centers(self):
        """List of per-axis center arrays for in-range global bins (global order)."""
        grids = np.meshgrid(*[0.5 * (e[1:] + e[:-1]) for e in self.edges], indexing="ij")
        return [g.ravel(order="F") for g in grids]
