import numpy as np

from xsec.binning import Binning
from xsec.targets import n_planes, tracker_nucleons
from xsec.unfold import unfold


def test_binning_index_and_inner():
    b = Binning([[0, 1, 2], [0, 10, 20, 30]], flow=True)
    assert b.size == 4 * 5
    g = b.index(np.array([0.5, 1.5, -1, 0.5]), np.array([5, 25, 5, 35]))
    # ix + 4*iy
    assert list(g) == [1 + 4 * 1, 2 + 4 * 3, 0 + 4 * 1, 1 + 4 * 4]
    assert b.inner_mask().sum() == 6
    w = b.widths()
    assert w.shape == (6,) and np.allclose(w, 10)
    cx, cy = b.centers()
    assert np.allclose(cx[:2], [0.5, 1.5]) and np.allclose(cy[::2], [5, 15, 25])


def test_unfold_identity_and_jacobian():
    rng = np.random.default_rng(1)
    nt, nr = 6, 7
    mig = rng.uniform(0, 1, (nr, nt)) + 5 * np.eye(nr, nt)
    truth = mig.sum(0)
    meas = mig.sum(1) * 1.1
    u, cov, J = unfold(meas, mig, 4, cov_measured=meas, return_jacobian=True)
    # finite-difference Jacobian
    eps = 1e-6
    Jn = np.zeros_like(J)
    for j in range(nr):
        d = np.zeros(nr)
        d[j] = eps * meas[j]
        Jn[:, j] = (unfold(meas + d, mig, 4) - unfold(meas - d, mig, 4)) / (2 * d[j])
    assert np.allclose(J, Jn, rtol=1e-5, atol=1e-8)
    assert np.allclose(cov, J @ np.diag(meas) @ J.T)
    # folding back the truth reproduces the measurement for a perfect prior
    u0 = unfold(mig.sum(1), mig, 10)
    assert np.allclose(u0, truth)


def test_unfold_multi_matches_single():
    rng = np.random.default_rng(2)
    mig = rng.uniform(0, 1, (5, 5)) + 3 * np.eye(5)
    M = rng.uniform(5, 10, (5, 3))
    multi = unfold(M, mig, 3)
    for k in range(3):
        assert np.allclose(multi[:, k], unfold(M[:, k], mig, 3))


def test_targets():
    assert n_planes(5980, 8422) == 108
    assert abs(tracker_nucleons() / 3.235e30 - 1) < 2e-3
