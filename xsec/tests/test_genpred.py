"""Generator-prediction bridge (xsec.genpred).  The normalisation test needs a nuflow file:
NUFLOW_TEST_FILE=<model>/runs/x.root venv/bin/python -m pytest xsec/tests/test_genpred.py"""
import os

import numpy as np
import pytest

from xsec.binning import Binning
from xsec.genpred import _from_beam, int_type
from xsec.kinematics import rotate_to_beam


def test_beam_rotation_inverse():
    rng = np.random.default_rng(1)
    p = rng.normal(size=(3, 100))
    q = rotate_to_beam(*_from_beam(*p))
    assert np.allclose(np.array(q), p)


def test_int_type():
    assert list(int_type(np.array([200, 300, 400, 500, 600, 100, 0]))) == [1, 8, 2, 3, 3, 4, 0]


class _AllEnu:
    tune_name = "none"
    true_binning = Binning([np.linspace(0, 100, 201)], ["enu"])

    def signal(self, ev):
        return np.ones(len(ev), dtype=bool)

    def truth(self, ev):
        return (ev["mc_Enu"] / 1e3,)


@pytest.mark.skipif(not os.environ.get("NUFLOW_TEST_FILE"), reason="needs NUFLOW_TEST_FILE")
def test_normalisation():
    from xsec.genpred import GeneratorSample
    g = GeneratorSample(os.environ["NUFLOW_TEST_FILE"], max_events=200000)
    m = _AllEnu()
    p = g.predict(m)
    total = (p * m.true_binning.widths()).sum()
    expect = g.run["flux_avg_xsec"] * 1e-38 / g.run["nucleons_per_target"]
    assert total == pytest.approx(expect, rel=1e-9)
    # the muon is the primary lepton of every CC event and the beam-frame angle is recovered
    t = g.ev.true_muon()
    assert np.all(g.df["mc_current"] == 1)
    assert np.allclose(np.cos(t["theta"]), g.df["t_pz_mu"] / g.df["t_p_mu"])


def _write_nuflow(path, events):
    """Minimal nuflow `events` tree from a list of (Enu [GeV], particles[(pdg, E, px, py, pz)] in GeV, beam frame)."""
    import awkward as ak
    import uproot
    n = len(events)
    cols = {k: [] for k in ("pdg", "status", "E", "px", "py", "pz")}
    for _, parts in events:
        cols["pdg"].append([p[0] for p in parts])
        cols["status"].append([1] * len(parts))
        for i, k in enumerate(("E", "px", "py", "pz")):
            cols[k].append([float(p[i + 1]) for p in parts])
    one = np.ones(n)
    import warnings
    with uproot.recreate(path) as f, warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)  # uproot >= 5.6: default RNTuple notice
        f["events"] = {
            "weight": one, "channel": np.full(n, 400, dtype=np.int32), "nu_pdg": np.full(n, 14, dtype=np.int32),
            "target_pdg": np.full(n, 1000060120, dtype=np.int64), "Enu": np.array([e[0] for e in events]),
            "Q2": one * 0.5, "q0": one, "q3": one * 1.2, "W": one * 1.3, "hit_pdg": np.full(n, 2212, dtype=np.int32),
            "pdg": ak.Array(cols["pdg"]), "status": ak.Array(cols["status"]), "E": ak.Array(cols["E"]),
            "px": ak.Array(cols["px"]), "py": ak.Array(cols["py"]), "pz": ak.Array(cols["pz"]),
        }


def _pi0(px, py, pz):
    m = 0.134977
    return (111, np.sqrt(m * m + px * px + py * py + pz * pz), px, py, pz)


def test_load_events_pi0_and_photon_slots(tmp_path):
    from xsec.genpred import M_PI0_DF, load_events
    mu = (13, 3.0, 0.1, 0.2, 2.99)
    lo, hi = _pi0(0.1, 0.0, 0.2), _pi0(0.0, 0.3, 0.5)
    gam_lo, gam_hi = (22, 0.05, 0.0, 0.0, 0.05), (22, 0.2, 0.2, 0.0, 0.0)
    events = [
        (5.0, [mu, (2212, 1.0, 0.0, 0.1, 0.3), lo, hi, gam_lo, gam_hi]),  # 2 pi0 (lower-T first) + 2 photons
        (5.0, [mu, (211, 0.5, 0.1, 0.1, 0.4), hi]),                        # 1 pi0, no photon
        (5.0, [mu, (2212, 1.0, 0.0, 0.1, 0.3)]),                           # no pi0
    ]
    path = str(tmp_path / "nuflow.root")
    _write_nuflow(path, events)
    df = load_events(path)
    assert list(df["tfs_npi0"]) == [2, 1, 0]
    # leading pi0 slot = higher kinetic energy, momenta in the detector frame (_from_beam), MeV
    for k, p in ((0, hi), (1, lo)):
        bx, by, bz = _from_beam(np.array([p[2]]), np.array([p[3]]), np.array([p[4]]))
        assert df[f"tfs_pi0_T_{k}"][0] == pytest.approx(p[1] * 1e3 - M_PI0_DF, abs=1e-6)
        assert df[f"tfs_pi0_px_{k}"][0] == pytest.approx(bx[0] * 1e3)
        assert df[f"tfs_pi0_py_{k}"][0] == pytest.approx(by[0] * 1e3)
        assert df[f"tfs_pi0_pz_{k}"][0] == pytest.approx(bz[0] * 1e3)
    # padding: NaN (as makedf); MINERvA-only Truth-tree pi0 columns are not provided
    assert np.isnan(df["tfs_pi0_T_1"][1]) and np.isnan(df["tfs_pi0_px_0"][2])
    assert not any(c.startswith(("t_pi0_", "t_gamma_")) for c in df.columns)
    # charged-pion slots are unaffected by pi0s
    assert df["tfs_pi_pdg_0"][1] == 211 and df["tfs_pi_pdg_0"][0] == 0
    # leading final-state photon by energy
    assert df["tfs_gamma_E_0"][0] == pytest.approx(200.0)
    bx, by, bz = _from_beam(np.array([0.2]), np.array([0.0]), np.array([0.0]))
    assert (df["tfs_gamma_px_0"][0], df["tfs_gamma_py_0"][0], df["tfs_gamma_pz_0"][0]) == pytest.approx(
        (bx[0] * 1e3, by[0] * 1e3, bz[0] * 1e3))
    assert np.isnan(df["tfs_gamma_E_0"][1]) and df["tfs_ngamma"][0] == 2
    # eavail counts the total pi0 energy (unchanged v1 definition)
    assert df["tfs_eavail"][1] == pytest.approx((0.5 - 0.13957) * 1e3 + hi[1] * 1e3, rel=1e-6)
