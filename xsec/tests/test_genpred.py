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
