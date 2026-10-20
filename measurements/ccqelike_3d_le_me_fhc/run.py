#!/usr/bin/env python3
"""QE-like 3D cross section in the LE/ME-paper binning (arXiv:2606.00745), ME beam.

No numerical release exists for this binning.  As a quantitative check the result is projected
onto bins common with PRL 129 021803 (pT 0-0.15 merged, p|| 1.5-3.5 and 3.5-4.5 GeV/c) and
compared with the PRL tables projected the same way (integrated cross sections per coarse bin,
with the PRL covariance propagated).  Qualitative comparison to the paper figures is left to plots.

  run.py --mc <DF dir> --data <DF dir> [-j 16] [-o results.pkl]  |  run.py --load results.pkl
"""
import argparse
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analysis import PT, PZ, QELike3DLEBinning, background_scales  # noqa: E402

from xsec.util import load_measurement_module  # noqa: E402

prl_run = load_measurement_module("ccqelike_3d_me_fhc", "run")
_prl = load_measurement_module("ccqelike_3d_me_fhc")
PRL_PT, PRL_PZ, SUMTP = _prl.PT, _prl.PZ, _prl.SUMTP

from xsec import compare, io  # noqa: E402
from xsec.extract import extract  # noqa: E402
from xsec.runner import run  # noqa: E402
from xsec.tunes import Tune  # noqa: E402

COARSE_PT = [0, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 2.5]
COARSE_PZ = [1.5, 3.5, 4.5]


def projector(sumtp, pt, pz):
    """Matrix mapping density bins (sumtp fastest) to integrated cross sections in coarse bins."""
    nx, ny, nz = len(sumtp) - 1, len(pt) - 1, len(pz) - 1
    cx, cy, cz = len(SUMTP) - 1, len(COARSE_PT) - 1, len(COARSE_PZ) - 1
    P = np.zeros((cx * cy * cz, nx * ny * nz))
    for k in range(nz):
        kc = np.searchsorted(COARSE_PZ, 0.5 * (pz[k] + pz[k + 1])) - 1
        if not 0 <= kc < cz:
            continue
        for j in range(ny):
            jc = np.searchsorted(COARSE_PT, 0.5 * (pt[j] + pt[j + 1])) - 1
            for i in range(nx):
                vol = (sumtp[i + 1] - sumtp[i]) * (pt[j + 1] - pt[j]) * (pz[k + 1] - pz[k])
                P[i + cx * (jc + cy * kc), i + nx * (j + ny * k)] = vol
    return P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", nargs="*", default=[])
    ap.add_argument("--data", nargs="*", default=[])
    ap.add_argument("-j", "--workers", type=int, default=16)
    ap.add_argument("-o", "--output", default=os.path.join(HERE, "results.pkl"))
    ap.add_argument("--load")
    args = ap.parse_args()
    meas = QELike3DLEBinning()
    if args.load:
        with open(args.load, "rb") as f:
            acc, info = pickle.load(f)
    else:
        acc, info = run(meas, io.df_files(*args.mc), io.df_files(*args.data), workers=args.workers)
        with open(args.output, "wb") as f:
            pickle.dump((dict(acc), {k: v for k, v in info.items() if k not in ("tune", "universes")}), f)
    info["tune"] = Tune(meas.tune_name, **meas.tune_options)
    info["universes"] = meas.universes(info["tune"])
    scales, _ = background_scales(acc, meas)
    res = extract(acc, meas, info, bkg_scales=scales)
    with open(args.output.replace(".pkl", "_xsec.pkl"), "wb") as f:
        pickle.dump(res, f)
    # project ours and PRL onto common coarse bins
    Po = projector(SUMTP, PT, PZ)
    Pr = projector(SUMTP, PRL_PT, PRL_PZ)
    ref, rc = prl_run.load_reference()
    xo, xr = Po @ res["xsec"], Pr @ ref
    co = Po @ res["cov_total"] @ Po.T
    cr = Pr @ rc["TotalCovariance"] @ Pr.T
    m = xr > 0
    print(compare.fmt(compare.summary(xo, xr, cr, co, m, "vs PRL (coarse bins)")))
    nx, ny = len(SUMTP) - 1, len(COARSE_PT) - 1
    for k in range(len(COARSE_PZ) - 1):
        sl = slice(k * nx * ny, (k + 1) * nx * ny)
        print(f"p|| [{COARSE_PZ[k]}, {COARSE_PZ[k + 1]}]: ratio {xo[sl].sum() / xr[sl].sum():.4f}")


if __name__ == "__main__":
    main()
