#!/usr/bin/env python3
"""Run the ME FHC CC-inclusive (pT, p||) measurement and compare with PRD 104, 092007.

  run.py --mc '<df dir or glob>' --data '<df dir or glob>' [-j 16] [-o results.pkl]
  run.py --load results.pkl            (re-extract / compare only)
"""
import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analysis import PT, PZ, CCInclusive2D  # noqa: E402

from xsec import compare, io  # noqa: E402
from xsec.extract import extract  # noqa: E402
from xsec.runner import run  # noqa: E402


def load_reference():
    ref = pd.read_csv(os.path.join(HERE, "reference", "data_result_ptpl_2D_minerva_inclusive_6GeV.txt"))
    ref.columns = ["pzbin", "ptbin", "xsec", "stat", "total"]
    npz, npt = len(PZ) - 1, len(PT) - 1
    # our order: pt fastest (Binning([PT, PZ]) inner bins)
    order = np.array([(iz - 1) * npt + (it - 1) for iz, it in zip(ref.pzbin, ref.ptbin)])
    gid = np.array([(iz - 1) + npz * (it - 1) for iz, it in zip(ref.pzbin, ref.ptbin)])
    x = np.zeros(npt * npz)
    x[order] = ref.xsec
    covs = {}
    for name in ("total", "stat", "flux", "energyscale"):
        c = pd.read_csv(os.path.join(HERE, "reference", f"cov_ptpl_minerva_inclusive_6GeV_{name}.txt"))
        c.columns = ["i", "j", "v"]
        M = np.zeros((npt * npz, npt * npz))
        g2o = dict(zip(gid, order))
        M[[g2o[i] for i in c.i], [g2o[j] for j in c.j]] = c.v
        covs[name] = M
    return x, covs


def report(res, meas):
    ref, rcov = load_reference()
    x = res["xsec"]
    mask = ref != 0
    print(f"flux integral {res['flux_integral'] * 1e-4:.4e} /cm^2/POT   targets {res['n_targets']:.4e}   POT {res['pot']:.4e}")
    s = compare.summary(x, ref, rcov["total"], res["cov_total"], mask, "total cov")
    print(compare.fmt(s))
    s = compare.summary(x, ref, rcov["stat"], res["cov_stat"], mask, "stat cov")
    print(compare.fmt(s))
    # uncertainty comparison
    for name, ours in (("stat", res["cov_stat"]), ("flux", res["cov_bands"].get("Flux")), ("total", res["cov_total"])):
        if ours is None:
            continue
        good = mask & (x > 0)
        fr_ours = np.sqrt(np.clip(np.diag(ours), 0, None))[good] / x[good]
        fr_ref = np.sqrt(np.diag(rcov[name]))[good] / ref[good]
        print(f"frac. unc. {name:<6} ours median {np.median(fr_ours):.4f}  ref median {np.median(fr_ref):.4f}  "
              f"median ratio {np.median(fr_ours / fr_ref):.3f}")
    npt = len(PT) - 1
    print("per-p|| bin ratio (sum over pT, reported bins):")
    for iz in range(len(PZ) - 1):
        sl = slice(iz * npt, (iz + 1) * npt)
        m = mask[sl]
        if m.any():
            print(f"  p|| [{PZ[iz]:5.1f},{PZ[iz + 1]:5.1f}]  ratio {x[sl][m].sum() / ref[sl][m].sum():.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", nargs="*", default=[])
    ap.add_argument("--data", nargs="*", default=[])
    ap.add_argument("-j", "--workers", type=int, default=16)
    ap.add_argument("-o", "--output", default=os.path.join(HERE, "results.pkl"))
    ap.add_argument("--load")
    ap.add_argument("--mat-signal", action="store_true", help="MINERvA-101 efficiency-numerator convention")
    ap.add_argument("--iter", type=int, default=None)
    args = ap.parse_args()
    meas = CCInclusive2D(mat_signal_convention=args.mat_signal)
    if args.load:
        with open(args.load, "rb") as f:
            acc, info = pickle.load(f)
        from xsec.tunes import Tune
        info["tune"] = Tune(meas.tune_name)
        info["universes"] = meas.universes(info["tune"])
    else:
        acc, info = run(meas, io.df_files(*args.mc), io.df_files(*args.data), workers=args.workers)
        with open(args.output, "wb") as f:
            pickle.dump((dict(acc), {k: v for k, v in info.items() if k not in ("tune", "universes")}), f)
    res = extract(acc, meas, info, n_iter=args.iter)
    report(res, meas)
    with open(args.output.replace(".pkl", "_xsec.pkl"), "wb") as f:
        pickle.dump({k: v for k, v in res.items()}, f)


if __name__ == "__main__":
    main()
