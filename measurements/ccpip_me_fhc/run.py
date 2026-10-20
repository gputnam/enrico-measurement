#!/usr/bin/env python3
"""Run the ME CC1pi+ measurements (all variables) and compare with arXiv:2605.24224.

  run.py --mc '<df dir>' --data '<df dir>' [-j 16] [--vars tpi pmu ...] [-o outdir]
  run.py --load-dir outdir
"""
import argparse
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analysis import VARIABLES, CCPip, background_scales  # noqa: E402
from reference_tables import load_cov, load_xsec  # noqa: E402

from xsec import compare, io  # noqa: E402
from xsec.extract import extract  # noqa: E402
from xsec.runner import run  # noqa: E402
from xsec.tunes import Tune  # noqa: E402


def report(var, res, ref_x, ref_c):
    meta = VARIABLES[var]
    x = res["xsec"] * 1e42
    r = ref_x[:, 2]
    mask = r > 0
    fstat = np.sqrt(np.clip(np.diag(res["cov_stat"]), 0, None)) * 1e42 / np.where(x != 0, np.abs(x), 1)
    fsys = np.sqrt(np.clip(np.diag(res["cov_sys"]), 0, None)) * 1e42 / np.where(x != 0, np.abs(x), 1)
    print(f"--- {var}: integral ratio {x[mask].sum() / r[mask].sum():.4f}")
    if meta["cov"] and (meta["cov"], "stat") in ref_c:  # tables in units of 1e-84 = (1e-42)^2
        cov = ref_c[(meta["cov"], "stat")] + ref_c[(meta["cov"], "sys")]
        ours = (res["cov_total"]) * 1e84
        c2, n = compare.chi2(x, r, cov + ours, mask)
        print(f"  chi2(ref+ours cov) = {c2:.1f} / {n}")
    print("  bin        ours      ref    ratio  fstat(ours/ref)  fsys(ours/ref)")
    for i in range(len(r)):
        rt = x[i] / r[i] if r[i] > 0 else np.nan
        print(f"  [{ref_x[i, 0]:g},{ref_x[i, 1]:g}]  {x[i]:9.4g} {r[i]:9.4g}  {rt:6.3f}   {fstat[i]:.3f}/{ref_x[i, 4]:.3f}"
              f"     {fsys[i]:.3f}/{ref_x[i, 5]:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", nargs="*", default=[])
    ap.add_argument("--data", nargs="*", default=[])
    ap.add_argument("-j", "--workers", type=int, default=16)
    ap.add_argument("--vars", nargs="*", default=list(VARIABLES))
    ap.add_argument("-o", "--outdir", default=os.path.join(HERE, "results"))
    ap.add_argument("--load-dir")
    ap.add_argument("--no-fit", action="store_true")
    args = ap.parse_args()
    ref_x, ref_c = load_xsec(), load_cov()
    os.makedirs(args.outdir, exist_ok=True)
    for var in args.vars:
        meas = CCPip(var)
        path = os.path.join(args.load_dir or args.outdir, f"{var}.pkl")
        if args.load_dir:
            with open(path, "rb") as f:
                acc, info = pickle.load(f)
        else:
            acc, info = run(meas, io.df_files(*args.mc), io.df_files(*args.data), workers=args.workers, verbose=False)
            with open(path, "wb") as f:
                pickle.dump((dict(acc), {k: v for k, v in info.items() if k not in ("tune", "universes")}), f)
        info["tune"] = Tune(meas.tune_name, **meas.tune_options)
        info["universes"] = meas.universes(info["tune"])
        scales = None
        if not args.no_fit:
            scales, sc = background_scales(acc, meas)
            print(f"{var}: W sideband fit scales (CV) wmid {sc[0, 0]:.3f}  whigh {sc[1, 0]:.3f}")
        res = extract(acc, meas, info, bkg_scales=scales)
        with open(os.path.join(args.outdir, f"{var}_xsec.pkl"), "wb") as f:
            pickle.dump(res, f)
        report(var, res, ref_x[VARIABLES[var]["table"]], ref_c)


if __name__ == "__main__":
    main()
