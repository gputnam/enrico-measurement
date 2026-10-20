#!/usr/bin/env python3
"""Generic driver for a corpus measurement: event loop on the DFs, optional background fit, unfolding.

Runs every entry of measurements/<m>/corpus.json (one Measurement instance per entry; entries sharing
class+kwargs share one event loop) and writes, for entry e,
    <results-root>/<e.result>            extracted cross section (*_xsec.pkl, input of `xsec.corpus add`)
    <results-root>/<e.result minus _xsec>.pkl   accumulators (re-extract with --load)
If analysis.py defines `background_scales(acc, meas) -> (scales, coef)` it is used as the background fit.
Fiducial volume for the target count: Measurement attribute `fv = (zmin, zmax, apothem)` (mm), default
tracker FV (5980, 8422, 850).

  run_measurement.py <m> --results-root R [--mc DIR] [--data DIR] [-j 16] [--max-files N] [--load] [--ids ...]

--max-files N uses only the first N MC and N data DF files per playlist (fast development; the
normalisation is still correct since POT is counted per file, but statistics are poor).
Closure: the printed closure ratio (extract with data := MC prediction, compared with the truth
prediction) must be ~1 in every bin with sizeable efficiency.
"""
import argparse
import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from xsec import io  # noqa: E402
from xsec.extract import extract  # noqa: E402
from xsec.runner import run  # noqa: E402
from xsec.tunes import Tune  # noqa: E402
from xsec.util import load_measurement_module  # noqa: E402

DF_MC = "/pnfs/sbnd/scratch/users/gputnam/MvA-Loop/df/v1/me_fhc_mc"
DF_DATA = "/pnfs/sbnd/scratch/users/gputnam/MvA-Loop/df/v1/me_fhc_data"


def limit_files(paths, n):
    if not n:
        return paths
    by = {}
    for p in paths:
        by.setdefault(io._playlist_from_path(p) or os.path.basename(p).split("_")[3], []).append(p)
    return sorted(q for v in by.values() for q in sorted(v)[:n])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("measurement")
    ap.add_argument("--results-root", required=True)
    ap.add_argument("--mc", default=os.environ.get("MVA_DF_MC", DF_MC))
    ap.add_argument("--data", default=os.environ.get("MVA_DF_DATA", DF_DATA))
    ap.add_argument("-j", "--workers", type=int, default=16)
    ap.add_argument("--max-files", type=int, default=0)
    ap.add_argument("--load", action="store_true", help="reuse saved accumulators")
    ap.add_argument("--ids", nargs="*")
    args = ap.parse_args()
    mod = load_measurement_module(args.measurement)
    entries = json.load(open(os.path.join(HERE, args.measurement, "corpus.json")))
    if args.ids:
        entries = [e for e in entries if e["id"] in args.ids]
    mc = limit_files(io.df_files(args.mc), args.max_files)
    data = limit_files(io.df_files(args.data), args.max_files)
    for e in entries:
        meas = getattr(mod, e["class"])(**e.get("kwargs", {}))
        out = os.path.join(args.results_root, e["result"])
        acc_path = out.replace("_xsec.pkl", ".pkl")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        tune = Tune(meas.tune_name, **getattr(meas, "tune_options", {}))
        if args.load:
            with open(acc_path, "rb") as f:
                acc, info = pickle.load(f)
            info["tune"], info["universes"] = tune, meas.universes(tune)
        else:
            print(f"== {e['id']}: {len(mc)} MC, {len(data)} data DF files", flush=True)
            acc, info = run(meas, mc, data, workers=args.workers)
            with open(acc_path, "wb") as f:
                pickle.dump((dict(acc), {k: v for k, v in info.items() if k not in ("tune", "universes")}), f)
        scales = None
        if hasattr(mod, "background_scales"):
            scales, _ = mod.background_scales(acc, meas)
        zmin, zmax, apothem = getattr(meas, "fv", (5980.0, 8422.0, 850.0))
        kw = dict(n_iter=None, bkg_scales=scales, zmin=zmin, zmax=zmax, apothem=apothem)
        res = extract(acc, meas, info, **kw)
        clo = extract(acc, meas, info, closure=True, **kw)
        with open(out, "wb") as f:
            pickle.dump(res, f)
        x, t, s = res["xsec"], res["truth_xsec"], np.sqrt(np.clip(np.diag(res["cov_total"]), 0, None))
        good = (clo["truth_xsec"] > 0) & (clo["eff"] > 0.02)
        cr = clo["xsec"][good] / clo["truth_xsec"][good]
        vol = meas.true_binning.widths()
        print(f"{e['id']}: {len(x)} bins, data/GENIE integral {(x * vol).sum() / (t * vol).sum():.3f}, "
              f"median frac. unc. {np.median(s[x > 0] / x[x > 0]):.3f}, "
              f"closure ratio median {np.median(cr):.4f} min {cr.min():.3f} max {cr.max():.3f} "
              f"(bins with eff > 2%: {good.sum()}), efficiency median {np.median(res['eff']):.3f}")
        print(f"   wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
