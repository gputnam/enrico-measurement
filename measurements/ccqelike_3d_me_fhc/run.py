#!/usr/bin/env python3
"""Run the ME FHC QE-like triple-differential measurement and compare with PRL 129, 021803.

  run.py --mc '<df dir/glob>' --data '<df dir/glob>' [-j 16] [-o results.pkl] [--no-fit]
  run.py --load results.pkl
"""
import argparse
import os
import pickle
import sys

import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analysis import PT, PZ, SUMTP, QELike3D, background_scales  # noqa: E402

from xsec import compare, io  # noqa: E402
from xsec.extract import extract  # noqa: E402
from xsec.runner import run  # noqa: E402
from xsec.tunes import Tune  # noqa: E402


def load_reference():
    f = uproot.open(os.path.join(HERE, "reference", "MINERvA_TripleDiffQELike_ptpzsumtp.root"))
    h = f["ptpzsumtp_data_cross_section_with_total_unc"]
    x = h.values().ravel(order="F")  # (sumtp, pt, pz) -> sumtp fastest
    covs = {}
    for name in ("TotalCovariance", "StatisticalCovariance", "FluxCovariance"):
        m = f[name]
        n = m.member("fNrows")
        full = np.array(m.member("fElements")).reshape(n, n)
        covs[name] = QELike3D.reco_binning.inner(QELike3D.reco_binning.inner(full, 0), 1)
    return x, covs


def report(res):
    ref, rc = load_reference()
    x = res["xsec"]
    mask = ref != 0
    print(f"flux integral {res['flux_integral'] * 1e-4:.4e} /cm^2/POT  targets {res['n_targets']:.4e}  POT {res['pot']:.4e}")
    print(compare.fmt(compare.summary(x, ref, rc["TotalCovariance"], res["cov_total"], mask, "total cov")))
    print(compare.fmt(compare.summary(x, ref, rc["StatisticalCovariance"], res["cov_stat"], mask, "stat cov")))
    for name, ours, theirs in (("stat", res["cov_stat"], rc["StatisticalCovariance"]),
                               ("flux", res["cov_bands"].get("Flux"), rc["FluxCovariance"]),
                               ("total", res["cov_total"], rc["TotalCovariance"])):
        if ours is None:
            continue
        fo = np.sqrt(np.clip(np.diag(ours), 0, None))[mask] / np.abs(x[mask])
        fr = np.sqrt(np.diag(theirs))[mask] / ref[mask]
        good = np.isfinite(fo) & (fo > 0)
        print(f"frac. unc. {name:<6} ours median {np.median(fo[good]):.4f}  ref median {np.median(fr):.4f}  "
              f"median ratio {np.median(fo[good] / fr[good]):.3f}")
    nx, ny, nz = len(SUMTP) - 1, len(PT) - 1, len(PZ) - 1
    X, Rf = x.reshape(nz, ny, nx), ref.reshape(nz, ny, nx)
    print("ratio integrated over SumTp and pT, per p|| bin:")
    for k in range(nz):
        print(f"  p|| [{PZ[k]:4.1f},{PZ[k + 1]:4.1f}]  {X[k].sum() / Rf[k].sum():.4f}")
    print("ratio integrated over p|| and pT, per SumTp bin:")
    for i in range(nx):
        print(f"  SumTp [{SUMTP[i]:.3f},{SUMTP[i + 1]:.3f}]  {X[:, :, i].sum() / Rf[:, :, i].sum():.4f}")
    print("ratio integrated over p|| and SumTp, per pT bin:")
    for j in range(ny):
        print(f"  pT [{PT[j]:.3f},{PT[j + 1]:.3f}]  {X[:, j].sum() / Rf[:, j].sum():.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", nargs="*", default=[])
    ap.add_argument("--data", nargs="*", default=[])
    ap.add_argument("-j", "--workers", type=int, default=16)
    ap.add_argument("-o", "--output", default=os.path.join(HERE, "results.pkl"))
    ap.add_argument("--load")
    ap.add_argument("--no-fit", action="store_true")
    ap.add_argument("--iter", type=int, default=None)
    args = ap.parse_args()
    meas = QELike3D()
    if args.load:
        with open(args.load, "rb") as f:
            acc, info = pickle.load(f)
        info["tune"] = Tune(meas.tune_name, **meas.tune_options)
        info["universes"] = meas.universes(info["tune"])
    else:
        acc, info = run(meas, io.df_files(*args.mc), io.df_files(*args.data), workers=args.workers)
        with open(args.output, "wb") as f:
            pickle.dump((dict(acc), {k: v for k, v in info.items() if k not in ("tune", "universes")}), f)
    scales = None
    if not args.no_fit:
        scales, coef = background_scales(acc, meas)
        cv = coef[:, 0, :]
        print("background fit (CV) scale factors, median over populated cells:",
              dict(zip(meas.fit_categories, np.round(np.median(cv[(cv != 1).any(1)], 0), 3))))
    res = extract(acc, meas, info, n_iter=args.iter, bkg_scales=scales)
    report(res)
    with open(args.output.replace(".pkl", "_xsec.pkl"), "wb") as f:
        pickle.dump(res, f)


if __name__ == "__main__":
    main()
