#!/usr/bin/env python3
"""Check makedf grid outputs against a job list and write a resubmission list.

  grid/check_outputs.py <joblist> <outdir> [--resubmit-list out.txt]

A job is OK when <outdir>/<outname>.root exists, opens, and its meta tree has one row per
input file of the job with POT > 0.  Prints total POT, rows and size.
"""
import argparse
import os
import sys

import uproot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("joblist")
    ap.add_argument("outdir")
    ap.add_argument("--resubmit-list")
    args = ap.parse_args()
    bad, pot, nreco, ntruth, size = [], 0.0, 0, 0, 0
    lines = [l.split() for l in open(args.joblist) if l.strip()]
    for name, inputs in lines:
        path = os.path.join(args.outdir, name + ".root")
        try:
            with uproot.open(path) as f:
                m = f["meta"].arrays(["pot_used", "n_reco_pass", "n_truth_pass"], library="np")
                ok = len(m["pot_used"]) == len(inputs.split(",")) and (m["pot_used"] > 0).all()
                if ok:
                    pot += m["pot_used"].sum()
                    nreco += m["n_reco_pass"].sum()
                    ntruth += m["n_truth_pass"].sum()
                    size += os.path.getsize(path)
        except Exception:
            ok = False
        if not ok:
            bad.append(f"{name} {inputs}")
    print(f"{len(lines) - len(bad)}/{len(lines)} OK   POT {pot:.4e}   reco rows {nreco}   truth rows {ntruth}   "
          f"size {size / 1e9:.2f} GB")
    if args.resubmit_list and bad:
        with open(args.resubmit_list, "w") as f:
            f.write("\n".join(bad) + "\n")
        print(f"wrote {len(bad)} jobs to {args.resubmit_list}")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
