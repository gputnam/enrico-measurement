#!/usr/bin/env python3
"""Build makedf job lists from the MINERvA Open Data directory tree.

Each output line is one grid job:  <outname> <url1>[,<url2>,...]
MC: one input file per job.  Data: files grouped up to --data-gb per job.

Example:
  grid/make_joblists.py --sample MediumEnergy_FHC --kind MC --playlists 1A,1D -o joblists/v1/me_fhc_mc.txt
"""
import argparse
import os
import subprocess
import sys

DOOR = "root://fndcadoor.fnal.gov:1095"
BASE = "/pnfs/fnal.gov/usr/minerva/persistent/OpenData"


def xrd_ls(path):
    out = subprocess.run(["xrdfs", DOOR, "ls", "-l", path], capture_output=True, text=True, check=True).stdout
    entries = []
    for line in out.splitlines():
        parts = line.split()
        if len(parts) < 5:
            continue
        entries.append((parts[-1], int(parts[-2]), parts[0].startswith("d")))
    return entries


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True, help="e.g. MediumEnergy_FHC, MediumEnergy_RHC, LowEnergy_FHC")
    ap.add_argument("--kind", required=True, choices=["MC", "Data"])
    ap.add_argument("--subsample", default="StandardMC", help="MC subsample directory")
    ap.add_argument("--playlists", default="", help="comma separated (e.g. 1A,1D); default all")
    ap.add_argument("--data-gb", type=float, default=5.0, help="max input GB per data job")
    ap.add_argument("--max-files", type=int, default=0, help="per playlist file cap (for small-scale tests)")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    top = f"{BASE}/{args.sample}/{args.kind}" + (f"/{args.subsample}" if args.kind == "MC" else "")
    pls = [os.path.basename(p) for p, _, d in xrd_ls(top) if d and os.path.basename(p).startswith("Playlist")]
    if args.playlists:
        want = {"Playlist" + p for p in args.playlists.split(",")}
        pls = [p for p in pls if p in want]
    if not pls:
        sys.exit(f"no playlists found under {top}")

    short = args.sample.replace("MediumEnergy", "ME").replace("LowEnergy", "LE").lower()
    lines = []
    for pl in sorted(pls):
        files = sorted((p, s) for p, s, d in xrd_ls(f"{top}/{pl}") if not d and p.endswith(".root"))
        if args.max_files:
            files = files[: args.max_files]
        tag = pl.replace("Playlist", "")
        if args.kind == "MC":
            for p, _ in files:
                run = os.path.basename(p).split("_run")[1].split("_")[0]
                lines.append(f"{short}_mc_{tag}_{run} {DOOR}/{p}")
        else:
            group, size, ig = [], 0, 0
            for p, s in files + [(None, 0)]:
                if p is None or (group and size + s > args.data_gb * 1e9):
                    if group:
                        lines.append(f"{short}_data_{tag}_{ig:04d} " + ",".join(f"{DOOR}/{x}" for x in group))
                        ig += 1
                    group, size = [], 0
                if p is not None:
                    group.append(p)
                    size += s
        print(f"{pl}: {len(files)} files", file=sys.stderr)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {len(lines)} jobs to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
