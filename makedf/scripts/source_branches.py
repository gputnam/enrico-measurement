#!/usr/bin/env python3
"""Write makedf/config/SOURCE_BRANCHES.md: the source AnaTuple branches that are NOT in the DFs.

Input: a leaf listing of one MC AnaTuple (and optionally one data AnaTuple), one leaf per line
    <tree> <leaf> <type> <static length> <count leaf or ->
made with makedf/scripts/list_leaves.C, e.g.
    root -l -b -q 'makedf/scripts/list_leaves.C("root://.../MasterAnaDev_mc_AnaTuple_run00110000_Playlist.root")' > mc.txt
Usage: makedf/scripts/source_branches.py mc.txt [data.txt] > makedf/config/SOURCE_BRANCHES.md
"""
import collections
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MAKEDF = os.path.dirname(HERE)

# (category, regex on the branch name), first match wins
CATEGORIES = [
    ("GENIE / reweight knobs and MC weights (MC only)", r"^(truth_genie_wgt_|mc_wgt|mc_cvweight|mc_ppfx|mc_hornCurrent|mc_gen|truth_.*_wgt|wgt_|.*_weight)"),
    ("MC event record and truth (mc_*, truth tree quantities)", r"^(mc_|truth_|genie_|FittedMichel_.*true|improved_michel_matched|has_michel_primary|truth)"),
    ("Particle-response / recoil-energy systematic breakdowns (part_response_*)", r"^part_response_"),
    ("Per-cluster / per-hit energy points and systematics (EnergyPoint*, cluster_*, Signal*, odlattice*)", r"^(EnergyPoint|cluster_|Signal|n_Signal|odlattice|n_odindices|slice_|n_slices|clusters_found)"),
    ("Recoil blobs by particle type and subdetector (blob_*_E_<det>_<particle>, MC truth split)", r"^blob_"),
    ("pi0 / photon / shower reconstruction", r"^(pi0_|gamma[12]_|shower_|anglescan|Good_Blobs|Cut_.*Pi0|is_Pi0|ExtraEnergyClusters|Cut_Blob)"),
    ("NuE-style prong / EM-shower ID (filled only for non-MINOS-matched events: -999 in the DF preselection)", r"^(prong_|n_prongs|ConeEnergy|ExtraEnergy|UpstreamInline|EMLike|VertexTrackMultiplicity|StartPointVertexMultiplicity|FailContainedProng|iso_prong|primary_prong|n_iso_trk_prongs)"),
    ("Michel electrons", r"(michel|Michel)"),
    ("Isolated / vertex / muon blobs", r"(iso_blobs|vtx_blobs|muon_iso|nonvtx_iso|vtx_iso)"),
    ("Proton candidates", r"(proton|Proton)"),
    ("Pion / hadron candidates", r"(pion|Pion|hadron|piFit|_pi_)"),
    ("Muon and MINOS", r"(muon|minos|MINOS|Minos|_mu_)"),
    ("Blob tracks (short tracks built from blobs, MasterAnaDev_BlobTrack*)", r"BlobTrack|EvtHasNBlobTracks|MCEnergyFrac"),
    ("ANN / ML vertexing and nuclear targets", r"(ANN|ML|target|Target|nuclTargs)"),
    ("Recoil / calorimetry", r"(recoil|Recoil|energy|Energy|_E$|_E_|Esum|calE|visE)"),
    ("Vertex, tracks and event topology", r"(vtx|vertex|Vertex|track|Track|trk|prong|n_|has_|is_|phys_)"),
    ("Beam, run and bookkeeping", r"^(ev_|numi_|batch|gate|run|subrun|n_hyps|physEvtNum|dead)"),
]


NOTES = {
    "pi0 / photon / shower reconstruction":
        "Checked on ME FHC MC: `*_Old`, `*_hcal`, `gamma?_end_vertex`, `Cut_BlobDirectionBad` and "
        "`shower_dist_fitline` are never filled; `pi0_E_Cal` == `pi0_E`, `gamma?_P` == `gamma?_E`.  "
        "Everything else of the pi0 reconstruction is in the DF (see measurements/AUTHORING.md).",
    "NuE-style prong / EM-shower ID (filled only for non-MINOS-matched events: -999 in the DF preselection)":
        "These are -999 / -1 for every event passing the DF preselection (MINOS-matched muon), so they "
        "cannot help a numu CC measurement even if added.",
    "Particle-response / recoil-energy systematic breakdowns (part_response_*)":
        "Recoil energy split by subdetector / particle category with response uncertainties (the DF "
        "keeps the main ones as pr_*, rerr_*, rpass_*).",
}


def used_branches():
    used = set()
    for f in glob.glob(os.path.join(MAKEDF, "config", "*.cols")):
        for line in open(f):
            line = line.split("#")[0].split()
            if len(line) >= 2:
                used.add(line[1].split("[")[0])
    for f in glob.glob(os.path.join(MAKEDF, "src", "*.cxx")):
        used |= set(re.findall(r'(?:require|request)\("([A-Za-z0-9_]+)"\)', open(f).read()))
    return used


def read(path):
    out = collections.OrderedDict()
    for line in open(path):
        t = line.split()
        if len(t) < 3 or t[0] not in ("MasterAnaDev", "Truth"):
            continue
        out.setdefault(t[1], (t[0], t[2], t[3] if len(t) > 3 else "1", t[4] if len(t) > 4 else "-"))
    return out


def compress(names):
    """Collapse names sharing a prefix (up to the last '_') into prefix_{a,b,...}."""
    groups = collections.OrderedDict()
    for n in sorted(names):
        p, _, s = n.rpartition("_")
        groups.setdefault(p if p else n, []).append(s if p else "")
    out = []
    for p, ss in groups.items():
        if len(ss) >= 3 and all(ss):
            out.append(f"{p}_{{{','.join(ss)}}}")
        else:
            out += [f"{p}_{s}" if s else p for s in ss]
    return out


def main():
    mc = read(sys.argv[1])
    data = read(sys.argv[2]) if len(sys.argv) > 2 else None
    used = used_branches()
    sz = {n for n in mc if n.endswith("_sz") and n[:-3] in mc}
    todo = [n for n in mc if n not in used and n not in sz and not (n.endswith("_sz") and n[:-3] in used)]
    cats = collections.OrderedDict((c, []) for c, _ in CATEGORIES)
    cats["Other"] = []
    for n in todo:
        for c, rx in CATEGORIES:
            if re.search(rx, n):
                cats[c].append(n)
                break
        else:
            cats["Other"].append(n)
    w = sys.stdout.write
    w("# Source AnaTuple branches NOT in the DFs\n\n")
    w("Generated by `makedf/scripts/source_branches.py` from the leaf list of one ME FHC MC AnaTuple\n"
      "(`MasterAnaDev` and `Truth` trees; data AnaTuples have the same reco branches without the MC ones).\n"
      f"{len(mc)} source leaves, {len(used)} read by makedf (column files + derived C++ code), "
      f"{len(todo)} not in the DFs (variable-length arrays' `_sz` counters omitted).\n\n")
    w("Use this to tell whether a measurement is blocked by the **DF schema** (the quantity exists below and\n"
      "could be added to `makedf/config/*.cols` with a DF reprocessing) or by the **data** (no such branch\n"
      "anywhere).  `{a,b}` is brace expansion.  `[N]` = fixed array length, `[cnt]` = variable length.\n"
      "Branches marked `(MC)` are absent from data AnaTuples.  Many reco branches are only filled for\n"
      "some event classes; check a branch's content before relying on it.\n\n")
    for c, names in cats.items():
        if not names:
            continue
        w(f"## {c} ({len(names)})\n\n")
        # annotate arrays and MC-only
        def ann(n):
            tree, typ, ln, cnt = mc[n]
            a = ""
            if cnt != "-":
                a = "[cnt]"
            elif ln not in ("1", "0"):
                a = f"[{ln}]"
            if tree == "Truth":
                a += " (Truth tree)"
            elif data is not None and n not in data:
                a += " (MC)"
            return a
        if c in NOTES:
            w(NOTES[c] + "\n\n")
        if len(names) > 300:  # huge systematic families: summarise by prefix
            fam = collections.Counter("_".join(n.split("_")[:4]) for n in names)
            w("- prefixes (count): " + ", ".join(f"`{k}_*` ({v})" for k, v in sorted(fam.items())) + "\n\n")
            continue
        by_ann = collections.OrderedDict()
        for n in names:
            by_ann.setdefault(ann(n), []).append(n)
        for a, ns in by_ann.items():
            label = a.strip() or "scalar"
            w(f"- {label}: " + ", ".join(f"`{x}`" for x in compress(ns)) + "\n")
        w("\n")


if __name__ == "__main__":
    main()
