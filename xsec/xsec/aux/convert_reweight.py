#!/usr/bin/env python3
"""Convert MINERvA reweight ROOT files (TH1/TH2/TGraph/TArrayD objects) to npz.

Usage (spack ROOT env): python3 -m xsec.aux.convert_reweight [aux_dir]
For each object <name> in file <stem>.root writes <aux>/npz/reweight/<stem>.npz with keys
  <name>__x, <name>__y            : TH1 / TH2 bin edges
  <name>__c                       : contents INCLUDING under/overflow, indexed like GetBinContent
  <name>__gx, <name>__gy          : TGraph points
  <name>__arr                     : TArrayD
"""
import glob
import os
import sys

import numpy as np
import ROOT

from . import mnvlib

FILES = [
    "MParamFiles/data/Reweight/outNievesRPAratio-*.root",
    "MParamFiles/data/Reweight/lowQ2pi_weights.root",
    "MParamFiles/data/Reweight/cohPi_weights_tracker_ME.root",
    "MParamFiles/data/Reweight/Z_Expansion_Reweight_v2126.root",
    "MParamFiles/data/Reweight/fsi_weight_qe_flat_1d50.root",
]


def edges(ax):
    return np.array([ax.GetBinLowEdge(i) for i in range(1, ax.GetNbins() + 2)])


def dump(obj, name, out):
    if not hasattr(obj, "InheritsFrom"):  # bare TArrayD (TH1D also derives from TArrayD)
        out[name + "__arr"] = np.array([obj.At(i) for i in range(obj.GetSize())])
    elif obj.ClassName() == "PlotUtils::MnvH1D":
        d = mnvlib.mnvh1d_to_dict(obj)
        out[name + "__x"] = d["edges"]
        c = np.zeros(len(d["cv"]) + 2)
        c[1:-1] = d["cv"]
        out[name + "__c"] = c
        for k, v in d.items():
            if k.startswith("band"):
                out[f"{name}__{k}"] = v
    elif obj.InheritsFrom("TH2"):
        nx, ny = obj.GetNbinsX(), obj.GetNbinsY()
        out[name + "__x"] = edges(obj.GetXaxis())
        out[name + "__y"] = edges(obj.GetYaxis())
        c = np.zeros((nx + 2, ny + 2))
        for i in range(nx + 2):
            for j in range(ny + 2):
                c[i, j] = obj.GetBinContent(i, j)
        out[name + "__c"] = c
    elif obj.InheritsFrom("TH1"):
        n = obj.GetNbinsX()
        out[name + "__x"] = edges(obj.GetXaxis())
        out[name + "__c"] = np.array([obj.GetBinContent(i) for i in range(n + 2)])
    elif obj.InheritsFrom("TGraph"):
        n = obj.GetN()
        out[name + "__gx"] = np.array([obj.GetX()[i] for i in range(n)])
        out[name + "__gy"] = np.array([obj.GetY()[i] for i in range(n)])


def main(aux):
    outdir = os.path.join(aux, "npz", "reweight")
    os.makedirs(outdir, exist_ok=True)
    mnvlib.load(os.path.join(aux, "MATFluxAndReweightFiles", "flux", "flux-gen2thin-pdg14-minervame1D_rearrangedUniverses.root"),
                os.path.join(aux, "mnvproj"))
    for pat in FILES:
        for fn in sorted(glob.glob(os.path.join(aux, pat))):
            f = ROOT.TFile.Open(fn)
            out = {}
            for key in f.GetListOfKeys():
                if key.GetClassName() == "TArrayD":
                    dump(f.Get(key.GetName()), key.GetName(), out)
                    continue
                obj = key.ReadObj()
                if obj:
                    dump(obj, key.GetName(), out)
            dst = os.path.join(outdir, os.path.basename(fn)[:-5] + ".npz")
            np.savez_compressed(dst, **out)
            print(dst, sorted(out.keys())[:12], "..." if len(out) > 12 else "")
            f.Close()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/exp/sbnd/data/users/gputnam/MvA-Loop/aux")
