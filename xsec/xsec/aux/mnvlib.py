"""Load MINERvA MnvH1D/MnvVertErrorBand classes for PyROOT without the MINERvA software.

The classes are generated from the files' own streamer info with ROOT's TFile::MakeProject,
so reading needs nothing beyond ROOT itself.
"""
import os

import ROOT

_LIB = None


def load(example_file, workdir):
    """Build (once) and load the generated class library."""
    global _LIB
    if _LIB:
        return
    so = os.path.join(workdir, "MnvLib", "MnvLib.so")
    if not os.path.exists(so):
        os.makedirs(workdir, exist_ok=True)
        cwd = os.getcwd()
        os.chdir(workdir)
        f = ROOT.TFile.Open(example_file)
        f.MakeProject("MnvLib", "*", "recreate++")
        f.Close()
        os.chdir(cwd)
    ROOT.gSystem.Load(so)
    _LIB = so


def th1_to_arrays(h):
    """Bin edges and contents (without under/overflow) of a TH1."""
    import numpy as np
    n = h.GetNbinsX()
    edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, n + 2)])
    vals = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
    return edges, vals


def mnvh1d_to_dict(h):
    """CV, edges, and all vertical error bands {name: (univ contents [N, nbins], weights)}."""
    import numpy as np
    edges, cv = th1_to_arrays(h)
    out = {"edges": edges, "cv": cv}
    for kv in h.fVertErrorBandMap:
        name, eb = str(kv.first), kv.second
        univ = np.array([th1_to_arrays(eb.fHists[i])[1] for i in range(eb.fHists.size())])
        w = np.array(list(eb.fUnivWgts)) if eb.fUnivWgts.size() else np.ones(len(univ))
        out[f"band_{name}"] = univ
        out[f"bandw_{name}"] = w
    return out
