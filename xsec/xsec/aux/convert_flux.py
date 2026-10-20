#!/usr/bin/env python3
"""Convert MINERvA flux MnvH1D files to npz.

Usage (spack ROOT env):  python3 -m xsec.aux.convert_flux <aux_dir>
Writes <aux_dir>/npz/flux/<file stem>__<hist>.npz with keys edges, cv, band_<name>, bandw_<name>.
"""
import glob
import os
import sys

import numpy as np
import ROOT

from . import mnvlib


def main(aux):
    fluxdir = os.path.join(aux, "MATFluxAndReweightFiles", "flux")
    outdir = os.path.join(aux, "npz", "flux")
    os.makedirs(outdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(fluxdir, "flux-*.root")) + glob.glob(os.path.join(fluxdir, "sys", "*.root")))
    mnvlib.load(files[0], os.path.join(aux, "mnvproj"))
    for fn in files:
        f = ROOT.TFile.Open(fn)
        for key in f.GetListOfKeys():
            obj = key.ReadObj()
            if not obj.InheritsFrom("TH1") or obj.InheritsFrom("TH2"):
                continue
            if obj.ClassName() == "PlotUtils::MnvH1D":
                d = mnvlib.mnvh1d_to_dict(obj)
            else:
                e, v = mnvlib.th1_to_arrays(obj)
                d = {"edges": e, "cv": v}
            stem = os.path.basename(fn)[:-5]
            out = os.path.join(outdir, f"{stem}__{key.GetName()}.npz")
            np.savez_compressed(out, **d)
            print(out, {k: np.shape(v) for k, v in d.items()})
        f.Close()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/exp/sbnd/data/users/gputnam/MvA-Loop/aux")
