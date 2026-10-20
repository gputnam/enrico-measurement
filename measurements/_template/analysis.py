"""TEMPLATE: copy this directory to measurements/<new_name>/ and edit.

Example: numu CC inclusive dsigma/dp_mu for theta_mu < 20 deg (MnvTune v1), selection as in
ccinc_2d_me_fhc.  See ../AUTHORING.md for the interface, the available DF columns and the rules.
Units in the event accessors: MeV, mm, rad.  Report observables in GeV / GeV/c / deg / GeV^2 and
say so in definition.md and in the corpus.json title.
"""
import numpy as np

from xsec.binning import Binning
from xsec.measurement import Measurement

ZMIN, ZMAX, APOTHEM = 5980.0, 8422.0, 850.0
PMU = [1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 15, 20]  # GeV/c


class TemplateMeasurement(Measurement):
    name = "_template"
    tune_name = "v1"          # MnvTune used for the efficiency/background model and the GENIE prediction
    n_iter = 4                 # D'Agostini iterations
    fv = (ZMIN, ZMAX, APOTHEM)  # fiducial volume for the target count
    samples = ("sel",)
    signal_sample = "sel"
    bkg_categories = ("NC", "wrong_sign", "nue", "out_of_phase_space")
    reco_binning = Binning([PMU], ["pmu"])
    true_binning = Binning([PMU], ["pmu"])
    extra_reco_columns = []    # reco-tree columns beyond MUON_COLUMNS (xsec/xsec/measurement.py)
    extra_truth_columns = []   # truth columns beyond TRUTH_COLUMNS (e.g. "tfs_np", "tfs_sumTp")

    def preselect(self, ev):
        """Cheap, shift-independent superset of the selection (reco rows)."""
        return ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (ev["minos_match"] == 1)

    def reco(self, ev, shift):
        """(tuple of reco observables, {sample: mask}) in the lateral systematic `shift`."""
        mu = ev.muon(shift)
        sel = (ev.reco_vertex_in_fv(ZMIN, ZMAX, APOTHEM) & (mu["theta"] < np.radians(20.0))
               & (ev["minos_match"] == 1) & (ev["dead"] <= 1) & (ev["minos_qp"] < 0))
        return (mu["p"] / 1000.0,), {"sel": sel}

    def truth(self, ev):
        """True observables, same order/units as reco().  Must only use generator-available truth
        columns (see AUTHORING.md) so the model prediction can be made."""
        t = ev.true_muon()
        return (t["p"] / 1000.0,)

    def signal(self, ev):
        """Signal definition (incl. phase space) on truth columns."""
        t = ev.true_muon()
        return ((ev["mc_incoming"] == 14) & (ev["mc_current"] == 1) & ev.true_vertex_in_fv(ZMIN, ZMAX, APOTHEM)
                & (t["theta"] < np.radians(20.0)))

    def category(self, ev):
        """Background category index (into bkg_categories) for non-signal selected events."""
        cat = np.full(len(ev), 3)
        cat[ev["mc_current"] != 1] = 0
        cc = ev["mc_current"] == 1
        cat[cc & (ev["mc_incoming"] == -14)] = 1
        cat[cc & (np.abs(ev["mc_incoming"]) == 12)] = 2
        return cat

# Optional background fit (see ccqelike_3d_me_fhc/analysis.py, ccpip_me_fhc/analysis.py):
# def background_scales(acc, meas):
#     return {category: (R, U) scale array}, coefficients
