"""QE-like triple-differential cross section in the binning of the LE/ME comparison paper.

Reference: D. Ruterbories et al., arXiv:2606.00745 (LE vs ME QE-like 3D).  Same selection and
method as PRL 129, 021803 (see ../ccqelike_3d_me_fhc), with p|| restricted to 1.5-4.5 GeV/c and
the background fit integrated over that range.

Binning: p|| and pT edges are from the paper text / figure panel labels; the SumTp edges are not
given in the paper, so the PRL SumTp binning is used (to be revisited if tables become public).
The LE beam configuration (playlists minerva1/7/9/13A-E, LE flux, GENIE 2.8.4-based tune without
the elastic-FSI fix) is prepared via `LE_CONFIG` but the LE open-data samples are not yet released.
"""
from xsec.binning import Binning
from xsec.util import load_measurement_module

_prl = load_measurement_module("ccqelike_3d_me_fhc")
SUMTP, QELike3D, fit_cells = _prl.SUMTP, _prl.QELike3D, _prl.fit_cells

PT = [0, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 2.5]
PZ = [1.5, 2.0, 2.5, 3.0, 3.5, 4.5]

LE_CONFIG = dict(
    playlists=["minerva1", "minerva7", "minerva9", "minerva13A", "minerva13B", "minerva13C", "minerva13D",
               "minerva13E"],
    flux_playlist="LE_FHC",
    tune_options=dict(fsi_elastic=False),  # LE GENIE 2.8.4 sample: elastic-FSI bug not fixed
    muon_energy_scale_err=0.02,
)


class QELike3DLEBinning(QELike3D):
    name = "ccqelike_3d_le_me_fhc"
    reco_binning = Binning([SUMTP, PT, PZ], ["sumtp", "pt", "pz"])
    true_binning = Binning([SUMTP, PT, PZ], ["sumtp", "pt", "pz"])


def background_scales(acc, meas):
    from xsec.sideband import fit_scales
    cell, keep = fit_cells(meas.reco_binning)  # integrates over the in-range p|| = [1.5, 4.5]
    scales, coef = fit_scales(acc, meas, cell, meas.fit_categories, keep=keep)
    return {c: v for c, v in scales.items() if c != "signal"}, coef
