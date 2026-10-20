#!/usr/bin/env python3
"""QE-like 3D in the LE/ME-paper binning (arXiv:2606.00745), ME FHC.

That paper has no numerical release, so the comparison is to PRL 129 021803 (same ME data) with both
results projected onto their common coarse bins (pT as in the LE/ME paper, p|| 1.5-3.5 and 3.5-4.5 GeV/c;
the projection of measurements/ccqelike_3d_le_me_fhc/run.py).

Figures (figures/ccqelike_3d_le_me_fhc/):
  xsec_coarse_pz<k>  d3sigma vs SumTp, one panel per coarse pT bin, for coarse p|| bin k, with ratio
"""
import numpy as np

import common as C
from xsec import compare

M = "ccqelike_3d_le_me_fhc"
REF_LABEL = "MINERvA PRL 129 021803, projected to common bins"
YSCALE = 1e-39
YLABEL = r"d$^3\sigma$/dp$_T$dp$_{||}d\Sigma T_p$  [10$^{-39}$ cm$^2$ / (GeV/c)$^2$ / GeV / nucleon]"


def main():
    args = C.parser(M, "qele_v1_xsec.pkl", __doc__).parse_args()
    run = C.run_module(M)
    SUMTP, CPT, CPZ = (np.array(e, float) for e in (run.SUMTP, run.COARSE_PT, run.COARSE_PZ))
    res = C.load_pickle(args.input)
    Po, Pr = run.projector(run.SUMTP, run.PT, run.PZ), run.projector(run.SUMTP, run.PRL_PT, run.PRL_PZ)
    ref, rc = run.prl_run.load_reference()
    # integrated cross sections per coarse bin -> densities
    nx, ny, nz = len(SUMTP) - 1, len(CPT) - 1, len(CPZ) - 1
    vol = np.array([(SUMTP[i + 1] - SUMTP[i]) * (CPT[j + 1] - CPT[j]) * (CPZ[k + 1] - CPZ[k])
                    for k in range(nz) for j in range(ny) for i in range(nx)])
    Do, Dr = Po / vol[:, None], Pr / vol[:, None]
    x, r = Do @ res["xsec"], Dr @ ref
    cs, ct = Do @ res["cov_stat"] @ Do.T, Do @ res["cov_total"] @ Do.T
    cr = Dr @ rc["TotalCovariance"] @ Dr.T
    mask = r > 0
    s = compare.summary(x, r, cr, ct, mask)
    integ = (x * vol)[mask].sum() / (r * vol)[mask].sum()
    stats = f"integral ratio = {integ:.3f}    χ² (sum of covariances) = {s['chi2_sumcov']:.0f} / {s['nbins']}"
    print(M, stats)
    es, et, er = C.err(cs), C.err(ct), C.err(cr)

    for k in range(nz):
        fig, pairs = C.ratio_figure(3, 4, (15, 11.5))
        for j, (a, b) in enumerate(pairs):
            if j >= ny:
                a.remove(), b.remove()
                continue
            idx = np.arange(nx) + nx * (j + ny * k)
            C.draw_ref(a, SUMTP, np.where(mask[idx], r[idx], np.nan) / YSCALE, er[idx] / YSCALE)
            C.draw_ours(a, SUMTP, x[idx] / YSCALE, es[idx] / YSCALE, et[idx] / YSCALE, ms=3.5)
            C.draw_ratio(b, SUMTP, x[idx], et[idx], r[idx], er[idx], ylim=(0.0, 2.0))
            a.set_title(f"{CPT[j]:g} < p$_T$ < {CPT[j + 1]:g} GeV/c")
            if j % 4 == 0:
                b.set_ylabel("ours / pub.")
            if j >= ny - 4:
                b.set_xlabel(r"$\Sigma T_p$ [GeV]")
        fig.supylabel(YLABEL, x=0.07)
        C.figure_legend(fig, REF_LABEL, y=0.965,
                        extra=f"LE/ME-paper binning (arXiv:2606.00745), ME FHC, {CPZ[k]:g} < p$_{{||}}$ < "
                              f"{CPZ[k + 1]:g} GeV/c        (all bins: {stats})")
        C.save(fig, args.outdir, f"xsec_coarse_pz{k}", args.format)


if __name__ == "__main__":
    main()
