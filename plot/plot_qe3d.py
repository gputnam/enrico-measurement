#!/usr/bin/env python3
"""QE-like d3sigma/dpT dp|| dSumTp (ME FHC) vs PRL 129 021803.

Figures (figures/ccqelike_3d_me_fhc/):
  xsec_pz<k>     d3sigma vs SumTp, one panel per pT bin, for p|| bin k, with ratio to published
  projections    dsigma/dp||, dsigma/dpT, dsigma/dSumTp (integrated over the other two axes)
  frac_unc       per-bin fractional uncertainties (stat, flux, total), ours vs published
"""
import numpy as np

import common as C
from xsec import compare

M = "ccqelike_3d_me_fhc"
REF_LABEL = "MINERvA, PRL 129 021803"
YSCALE = 1e-39
YLABEL = r"d$^3\sigma$/dp$_T$dp$_{||}d\Sigma T_p$  [10$^{-39}$ cm$^2$ / (GeV/c)$^2$ / GeV / nucleon]"
AXES = (("sumtp", r"$\Sigma T_p$ [GeV]", "GeV"), ("pt", r"muon p$_T$ [GeV/c]", "GeV/c"),
        ("pz", r"muon p$_{||}$ [GeV/c]", "GeV/c"))


def axis_projection(edges, vol, mask, axis):
    """P: density bins (sumtp fastest, then pt, pz) -> dsigma/d(axis), integrating the other axes over
    published (mask) bins."""
    shape = tuple(len(e) - 1 for e in edges)[::-1]  # (nz, ny, nx)
    ii = np.indices(shape).reshape(3, -1)[::-1]  # (ix, iy, iz) per flat bin
    w = np.diff(edges[axis])
    n = shape[::-1][axis]
    return C.projection([(np.flatnonzero((ii[axis] == k) & mask), vol[(ii[axis] == k) & mask] / w[k])
                         for k in range(n)], len(vol))


def main():
    args = C.parser(M, "qe3d_v1_xsec.pkl", __doc__).parse_args()
    run = C.run_module(M)
    SUMTP, PT, PZ = (np.array(e, float) for e in (run.SUMTP, run.PT, run.PZ))
    vol = run.QELike3D.reco_binning.widths()
    res = C.load_pickle(args.input)
    ref, rc = run.load_reference()
    x = res["xsec"]
    mask = ref != 0
    s = compare.summary(x, ref, rc["TotalCovariance"], res["cov_total"], mask)
    integ = (x * vol)[mask].sum() / (ref * vol)[mask].sum()
    stats = (f"integral ratio (Σσ·ΔV) = {integ:.3f}    χ² (sum of published + our cov.) = "
             f"{s['chi2_sumcov']:.0f} / {s['nbins']}")
    print(M, stats)
    es, et, er = C.err(res["cov_stat"]), C.err(res["cov_total"]), C.err(rc["TotalCovariance"])
    nx, ny, nz = len(SUMTP) - 1, len(PT) - 1, len(PZ) - 1

    # --- one figure per p|| bin, panels in pT, x = SumTp
    for k in range(nz):
        fig, pairs = C.ratio_figure(3, 4, (15, 11.5))
        for j, (a, b) in enumerate(pairs):
            if j >= ny:
                a.remove(), b.remove()
                continue
            idx = np.arange(nx) + nx * (j + ny * k)
            m = mask[idx]
            if not m.any():
                a.text(0.5, 0.5, "not reported", transform=a.transAxes, ha="center", color=C.MUTED)
            C.draw_ref(a, SUMTP, np.where(m, ref[idx], np.nan) / YSCALE, er[idx] / YSCALE)
            C.draw_ours(a, SUMTP, x[idx] / YSCALE, es[idx] / YSCALE, et[idx] / YSCALE, ms=3.5)
            C.draw_ratio(b, SUMTP, x[idx], et[idx], ref[idx], er[idx], ylim=(0.0, 2.0))
            a.set_title(f"{PT[j]:g} < p$_T$ < {PT[j + 1]:g} GeV/c")
            if j % 4 == 0:
                b.set_ylabel("ours / pub.")
            if j >= ny - 4:
                b.set_xlabel(AXES[0][1])
                b.tick_params(labelbottom=True)
        fig.supylabel(YLABEL, x=0.07)
        C.figure_legend(fig, REF_LABEL, y=0.965,
                        extra=f"{PZ[k]:g} < p$_{{||}}$ < {PZ[k + 1]:g} GeV/c        (all bins: {stats})")
        C.save(fig, args.outdir, f"xsec_pz{k}", args.format)

    # --- 1D projections
    edges = [SUMTP, PT, PZ]
    fig, pairs = C.ratio_figure(1, 3, (16, 5.5))
    for axis, (a, b) in enumerate(pairs):
        name, xl, unit = AXES[axis]
        P = axis_projection(edges, vol, mask, axis)
        px, pr = P @ x, P @ ref
        pes, pet = C.err(P @ res["cov_stat"] @ P.T), C.err(P @ res["cov_total"] @ P.T)
        per = C.err(P @ rc["TotalCovariance"] @ P.T)
        C.draw_ref(a, edges[axis], pr / YSCALE, per / YSCALE, label=REF_LABEL)
        C.draw_ours(a, edges[axis], px / YSCALE, pes / YSCALE, pet / YSCALE, label=C.OURS_LABEL)
        C.draw_ratio(b, edges[axis], px, pet, pr, per, ylim=(0.7, 1.3))
        a.set_yscale("log")
        if name == "pz":
            a.set_xscale("log")
        a.set_ylabel(rf"d$\sigma$/d{xl.split(' [')[0].replace('muon ', '')}  [10$^{{-39}}$ cm$^2$/{unit}/N]")
        b.set_xlabel(xl)
        b.set_ylabel("ours / pub.")
    pairs[0][0].legend(C.legend_handles(), [C.OURS_LABEL, REF_LABEL], loc="lower left")
    fig.suptitle(f"ME FHC QE-like, projections of the 3D result over published bins    ({stats})", y=0.97)
    C.save(fig, args.outdir, "projections", args.format)

    # --- per-bin fractional uncertainties, ours vs published
    fig, axs = C.plt.subplots(1, 3, figsize=(15, 5))
    for ax, (lab, ours, theirs) in zip(axs, (("statistical", res["cov_stat"], rc["StatisticalCovariance"]),
                                             ("flux", res["cov_bands"]["Flux"], rc["FluxCovariance"]),
                                             ("total", res["cov_total"], rc["TotalCovariance"]))):
        fo = C.err(ours)[mask] / np.abs(x[mask])
        fr = C.err(theirs)[mask] / ref[mask]
        good = np.isfinite(fo) & (fo > 0)
        ax.scatter(fr[good], fo[good], s=8, color=C.OURS, alpha=0.5, lw=0)
        lim = np.quantile(np.concatenate([fr[good], fo[good]]), 0.98) * 1.1
        ax.plot([0, lim], [0, lim], color=C.REF, lw=1.2)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"published fractional {lab} uncertainty")
        ax.set_ylabel(f"our fractional {lab} uncertainty")
        C.annotate(ax, f"median ours / published = {np.median(fo[good] / fr[good]):.2f}\n"
                       f"{good.sum()} bins; axes clipped at the 98th percentile", loc="upper left")
    fig.suptitle("ME FHC QE-like 3D: per-bin uncertainties, this reproduction vs PRL 129 021803")
    fig.tight_layout()
    C.save(fig, args.outdir, "frac_unc", args.format)


if __name__ == "__main__":
    main()
