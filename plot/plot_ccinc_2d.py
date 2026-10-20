#!/usr/bin/env python3
"""CC-inclusive d2sigma/dpT dp|| (ME FHC) vs PRD 104 092007.

Figures (figures/ccinc_2d_me_fhc/):
  xsec_pt_in_pz      d2sigma/dpT dp|| vs pT, one panel per p|| bin, with ratio to published
  xsec_pz_projection dsigma/dp|| (summed over the reported pT bins) with ratio
  frac_unc           stat / flux / total fractional uncertainties, ours vs published, per p|| bin
"""
import numpy as np

import common as C
from xsec import compare

M = "ccinc_2d_me_fhc"
REF_LABEL = "MINERvA, PRD 104 092007"
YSCALE = 1e-39


def main():
    ap = C.parser(M, "ccinc_v1_xsec.pkl", __doc__)
    ap.add_argument("--variant", choices=["mat", "strict"], default="mat",
                    help="mat: published MINERvA-101 signal convention (ccinc_v1); strict: ccinc_v1_strict")
    args = ap.parse_args()
    if args.variant == "strict":
        args.input = args.input.replace("ccinc_v1_xsec.pkl", "ccinc_v1_strict_xsec.pkl")
    sfx = "" if args.variant == "mat" else "_strict"
    run = C.run_module(M)
    PT, PZ = np.array(run.PT, float), np.array(run.PZ, float)
    res = C.load_pickle(args.input)
    ref, rcov = run.load_reference()
    x = res["xsec"]
    mask = ref != 0
    s = compare.summary(x, ref, rcov["total"], res["cov_total"], mask)
    stats = (f"χ² (published total cov.) = {s['chi2_refcov']:.1f} / {s['nbins']}    "
             f"integral ratio = {s['integral_ratio']:.3f}")
    print(M, args.variant, stats)
    npt, npz = len(PT) - 1, len(PZ) - 1
    es, et = C.err(res["cov_stat"]), C.err(res["cov_total"])
    er = C.err(rcov["total"])

    # --- per-p|| panels
    fig, pairs = C.ratio_figure(4, 4, (15, 15))
    for iz, (a, b) in enumerate(pairs):
        sl = slice(iz * npt, (iz + 1) * npt)
        rep = np.flatnonzero(mask[sl])
        lo, hi = rep[0], rep[-1] + 1
        idx = np.arange(sl.start + lo, sl.start + hi)
        e = PT[lo:hi + 1]
        C.draw_ref(a, e, ref[idx] / YSCALE, er[idx] / YSCALE)
        C.draw_ours(a, e, x[idx] / YSCALE, es[idx] / YSCALE, et[idx] / YSCALE, ms=3.5)
        C.draw_ratio(b, e, x[idx], et[idx], ref[idx], er[idx], ylim=(0.8, 1.2))
        a.set_ylim(bottom=0)
        a.set_title(f"{PZ[iz]:g} < p$_{{||}}$ < {PZ[iz + 1]:g} GeV/c")
        if iz % 4 == 0:
            b.set_ylabel("ours / pub.")
        if iz >= 12:
            b.set_xlabel(r"muon p$_T$ [GeV/c]")
    fig.supylabel(r"d$^2\sigma$/dp$_T$dp$_{||}$  [10$^{-39}$ cm$^2$ / (GeV/c)$^2$ / nucleon]", x=0.07)
    C.figure_legend(fig, REF_LABEL, y=0.955, extra=stats)
    C.save(fig, args.outdir, "xsec_pt_in_pz" + sfx, args.format)

    # --- dsigma/dp|| projection (sum over reported pT bins)
    wpt = np.diff(PT)
    P = C.projection([(np.arange(iz * npt, (iz + 1) * npt)[mask[iz * npt:(iz + 1) * npt]],
                       wpt[mask[iz * npt:(iz + 1) * npt]]) for iz in range(npz)], len(x))
    px, pr = P @ x, P @ ref
    pes, pet = C.err(P @ res["cov_stat"] @ P.T), C.err(P @ res["cov_total"] @ P.T)
    per = C.err(P @ rcov["total"] @ P.T)
    fig, [(a, b)] = C.ratio_figure(1, 1, (6.5, 5.5))
    C.draw_ref(a, PZ, pr / YSCALE, per / YSCALE, label=REF_LABEL)
    C.draw_ours(a, PZ, px / YSCALE, pes / YSCALE, pet / YSCALE, label=C.OURS_LABEL)
    C.draw_ratio(b, PZ, px, pet, pr, per, ylim=(0.85, 1.15))
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_ylabel(r"d$\sigma$/dp$_{||}$ [10$^{-39}$ cm$^2$/(GeV/c)/N]")
    b.set_xlabel(r"muon p$_{||}$ [GeV/c]")
    b.set_ylabel("ours / pub.")
    a.legend(C.legend_handles(), [C.OURS_LABEL, REF_LABEL], loc="lower left")
    a.set_title(r"ME FHC $\nu_\mu$ CC inclusive, summed over reported p$_T$ bins")
    C.save(fig, args.outdir, "xsec_pz_projection" + sfx, args.format)

    # --- fractional uncertainties
    fig, axs = C.plt.subplots(4, 4, figsize=(15, 13), sharey=True)
    fo = {k: C.err(c) / np.where(x != 0, np.abs(x), 1) for k, c in
          (("stat", res["cov_stat"]), ("flux", res["cov_bands"]["Flux"]), ("total", res["cov_total"]))}
    fr = {k: C.err(rcov[k]) / np.where(ref != 0, ref, 1) for k in ("stat", "flux", "total")}
    for iz, ax in enumerate(axs.flat):
        sl = slice(iz * npt, (iz + 1) * npt)
        rep = np.flatnonzero(mask[sl])
        lo, hi = rep[0], rep[-1] + 1
        idx = np.arange(sl.start + lo, sl.start + hi)
        C.frac_unc_panel(ax, PT[lo:hi + 1], [(k, fo[k][idx], fr[k][idx], ls) for k, ls in
                                              (("total", "-"), ("flux", "--"), ("stat", ":"))],
                         r"muon p$_T$ [GeV/c]" if iz >= 12 else "")
        if iz % 4:
            ax.set_ylabel("")
        ax.set_title(f"{PZ[iz]:g} < p$_{{||}}$ < {PZ[iz + 1]:g} GeV/c")
    axs.flat[0].set_ylim(0, 0.3)
    h, lab = axs.flat[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 0.94))
    fig.subplots_adjust(hspace=0.3, wspace=0.08)
    C.save(fig, args.outdir, "frac_unc" + sfx, args.format)


if __name__ == "__main__":
    main()
