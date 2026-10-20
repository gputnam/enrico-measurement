#!/usr/bin/env python3
"""CC1pi+ single-differential cross sections (ME FHC) vs arXiv:2605.24224.

Figures (figures/ccpip_me_fhc/):
  xsec_<var>  dsigma/d<var> with ratio to published, for each of the 7 variables
  summary     all variables on one page
"""
import os

import numpy as np

import common as C
from xsec import compare

M = "ccpip_me_fhc"
REF_LABEL = "MINERvA, arXiv:2605.24224"
LABELS = {  # axis label, unit of the differential
    "tpi": (r"pion kinetic energy T$_\pi$ [MeV]", "MeV", r"T$_\pi$"),
    "thetapi": (r"pion angle $\theta_\pi$ [deg]", "deg", r"$\theta_\pi$"),
    "pmu": (r"muon momentum p$_\mu$ [GeV/c]", "GeV/c", r"p$_\mu$"),
    "pzmu": (r"muon p$_{||}$ [GeV/c]", "GeV/c", r"p$_{||}$"),
    "ptmu": (r"muon p$_T$ [GeV/c]", "GeV/c", r"p$_T$"),
    "thetamu": (r"muon angle $\theta_\mu$ [deg]", "deg", r"$\theta_\mu$"),
    "q2": (r"Q$^2$ [GeV$^2$]", "GeV$^2$", r"Q$^2$"),
}


def load(var, indir, ref_x, ref_c, meta):
    res = C.load_pickle(os.path.join(indir, f"{var}_xsec.pkl"))
    x = res["xsec"] * 1e42
    t = ref_x[meta["table"]]
    edges = np.append(t[:, 0], t[-1, 1])
    r, er = t[:, 2], t[:, 3]
    mask = r > 0
    txt = f"integral ratio {x[mask].sum() / r[mask].sum():.3f}"
    if meta["cov"] and (meta["cov"], "stat") in ref_c:  # tables in units of 1e-84 = (1e-42)^2
        cov = ref_c[(meta["cov"], "stat")] + ref_c[(meta["cov"], "sys")]
        c2, n = compare.chi2(x, r, cov + res["cov_total"] * 1e84, mask)
        txt = f"χ² (published + our cov.) = {c2:.1f} / {n}\n" + txt
    return dict(edges=edges, x=x, es=C.err(res["cov_stat"]) * 1e42, et=C.err(res["cov_total"]) * 1e42, r=r, er=er,
                txt=txt)


def draw(a, b, d, var):
    C.draw_ref(a, d["edges"], d["r"], d["er"])
    C.draw_ours(a, d["edges"], d["x"], d["es"], d["et"])
    C.draw_ratio(b, d["edges"], d["x"], d["et"], d["r"], d["er"], ylim=(0.0, 2.0))
    xl, unit, sym = LABELS[var]
    a.set_ylim(bottom=0)
    a.set_ylabel(rf"d$\sigma$/d{sym} [10$^{{-42}}$ cm$^2$/{unit}/N]")
    b.set_xlabel(xl)
    b.set_ylabel("ours / pub.")
    a.set_title(d["txt"].replace("\n", "    "), fontsize=8.5, color=C.MUTED)


def main():
    ap = C.parser(M, "", __doc__)
    args = ap.parse_args()
    indir = args.input if os.path.isdir(args.input) else os.path.dirname(args.input)
    an, tables = C.run_module(M, "analysis"), C.run_module(M, "reference_tables")
    ref_x, ref_c = tables.load_xsec(), tables.load_cov()
    data = {v: load(v, indir, ref_x, ref_c, meta) for v, meta in an.VARIABLES.items()}
    for v, d in data.items():
        print(M, v, d["txt"].replace("\n", "   "))
        fig, [(a, b)] = C.ratio_figure(1, 1, (6.5, 5.5))
        draw(a, b, d, v)
        a.legend(C.legend_handles(), [C.OURS_LABEL, REF_LABEL], loc="lower left", bbox_to_anchor=(0, 1.08))
        C.save(fig, args.outdir, f"xsec_{v}", args.format)

    fig, pairs = C.ratio_figure(2, 4, (19, 10))
    for (a, b), (v, d) in zip(pairs, data.items()):
        draw(a, b, d, v)
    pairs[-1][0].remove(), pairs[-1][1].remove()
    C.figure_legend(fig, REF_LABEL, y=0.985, extra=r"ME FHC $\nu_\mu$ CC1$\pi^+$ on CH")
    C.save(fig, args.outdir, "summary", args.format)


if __name__ == "__main__":
    main()
