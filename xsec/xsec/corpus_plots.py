"""Human-facing comparison plots for corpus entries: data vs MINERvA GENIE vs the generator model.

Not read by the loop agents.  For every entry:
  projections.png   dsigma/d(axis) for every binning axis (other axes integrated), with ratio to data
  slices_<k>.png    (ndim >= 2) panels of axis 0 in bins of axis 1, one figure per bin of the
                    remaining axes
"""
import itertools
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

DATA, GENIE, MODEL, MUTED, GRID = "#2b2b2b", "#eb6834", "#2a78d6", "#6b6b6b", "#e4e4e0"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "xtick.direction": "in", "ytick.direction": "in",
                     "legend.frameon": False, "savefig.dpi": 120})


def _draw(ax, rax, edges, y, ey, ys, ey_stat=None):
    """ys: list of (label, color, values)."""
    edges = np.asarray(edges, float)
    c = 0.5 * (edges[1:] + edges[:-1])
    ax.errorbar(c, y, yerr=ey, fmt="o", color=DATA, ms=3.5, elinewidth=1.0, label="data (total unc.)", zorder=3)
    if ey_stat is not None:
        ax.errorbar(c, y, yerr=ey_stat, fmt="none", ecolor=DATA, elinewidth=1.8, capsize=2, zorder=3)
    good = y != 0
    rr = np.where(good, ey / np.where(good, np.abs(y), 1), np.nan)
    rax.fill_between(edges, 1 - np.append(rr, rr[-1]), 1 + np.append(rr, rr[-1]), step="post", color=DATA,
                     alpha=0.15, lw=0)
    rax.axhline(1, color=DATA, lw=0.8)
    for lab, col, v in ys:
        if not np.all(np.isfinite(v)):
            continue
        ax.step(edges, np.append(v, v[-1]), where="post", color=col, lw=1.4, label=lab)
        r = np.where(good, v / np.where(good, y, 1), np.nan)
        rax.step(edges, np.append(r, r[-1]), where="post", color=col, lw=1.2)
    rax.set_ylim(0.0, 2.0)
    rax.set_ylabel("pred / data", fontsize=8)


def _pair(fig, spec):
    inner = spec.subgridspec(2, 1, height_ratios=(3, 1.2), hspace=0.05)
    a = fig.add_subplot(inner[0])
    b = fig.add_subplot(inner[1], sharex=a)
    a.tick_params(labelbottom=False)
    return a, b


def _stats_line(met, model_version):
    out = []
    for k, lab in (("model", f"model {model_version}"), ("genie", "GENIE")):
        if k in met:
            m = met[k]
            out.append(f"{lab}: χ² {m['chi2']:.0f}/{m['ndf']}, norm {m['norm_ratio']:.2f}")
    return "    ".join(out)


def plot_entry(entry, TB, x, e_stat, e_tot, genie, model, cov, met, model_version, outdir):
    from .corpus import axis_projection
    os.makedirs(outdir, exist_ok=True)
    title = f"{entry['id']}: {entry.get('title', '')}\n{_stats_line(met, model_version)}"
    preds = lambda P, g, m: [(met.get("genie_label", "GENIE"), GENIE, P @ g), (f"model {model_version}", MODEL, P @ m)]
    # projections
    nd = TB.ndim
    fig = plt.figure(figsize=(4.6 * nd, 4.4))
    gs = fig.add_gridspec(1, nd, wspace=0.3)
    for axn in range(nd):
        P = axis_projection(TB, axn)
        a, b = _pair(fig, gs[0, axn])
        _draw(a, b, TB.edges[axn], P @ x, np.sqrt(np.clip(np.diag(P @ cov @ P.T), 0, None)), preds(P, genie, model))
        a.set_ylabel(f"dσ/d{TB.names[axn]} [cm²/nucleon/unit]")
        b.set_xlabel(TB.names[axn])
        if axn == 0:
            a.legend(fontsize=7)
    fig.suptitle(title, fontsize=9, y=1.02)
    fig.savefig(os.path.join(outdir, "projections.png"), bbox_inches="tight")
    plt.close(fig)
    if nd < 2:
        return
    # slices: x = axis 0, panels = axis 1, one figure per combination of axes >= 2
    nb = TB.nb
    idx = np.indices(nb[::-1]).reshape(nd, -1)[::-1]
    rest = list(itertools.product(*[range(n) for n in nb[2:]]))
    for fi, combo in enumerate(rest):
        n1 = nb[1]
        ncol = min(4, n1)
        nrow = int(np.ceil(n1 / ncol))
        fig = plt.figure(figsize=(4.0 * ncol, 3.6 * nrow))
        gs = fig.add_gridspec(nrow, ncol, hspace=0.45, wspace=0.3)
        sub = " ".join(f"{TB.names[2 + k]}∈[{TB.edges[2 + k][j]:.3g},{TB.edges[2 + k][j + 1]:.3g}]"
                       for k, j in enumerate(combo))
        for j in range(n1):
            sel = idx[1] == j
            for k, jj in enumerate(combo):
                sel &= idx[2 + k] == jj
            ii = np.flatnonzero(sel)
            ii = ii[np.argsort(idx[0][ii])]
            a, b = _pair(fig, gs[j // ncol, j % ncol])
            _draw(a, b, TB.edges[0], x[ii], e_tot[ii], [(met.get("genie_label", "GENIE"), GENIE, genie[ii]),
                                                         (f"model {model_version}", MODEL, model[ii])], e_stat[ii])
            a.set_title(f"{TB.names[1]}∈[{TB.edges[1][j]:.3g},{TB.edges[1][j + 1]:.3g}]", fontsize=8)
            b.set_xlabel(TB.names[0])
            if j == 0:
                a.legend(fontsize=6)
        fig.suptitle(title + ("\n" + sub if sub else ""), fontsize=9, y=1.0)
        fig.savefig(os.path.join(outdir, f"slices_{fi}.png"), bbox_inches="tight")
        plt.close(fig)
