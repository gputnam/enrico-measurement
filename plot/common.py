"""Shared helpers for the comparison plots: inputs, styling, drawing and saving.

Convention: this reproduction = blue points (inner tick = stat, full bar = total uncertainty);
published MINERvA result = orange step line with a band of its total uncertainty.
"""
import argparse
import os
import pickle

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from xsec.util import load_measurement_module  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUTPUTS = os.path.join(ROOT, "outputs")
FIGURES = os.path.join(ROOT, "figures")

OURS = "#2a78d6"
REF = "#eb6834"
INK = "#2b2b2b"
MUTED = "#6b6b6b"
GRID = "#e4e4e0"
OURS_LABEL = "This reproduction (Open Data)"

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10, "legend.fontsize": 9,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "xtick.direction": "in", "ytick.direction": "in",
    "legend.frameon": False, "savefig.dpi": 150, "figure.dpi": 100,
})


def parser(measurement, default_input, doc=None):
    ap = argparse.ArgumentParser(description=doc)
    ap.add_argument("--input", default=os.path.join(OUTPUTS, measurement, "v1", default_input),
                    help="extracted cross-section pickle (run.py *_xsec.pkl output)")
    ap.add_argument("--outdir", default=os.path.join(FIGURES, measurement))
    ap.add_argument("--format", nargs="+", default=["pdf", "png"])
    return ap


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def run_module(measurement, module="run"):
    return load_measurement_module(measurement, module)


def err(cov):
    return np.sqrt(np.clip(np.diag(cov), 0, None))


def projection(blocks, n):
    """Matrix P (len(blocks) x n) summing the listed index groups with per-index weights.
    blocks: list of (indices, weights)."""
    P = np.zeros((len(blocks), n))
    for r, (idx, w) in enumerate(blocks):
        P[r, idx] = w
    return P


def draw_ref(ax, edges, y, e, label=None):
    """Published result: step line with a filled +-1 sigma total-uncertainty band."""
    edges = np.asarray(edges, float)
    yy, ee = np.append(y, y[-1]), np.append(e, e[-1])
    ax.fill_between(edges, yy - ee, yy + ee, step="post", color=REF, alpha=0.25, lw=0)
    ax.step(edges, yy, where="post", color=REF, lw=1.5, label=label)


def draw_ours(ax, edges, y, e_stat, e_tot, label=None, ms=4):
    """This reproduction: markers at bin centers, total (thin) and statistical (capped) errors."""
    edges = np.asarray(edges, float)
    c = 0.5 * (edges[1:] + edges[:-1])
    ax.errorbar(c, y, yerr=e_tot, fmt="none", ecolor=OURS, elinewidth=1.0)
    ax.errorbar(c, y, yerr=e_stat, fmt="o", color=OURS, ms=ms, mec="white", mew=0.6, ecolor=OURS,
                elinewidth=1.6, capsize=2, label=label)


def draw_ratio(ax, edges, ours, ours_tot, ref, ref_tot, ylim=(0.5, 1.5)):
    """ours / published with published fractional uncertainty band around 1."""
    edges = np.asarray(edges, float)
    good = ref > 0
    rr = np.where(good, ref_tot / np.where(good, ref, 1), np.nan)
    ax.fill_between(edges, 1 - np.append(rr, rr[-1]), 1 + np.append(rr, rr[-1]), step="post", color=REF,
                    alpha=0.25, lw=0)
    ax.axhline(1, color=REF, lw=1.2)
    c = 0.5 * (edges[1:] + edges[:-1])
    r = np.where(good, ours / np.where(good, ref, 1), np.nan)
    re = np.where(good, ours_tot / np.where(good, ref, 1), np.nan)
    ax.errorbar(c, r, yerr=re, fmt="o", color=OURS, ms=3.5, mec="white", mew=0.5, elinewidth=1.0)
    # markers outside the range: arrows at the edge
    lo, hi = ylim
    for xc, v in zip(c, r):
        if np.isfinite(v) and (v > hi or v < lo):
            y = hi - 0.04 * (hi - lo) if v > hi else lo + 0.04 * (hi - lo)
            ax.annotate("", (xc, y), (xc, y - np.sign(v - 1) * 0.12 * (hi - lo)),
                        arrowprops=dict(arrowstyle="->", color=OURS, lw=1.0))
    ax.set_ylim(*ylim)
    q = 0.25 * (hi - lo)
    ax.set_yticks([lo + q, 1, hi - q])


def legend_handles():
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    return [Line2D([], [], marker="o", color=OURS, ms=5, mec="white", lw=1.2),
            (Patch(color=REF, alpha=0.25, lw=0), Line2D([], [], color=REF, lw=1.5))]


def figure_legend(fig, ref_label, y=0.995, extra=None):
    h = legend_handles()
    labels = [OURS_LABEL + " (capped bar: stat., full bar: total)", ref_label + " (band: total unc.)"]
    fig.legend(h, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, y))
    if extra:
        fig.text(0.5, y - 0.025, extra, ha="center", va="top", color=MUTED, fontsize=9)


def annotate(ax, text, loc="upper right"):
    x, ha = (0.97, "right") if "right" in loc else (0.03, "left")
    y, va = (0.95, "top") if "upper" in loc else (0.05, "bottom")
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, color=MUTED, fontsize=8)


def save(fig, outdir, name, formats):
    os.makedirs(outdir, exist_ok=True)
    for fmt in formats:
        fig.savefig(os.path.join(outdir, f"{name}.{fmt}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {os.path.join(outdir, name)}.{{{','.join(formats)}}}")


def ratio_figure(nrow_panels, ncol, figsize, height_ratios=(3, 1.2)):
    """Grid of (panel, ratio) pairs: returns fig, list of (ax_main, ax_ratio) row-major."""
    fig = plt.figure(figsize=figsize)
    outer = fig.add_gridspec(nrow_panels, ncol, hspace=0.42, wspace=0.28)
    pairs = []
    for r in range(nrow_panels):
        for c in range(ncol):
            inner = outer[r, c].subgridspec(2, 1, height_ratios=height_ratios, hspace=0.05)
            a = fig.add_subplot(inner[0])
            b = fig.add_subplot(inner[1], sharex=a)
            a.tick_params(labelbottom=False)
            pairs.append((a, b))
    return fig, pairs


def frac_unc_panel(ax, edges, items, xlabel, title=None, logx=False):
    """items: list of (label, ours_frac, ref_frac, linestyle). Ours solid blue, published dashed orange."""
    for lab, fo, fr, ls in items:
        draw_frac(ax, edges, fo, OURS, ls, f"{lab} (ours)")
        draw_frac(ax, edges, fr, REF, ls, f"{lab} (published)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("fractional uncertainty")
    ax.set_ylim(bottom=0)
    if logx:
        ax.set_xscale("log")
    if title:
        ax.set_title(title)


def draw_frac(ax, edges, f, color, ls, label):
    edges = np.asarray(edges, float)
    ax.step(edges, np.append(f, f[-1]), where="post", color=color, ls=ls, lw=1.4, label=label)
