"""The measurement corpus: measured cross sections compared with a generator model.

A corpus is a directory (the data model of the closed loop in orchestrate/).  It is read by the loop
agents and therefore contains ONLY data and the model under study -- no prediction of any other
generator (in particular not the MINERvA GENIE/MnvTune prediction):

    <corpus>/index.json            entries (id, measurement dir, class, kwargs, title, origin)
    <corpus>/summary.txt           one line per entry: model chi2 / pulls / normalisation
    <corpus>/<id>/definition.md    signal definition, observables, binning, motivation
    <corpus>/<id>/result.npz       data, covariances, binning (machine readable)
    <corpus>/<id>/comparison.txt   per bin: edges, data, errors, model, ratio, pull
    <corpus>/<id>/projections.txt  1D projections on every axis (data +- total, model)
    <corpus>/<id>/metrics.json     chi2 (full covariance and diagonal), pulls, normalisation
    <corpus>/<id>/model.npy        current model prediction

The MINERvA GENIE prediction (`truth_xsec` of a result) is kept only in a *private* directory
(`--private`, human-only, e.g. inside the plots area) and is used only for the human plots
(<plots>/<id>/*.png), which also must live outside the corpus.

Each measurement directory measurements/<m>/ declares its corpus entries in corpus.json:

    [{"id": "ccpip_me_fhc_tpi", "class": "CCPip", "kwargs": {"variable": "tpi"},
      "title": "CC1pi+ dsigma/dTpi", "result": "ccpip_me_fhc/v1/tpi_xsec.pkl"}]

("result" is relative to the results root given on the command line.)

    python -m xsec.corpus add     --corpus C --measurement <m> --results-root R [--origin seed] [--private P]
    python -m xsec.corpus predict --corpus C --model model.root --model-version v3 [--plots D --private P] [--only id ...]
    python -m xsec.corpus summary --corpus C
"""
import argparse
import datetime
import json
import os
import pickle

import numpy as np

from .compare import chi2
from .util import MEASUREMENTS, load_measurement_module

MAX_FRAC_UNC = 1.0  # bins with total fractional uncertainty >= this are not used in chi2 / pulls
GENIE_LABEL = "MINERvA GENIE (MnvTune {tune})"


# ---------------------------------------------------------------- helpers
def measurement_entries(measurement):
    with open(os.path.join(MEASUREMENTS, measurement, "corpus.json")) as f:
        return json.load(f)


def make_measurement(entry):
    mod = load_measurement_module(entry["measurement"])
    return getattr(mod, entry["class"])(**entry.get("kwargs", {}))


def load_index(corpus):
    p = os.path.join(corpus, "index.json")
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return json.load(f)


def save_index(corpus, index):
    with open(os.path.join(corpus, "index.json"), "w") as f:
        json.dump(index, f, indent=2)


def err(cov):
    return np.sqrt(np.clip(np.diag(cov), 0, None))


def valid_bins(x, cov):
    s = err(cov)
    return (s > 0) & np.isfinite(x) & (s < MAX_FRAC_UNC * np.abs(x))


def axis_projection(binning, axis, mask=None):
    """Matrix P (n_axis x n_inner): density -> dsigma/d(axis), integrating the other axes."""
    nb = binning.nb
    idx = np.indices(nb[::-1]).reshape(len(nb), -1)[::-1]  # per inner bin (x fastest): index per axis
    vol = binning.widths()
    w = np.diff(binning.edges[axis])
    m = np.ones(vol.size, dtype=bool) if mask is None else mask
    P = np.zeros((nb[axis], vol.size))
    for k in range(nb[axis]):
        sel = (idx[axis] == k) & m
        P[k, sel] = vol[sel] / w[k]
    return P


def bin_edges_table(binning):
    nb = binning.nb
    idx = np.indices(nb[::-1]).reshape(len(nb), -1)[::-1]
    lo = [binning.edges[d][idx[d]] for d in range(len(nb))]
    hi = [binning.edges[d][idx[d] + 1] for d in range(len(nb))]
    return lo, hi


# ---------------------------------------------------------------- private (human-only) GENIE store
def genie_path(private, eid):
    return os.path.join(private, f"{eid}_genie.npz")


def store_genie(private, eid, genie, label):
    os.makedirs(private, exist_ok=True)
    np.savez(genie_path(private, eid), genie=np.asarray(genie), label=label)


def load_genie(private, eid):
    p = genie_path(private, eid) if private else None
    if not p or not os.path.exists(p):
        return None, None
    d = np.load(p)
    return d["genie"], str(d["label"])


# ---------------------------------------------------------------- add
def add(corpus, measurement, results_root, origin="seed", ids=None, private=None):
    os.makedirs(corpus, exist_ok=True)
    index = load_index(corpus)
    mdir = os.path.join(MEASUREMENTS, measurement)
    added = []
    for e in measurement_entries(measurement):
        if ids and e["id"] not in ids:
            continue
        entry = dict(e, measurement=measurement, origin=origin,
                     added=datetime.datetime.now().isoformat(timespec="seconds"))
        with open(os.path.join(results_root, e["result"]), "rb") as f:
            res = pickle.load(f)
        meas = make_measurement(entry)
        TB = meas.true_binning
        n = int(TB.inner_mask().sum())
        assert len(res["xsec"]) == n, f"{e['id']}: result has {len(res['xsec'])} bins, binning {n}"
        if private and "truth_xsec" in res:
            store_genie(private, e["id"], res["truth_xsec"], GENIE_LABEL.format(tune=meas.tune_name))
        d = os.path.join(corpus, e["id"])
        os.makedirs(d, exist_ok=True)
        np.savez_compressed(os.path.join(d, "result.npz"), xsec=res["xsec"], cov_total=res["cov_total"],
                            cov_stat=res["cov_stat"], cov_sys=res["cov_sys"], names=np.array(TB.names),
                            **{f"edges_{i}": np.asarray(ed) for i, ed in enumerate(TB.edges)},
                            flux_integral=res.get("flux_integral", np.nan), n_targets=res.get("n_targets", np.nan),
                            pot=res.get("pot", np.nan))
        src = next((os.path.join(mdir, f) for f in ("definition.md", "README.md") if os.path.exists(os.path.join(mdir, f))), None)
        with open(os.path.join(d, "definition.md"), "w") as f:
            f.write(f"# {e['id']}: {e.get('title', '')}\n\n")
            f.write(f"- measurement code: measurement/measurements/{measurement}/ (class `{e['class']}`, kwargs {e.get('kwargs', {})})\n")
            f.write(f"- origin: {origin}\n- observables (binning axes): {', '.join(TB.names)}; {n} bins\n\n")
            if e.get("motivation"):
                f.write(f"## Motivation\n\n{e['motivation']}\n\n")
            if src:
                f.write(f"## Measurement description ({os.path.basename(src)})\n\n")
                f.write(open(src).read())
        index = [i for i in index if i["id"] != e["id"]] + [entry]
        added.append(e["id"])
    save_index(corpus, index)
    print(f"added {added} to {corpus}")
    return added


# ---------------------------------------------------------------- predict / compare
def compare_entry(corpus, entry, model=None, model_version="none", plots=None, private=None):
    d = os.path.join(corpus, entry["id"])
    r = dict(np.load(os.path.join(d, "result.npz")))
    meas = make_measurement(entry)
    TB = meas.true_binning
    x, C, Cs = r["xsec"], r["cov_total"], r["cov_stat"]
    if model is not None:
        p = model.predict(meas)
        np.save(os.path.join(d, "model.npy"), p)
    elif os.path.exists(os.path.join(d, "model.npy")):
        p = np.load(os.path.join(d, "model.npy"))
    else:
        p = np.full_like(x, np.nan)
    with open(os.path.join(d, "model_version.txt"), "w") as f:
        f.write(model_version + "\n")
    s, ss = err(C), err(Cs)
    used = valid_bins(x, C)
    vol = TB.widths()
    met = dict(id=entry["id"], model_version=model_version, nbins=int(len(x)), nbins_used=int(used.sum()))
    if np.all(np.isfinite(p[used])):
        c2, ndf = chi2(p, x, C, used)
        z = (p - x)[used] / s[used]
        tot_d, tot_p = (x * vol)[used].sum(), (p * vol)[used].sum()
        tot_err = float(np.sqrt(vol[used] @ C[np.ix_(used, used)] @ vol[used]))
        met["model"] = dict(chi2=c2, ndf=ndf, chi2_diag=float((z ** 2).sum()), median_abs_pull=float(np.median(np.abs(z))),
                            frac_within_1sigma=float(np.mean(np.abs(z) < 1)), frac_within_2sigma=float(np.mean(np.abs(z) < 2)),
                            norm_ratio=float(tot_p / tot_d), norm_ratio_err=float(tot_p / tot_d * tot_err / tot_d))
    with open(os.path.join(d, "metrics.json"), "w") as f:
        json.dump(met, f, indent=2)
    lo, hi = bin_edges_table(TB)
    with np.errstate(divide="ignore", invalid="ignore"):
        rm, zm = p / x, (p - x) / s
    with open(os.path.join(d, "comparison.txt"), "w") as f:
        f.write(f"# {entry['id']}: {entry.get('title', '')}\n# model version: {model_version}\n")
        f.write("# units: cross section in cm^2/nucleon per unit bin volume of the axes; pull = (model - data)/sigma_tot\n")
        f.write(f"# used=1: bin enters chi2/pulls (total fractional uncertainty < {MAX_FRAC_UNC})\n")
        if "model" in met:
            m = met["model"]
            f.write(f"# model: chi2 {m['chi2']:.1f}/{m['ndf']} (full cov), chi2_diag {m['chi2_diag']:.1f}, "
                    f"median|pull| {m['median_abs_pull']:.2f}, norm ratio {m['norm_ratio']:.3f} +- {m['norm_ratio_err']:.3f}\n")
        cols = ["bin"] + [f"{n}_lo {n}_hi" for n in TB.names] + ["data", "stat_err", "tot_err", "model", "model/data",
                                                                 "pull", "used"]
        f.write("# " + " ".join(cols) + "\n")
        for i in range(len(x)):
            ed = " ".join(f"{lo[k][i]:.5g} {hi[k][i]:.5g}" for k in range(TB.ndim))
            f.write(f"{i} {ed} {x[i]:.4e} {ss[i]:.3e} {s[i]:.3e} {p[i]:.4e} {rm[i]:.3f} {zm[i]:.2f} {int(used[i])}\n")
    with open(os.path.join(d, "projections.txt"), "w") as f:
        f.write(f"# {entry['id']} 1D projections (integrated over the other axes, all bins); model version {model_version}\n")
        for ax in range(TB.ndim):
            P = axis_projection(TB, ax)
            px, pe, pp = P @ x, np.sqrt(np.clip(np.diag(P @ C @ P.T), 0, None)), P @ p
            e = TB.edges[ax]
            f.write(f"\n## dsigma/d{TB.names[ax]}\n# lo hi data tot_err model model/data\n")
            for k in range(len(px)):
                f.write(f"{e[k]:.5g} {e[k + 1]:.5g} {px[k]:.4e} {pe[k]:.3e} {pp[k]:.4e} "
                        f"{pp[k] / px[k] if px[k] else np.nan:.3f}\n")
    if plots:
        from .corpus_plots import plot_entry
        g, glabel = load_genie(private, entry["id"])
        plot_entry(entry, TB, x, ss, s, g, glabel, p, C, met, model_version, os.path.join(plots, entry["id"]))
    return met


def summary(corpus):
    index = load_index(corpus)
    lines = ["# Corpus summary: model vs data, per measurement",
             "# chi2: full covariance over used bins; med|z|: median |pull|; norm: integrated model/data (+- data unc.)",
             f"# {'id':<28} {'nbin':>5} {'ver':>5} | {'chi2/ndf':>16} {'med|z|':>6} {'norm':>13} | origin"]
    for e in index:
        p = os.path.join(corpus, e["id"], "metrics.json")
        if not os.path.exists(p):
            lines.append(f"  {e['id']:<28} (no comparison yet)")
            continue
        m = json.load(open(p))
        if "model" in m:
            q = m["model"]
            part = (f"{q['chi2']:>9.1f}/{q['ndf']:<6d} {q['median_abs_pull']:>6.2f} "
                    f"{q['norm_ratio']:>6.3f}+-{q['norm_ratio_err']:<5.3f}")
        else:
            part = f"{'n/a':>16} {'':>6} {'':>13}"
        lines.append(f"  {e['id']:<28} {m['nbins_used']:>5} {m['model_version']:>5} | {part} | {e.get('origin', '')}")
    txt = "\n".join(lines) + "\n"
    with open(os.path.join(corpus, "summary.txt"), "w") as f:
        f.write(txt)
    return txt


def predict(corpus, model_path, model_version, plots=None, only=None, private=None):
    from .genpred import GeneratorSample
    model = GeneratorSample(model_path) if model_path else None
    for e in load_index(corpus):
        if only and e["id"] not in only:
            continue
        met = compare_entry(corpus, e, model, model_version, plots, private)
        mm = met.get("model", {})
        print(f"  {e['id']:<28} model chi2 {mm.get('chi2', np.nan):.1f}/{mm.get('ndf', 0)} norm {mm.get('norm_ratio', np.nan):.3f}")
    print(summary(corpus))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    priv = dict(default=os.environ.get("MVA_PRIVATE"),
                help="human-only directory for the GENIE predictions (default $MVA_PRIVATE; none: GENIE dropped)")
    a = sub.add_parser("add")
    a.add_argument("--corpus", required=True)
    a.add_argument("--measurement", required=True, help="directory name under measurements/")
    a.add_argument("--results-root", required=True)
    a.add_argument("--origin", default="seed")
    a.add_argument("--ids", nargs="*")
    a.add_argument("--private", **priv)
    p = sub.add_parser("predict")
    p.add_argument("--corpus", required=True)
    p.add_argument("--model", help="nuflow ROOT file (omit: reuse stored model predictions)")
    p.add_argument("--model-version", default="none")
    p.add_argument("--plots", help="directory for human-only plots")
    p.add_argument("--private", **priv)
    p.add_argument("--only", nargs="*")
    s = sub.add_parser("summary")
    s.add_argument("--corpus", required=True)
    args = ap.parse_args()
    if args.cmd == "add":
        add(args.corpus, args.measurement, args.results_root, args.origin, args.ids, args.private)
    elif args.cmd == "predict":
        predict(args.corpus, args.model, args.model_version, args.plots, args.only, args.private)
    else:
        print(summary(args.corpus))


if __name__ == "__main__":
    main()
