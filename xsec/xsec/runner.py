"""Event loop over DF files: fills data, MC (per universe) and truth accumulators.

A measurement object provides
    tune_name, minos_p_offset, reco_binning, true_binning, samples, signal_sample,
    bkg_categories, reco_columns, truth_columns,
    universes(tune) -> list[Universe]
    preselect(ev) -> mask                (cheap, CV-independent superset; reco rows)
    reco(ev, shift) -> (vars, {sample: mask})
    truth(ev) -> vars                    (true kinematics; reco and truth rows)
    signal(ev) -> mask                   (signal definition incl. phase space)
    category(ev) -> int array            (background category index for non-signal)
MC histograms are scaled per playlist to the data POT of that playlist.
"""
import multiprocessing as mp
import time

import numpy as np

from . import io
from .histo import Accumulator, fill, fill2d
from .kinematics import CV, Event
from .tunes import Tune

_STATE = {}


def _ratio_matrix(ev, tune, universes, truth):
    """(n, U) ratio for every universe evaluated on CV kinematics (lateral: extra ratio only)."""
    n = len(ev)
    R = np.ones((n, len(universes)))
    for j, u in enumerate(universes):
        if u.ratio is None or (truth and u.reco_only):
            continue
        R[:, j] = u.ratio(ev, tune)
    return R


def process_mc(path, playlist, scale):
    meas, tune, universes = _STATE["meas"], _STATE["tune"], _STATE["univ"]
    RB, TB = meas.reco_binning, meas.true_binning
    U = len(universes)
    acc = Accumulator()
    # ---------------- reco rows
    df = io.load(path, "reco", meas.reco_columns + meas.mc_reco_columns + meas.truth_columns)
    ev = Event(df, True, playlist, meas.minos_p_offset)
    if len(ev):
        pre = meas.preselect(ev)
        ev = Event(df[pre].reset_index(drop=True), True, playlist, meas.minos_p_offset)
    if len(ev):
        sig = meas.signal(ev)
        cat = meas.category(ev)
        t_idx = TB.index(*meas.truth(ev))
        lat = [j for j, u in enumerate(universes) if u.kind == "lat"]
        notlat = [j for j, u in enumerate(universes) if u.kind != "lat"]
        R = _ratio_matrix(ev, tune, universes, truth=False)
        w_cv = tune.cv_weight(ev, CV) * scale
        groups = [(CV, notlat)] + [(universes[j].shift, [j]) for j in lat]
        for shift, cols in groups:
            rv, samples = meas.reco(ev, shift)
            r_idx = RB.index(*rv)
            w = (w_cv if shift is CV else tune.cv_weight(ev, shift) * scale)[:, None] * R[:, cols]
            for s, m in samples.items():
                ms = m & sig
                h = fill(RB.size, r_idx[ms], w[ms])
                _scatter(acc, f"sig/{s}", h, cols, U)
                for ic, cname in enumerate(meas.bkg_categories):
                    mb = m & ~sig & (cat == ic)
                    _scatter(acc, f"bkg/{s}/{cname}", fill(RB.size, r_idx[mb], w[mb]), cols, U)
                if s == meas.signal_sample:
                    _scatter(acc, "effnum", fill(TB.size, t_idx[ms], w[ms]), cols, U)
                    if shift is CV:
                        acc.add("mig", fill2d(RB.size, TB.size, r_idx[ms], t_idx[ms], w_cv[ms]))
    # ---------------- truth rows
    dft = io.load(path, "truth", meas.truth_columns)
    evt = Event(dft, True, playlist, is_truth=True)
    if len(evt):
        sig = meas.signal(evt)
        evt = Event(dft[sig].reset_index(drop=True), True, playlist, is_truth=True)
    if len(evt):
        t_idx = TB.index(*meas.truth(evt))
        R = _ratio_matrix(evt, tune, universes, truth=True)
        w = tune.cv_weight(evt)[:, None] * R * scale
        acc.add("effden", fill(TB.size, t_idx, w))
    return acc


def _scatter(acc, key, h, cols, U):
    full = np.zeros((h.shape[0], U))
    full[:, cols] = h
    acc.add(key, full)


def process_data(path, playlist):
    meas = _STATE["meas"]
    RB = meas.reco_binning
    acc = Accumulator()
    df = io.load(path, "reco", meas.reco_columns)
    ev = Event(df, False, playlist, meas.minos_p_offset)
    if len(ev):
        ev = Event(df[meas.preselect(ev)].reset_index(drop=True), False, playlist, meas.minos_p_offset)
    if len(ev):
        rv, samples = meas.reco(ev, CV)
        r_idx = RB.index(*rv)
        for s, m in samples.items():
            acc.add(f"data/{s}", fill(RB.size, r_idx[m], np.ones(m.sum())))
    return acc


def _task(args):
    kind, path, playlist, scale = args
    t0 = time.time()
    acc = process_mc(path, playlist, scale) if kind == "mc" else process_data(path, playlist)
    return acc, (kind, path, time.time() - t0)


def run(meas, mc_files, data_files, workers=8, verbose=True):
    """Returns (Accumulator, info dict)."""
    tune = Tune(meas.tune_name, **getattr(meas, "tune_options", {}))
    universes = meas.universes(tune)
    _STATE.update(meas=meas, tune=tune, univ=universes)
    mc = io.file_summary(mc_files)
    da = io.file_summary(data_files)
    pot_data = da.groupby("playlist").pot.sum()
    pot_mc = mc.groupby("playlist").pot.sum()
    missing = set(pot_data.index) - set(pot_mc.index)
    if missing:
        raise ValueError(f"no MC for data playlists {missing}")
    mc = mc[mc.playlist.isin(pot_data.index)]
    scales = (pot_data / pot_mc).dropna()
    tasks = [("mc", r.path, r.playlist, float(scales[r.playlist])) for r in mc.itertuples()]
    tasks += [("data", r.path, r.playlist, 1.0) for r in da.itertuples()]
    tasks.sort(key=lambda t: t[0] != "mc")  # long MC tasks first
    acc = Accumulator()
    t0 = time.time()
    ctx = mp.get_context("fork")
    with ctx.Pool(workers) as pool:
        for i, (a, (kind, path, dt)) in enumerate(pool.imap_unordered(_task, tasks)):
            acc.merge(a)
            if verbose and (i % 10 == 0 or i == len(tasks) - 1):
                print(f"[{i + 1}/{len(tasks)}] {kind} {dt:.1f}s  elapsed {time.time() - t0:.0f}s", flush=True)
    info = dict(pot_data=pot_data.to_dict(), pot_mc=pot_mc.to_dict(), universes=universes, tune=tune,
                n_mc_files=len(mc), n_data_files=len(da))
    return acc, info
