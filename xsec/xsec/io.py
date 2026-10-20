"""Reading makedf dataframe files."""
import glob
import os
import re

import numpy as np
import pandas as pd
import uproot


def df_files(*patterns):
    out = []
    for p in patterns:
        if os.path.isdir(p):
            p = os.path.join(p, "**", "*.root")
        out += glob.glob(p, recursive=True)
    return sorted(set(out))


def _playlist_from_path(path):
    m = re.search(r"Playlist([0-9][A-Z](?:_Water(?:Full|Empty))?)", path)
    return m.group(1) if m else None


def read_meta(path):
    """Per-input-file metadata of one DF file as a DataFrame (adds 'playlist')."""
    with uproot.open(path) as f:
        m = f["meta"].arrays(["file_id", "is_mc", "pot_used", "pot_total", "n_reco", "n_reco_pass", "n_truth",
                              "n_truth_pass", "file"], library="np")
    df = pd.DataFrame({k: v for k, v in m.items() if k != "file"})
    df["file"] = [str(x) for x in m["file"]]
    df["playlist"] = [_playlist_from_path(x) for x in df["file"]]
    df["df_file"] = path
    return df


def file_summary(paths):
    """One row per DF file: is_mc, playlist, POT."""
    rows = []
    for p in paths:
        m = read_meta(p)
        pls = set(m.playlist)
        if len(pls) != 1:
            raise ValueError(f"{p}: DF file spans playlists {pls}")
        rows.append(dict(path=p, is_mc=bool(m.is_mc.iloc[0]), playlist=pls.pop(), pot=m.pot_used.sum(),
                         n_reco=m.n_reco_pass.sum(), n_truth=m.n_truth_pass.sum()))
    return pd.DataFrame(rows)


def load(path, tree, columns=None, cut=None):
    """Load a tree into pandas.  `columns` may contain names missing from the tree (ignored)."""
    with uproot.open(path) as f:
        if tree not in f:
            return pd.DataFrame()
        t = f[tree]
        keys = set(t.keys())
        cols = [c for c in columns if c in keys] if columns else list(keys)
        arr = t.arrays(cols, library="np")
    df = pd.DataFrame(arr)
    if cut is not None:
        df = df[cut(df)].reset_index(drop=True)
    return df
