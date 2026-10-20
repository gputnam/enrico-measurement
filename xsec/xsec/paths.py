"""Locations of auxiliary inputs."""
import functools
import os

import numpy as np

AUX = os.environ.get("MVA_AUX", "/exp/sbnd/data/users/gputnam/MvA-Loop/aux")


def aux(*p):
    return os.path.join(AUX, *p)


@functools.lru_cache(maxsize=None)
def npz(*p):
    """Load (and cache) an npz under aux/npz as a plain dict of arrays."""
    with np.load(aux("npz", *p)) as d:
        return {k: d[k] for k in d.files}
