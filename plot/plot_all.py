#!/usr/bin/env python3
"""Make all comparison figures (figures/<measurement>/).  Extra arguments (e.g. --format png) are passed on."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["plot_ccinc_2d.py", "plot_qe3d.py", "plot_qe3d_leme.py", "plot_ccpip.py"]

if __name__ == "__main__":
    failed = []
    for s in SCRIPTS:
        print(f"== {s}")
        # separate processes: several measurements' run.py modules import their own `analysis`
        if subprocess.call([sys.executable, os.path.join(HERE, s)] + sys.argv[1:]) != 0:
            failed.append(s)
    if failed:
        sys.exit(f"failed: {' '.join(failed)}")
