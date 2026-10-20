"""Parse the CC1pi+ (arXiv:2605.24224) cross-section and covariance tables from the TeX sources.

xsec tables: label tbl:CrossSection<var>, units 1e-42 cm^2/<unit>/nucleon.
covariances (CovarianceTables.tex): tbl:{stat,sys}cov_<var>_scintillator, units 1e-84.
"""
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def _tables(path):
    txt = open(path).read()
    for m in re.finditer(r"\\begin\{table\}.*?\\end\{table\}", txt, re.S):
        block = m.group(0)
        lab = re.search(r"\\label\{([^}]*)\}", block)
        yield (lab.group(1) if lab else None), block


def _rows(block):
    body = block.split("\\hline")[-1].split("\\end{tabular}")[0]
    for line in body.split("\\\\"):
        line = line.strip()
        if not line or line.startswith("\\end"):
            continue
        cells = [c.strip() for c in line.split("&")]
        yield cells


def _edges(label):
    lo, hi = re.split(r"\s*-{1,2}\s*", label.strip(), maxsplit=1)
    return float(lo), float(hi)


def load_xsec():
    out = {}
    for lab, block in _tables(os.path.join(HERE, "reference", "main.tex")):
        if not lab or not lab.startswith("tbl:CrossSection"):
            continue
        var = lab[len("tbl:CrossSection"):]
        rows = []
        for c in _rows(block):
            if len(c) < 6:
                continue
            lo, hi = _edges(c[0])
            rows.append([lo, hi] + [float(x) for x in c[1:6]])
        out[var] = np.array(rows)  # lo, hi, xsec, abs_unc, frac_stat, frac_sys, frac_flux
    return out


def load_cov():
    out = {}
    for lab, block in _tables(os.path.join(HERE, "reference", "CovarianceTables.tex")):
        if not lab:
            continue
        m = re.match(r"tbl:(stat|sys)cov_(\w+)_scintillator", lab)
        if not m:
            continue
        rows = [[float(x) for x in c[1:]] for c in _rows(block) if len(c) > 1]
        out[(m.group(2), m.group(1))] = np.array(rows)
    return out


if __name__ == "__main__":
    xs = load_xsec()
    for k, v in xs.items():
        print(k, v.shape, "edges", list(v[:, 0]) + [v[-1, 1]])
    for k, v in load_cov().items():
        print(k, v.shape)
