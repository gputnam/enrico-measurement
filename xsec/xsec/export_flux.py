"""Export the measurement flux (data-POT-weighted CV over playlists) for a generator.

The cross sections are normalised by this flux integrated over [0, 100] GeV (extract.flux_integrals),
so a generator prediction must use the same flux shape over the same range.  Writes a CSV in the
nuflow `csv` histogram format (BinLow(GeV),BinHigh(GeV),Flux(neutrinos/m^2/POT/GeV), density per_GeV).

    python -m xsec.export_flux -o flux.csv [--tune v1] [--pdg 14] [--emax 100]
"""
import argparse

import numpy as np

from .tunes import Tune

# ME FHC data POT per playlist of the v1 DFs (identical for all ME FHC measurements)
ME_FHC_POT = {'1A': 8.968946226551533e+19, '1B': 1.8673986975553036e+19, '1C': 4.294376238504553e+19,
              '1D': 1.4404452962278605e+20, '1E': 1.0296200958633948e+20, '1F': 1.6685546978706168e+20,
              '1G': 1.3756941427380455e+20, '1L': 1.3370718209243894e+19, '1M': 1.581268994168384e+20,
              '1N': 1.0656505761839271e+20, '1O': 2.979662817484796e+19, '1P': 4.675822390412009e+19}


def measurement_flux(tune_name="v1", pot_by_playlist=None, pdg_sign=1):
    """(edges [GeV], flux [nu/m^2/POT/GeV]) POT-weighted over playlists."""
    pot = pot_by_playlist or ME_FHC_POT
    tune = Tune(tune_name)
    tot = sum(pot.values())
    edges, cv = None, 0.0
    for pl, p in pot.items():
        f = tune.flux_for(pl, pdg_sign)
        if edges is None:
            edges = f.edges
        assert np.allclose(edges, f.edges)
        cv = cv + f.cv * p / tot
    return np.asarray(edges, float), np.asarray(cv, float)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--tune", default="v1")
    ap.add_argument("--emax", type=float, default=100.0, help="last bin kept: the one containing emax")
    args = ap.parse_args()
    e, cv = measurement_flux(args.tune)
    b1 = min(np.searchsorted(e, args.emax, side="right") - 1, len(e) - 2)
    with open(args.output, "w") as f:
        f.write("BinLow(GeV),BinHigh(GeV),Flux(neutrinos/m^2/POT/GeV)\n")
        for i in range(b1 + 1):
            f.write(f"{e[i]:.6e},{e[i + 1]:.6e},{cv[i]:.6e}\n")
    integ = (cv[:b1 + 1] * np.diff(e)[:b1 + 1]).sum()
    print(f"wrote {args.output}: {b1 + 1} bins, [{e[0]}, {e[b1 + 1]}] GeV, integral {integ:.5e} nu/m^2/POT")


if __name__ == "__main__":
    main()
