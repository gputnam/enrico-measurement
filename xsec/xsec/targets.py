"""Number of scattering targets in the MINERvA tracker (port of MAT TargetUtils)."""
import numpy as np

from .paths import aux

AVOGADRO = 6.0221412927e23
AREAL_DENSITY = 1.9872 / 100.0  # g/mm^2 (scintillator plane, NX MC and data)
ATOMIC_MASS = dict(C=12.0107, H=1.00794, O=15.9994, Ti=47.867, Al=26.982, Si=28.0855, Cl=35.453, N=14.007)
PROTONS = dict(C=6, H=1, O=8, Ti=22, Al=13, Si=14, Cl=17, N=7)
MASS_FRACTION = {
    "NX_MC": dict(C=0.8896, H=0.07533, O=0.02432, Ti=0.00697, Al=0.001613, Si=0.001613, Cl=0.000062, N=0.000651),
    "Eroica_MC": dict(C=0.8784, H=0.0819, O=0.0286, Ti=0.0046, Al=0.0013, Si=0.0014, Cl=0.0038, N=0.0),
    "data": dict(C=0.8851, H=0.0818, O=0.0250, Ti=0.0047, Al=0.0007, Si=0.0007, Cl=0.0020, N=0.0),
}
TARGET_MASS_ERR_CH = 0.014  # fractional (MAT NSFDefaults::ch_err)


def plane_z():
    z = []
    for line in open(aux("MATFluxAndReweightFiles", "Minerva_Planes.txt")):
        p = line.split()
        if len(p) >= 6:
            z.append(float(p[5]))
    return np.array(z)


def n_planes(zmin, zmax):
    z = plane_z()
    return int(((z >= zmin) & (z <= zmax)).sum())


def hex_area(apothem=850.0):
    return apothem ** 2 * 6.0 / np.sqrt(3.0)


def tracker_mass(zmin, zmax, apothem=850.0):
    """grams"""
    return n_planes(zmin, zmax) * hex_area(apothem) * AREAL_DENSITY


def tracker_nucleons(zmin=5980.0, zmax=8422.0, apothem=850.0, fractions="NX_MC"):
    mass = tracker_mass(zmin, zmax, apothem)
    f = MASS_FRACTION[fractions]
    # protons + neutrons, with neutrons = A - Z per atom (MAT convention)
    return sum(mass * f[e] * AVOGADRO / ATOMIC_MASS[e] * ATOMIC_MASS[e] for e in f)


def tracker_atoms(element, zmin=5980.0, zmax=8422.0, apothem=850.0, fractions="NX_MC"):
    mass = tracker_mass(zmin, zmax, apothem)
    return mass * MASS_FRACTION[fractions][element] * AVOGADRO / ATOMIC_MASS[element]
