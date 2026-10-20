"""MINERvA open-data cross-section framework (python stage).

Modules
-------
paths      : locations of auxiliary inputs (flux / reweight npz)
io         : reading makedf dataframe files
flux       : flux model (CV reweight, nu-e constraint, universes, integrals)
reweight   : model weights (2p2h, RPA, MINOS efficiency, non-res pion, ...)
tunes      : MINERvA tune compositions
kinematics : muon / hadron reconstruction quantities with lateral shifts
universes  : systematic universe registry
binning    : N-d binning with linearized global bins
histo      : universe-aware histogram accumulation
unfold     : D'Agostini unfolding with analytic error propagation
extract    : background subtraction, unfolding, efficiency, normalization
targets    : number of target nucleons in the tracker fiducial volume
"""
