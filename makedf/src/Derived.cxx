#include "Derived.h"

#include <algorithm>
#include <cmath>

#include "TTree.h"

namespace makedf {

bool InHexagon(double x, double y, double apothem) {
  // regular hexagon with flat sides at |x| = apothem (MAT CCInclusiveCuts convention)
  if (std::fabs(x) >= apothem) return false;
  const double lenOfSide = apothem * 2. / std::sqrt(3.);
  return std::fabs(y) < lenOfSide - std::fabs(x) / std::sqrt(3.);
}

double Theta3D(double tx, double ty) {
  const double sx = 1. / std::cos(tx), sy = 1. / std::cos(ty);
  return std::acos(1. / std::sqrt(sx * sx + sy * sy - 1.));
}

// ---------------------------------------------------------------------------
void TruthSummary::request(BranchSet& bs) {
  fN = bs.require("mc_nFSPart");
  fPDG = bs.require("mc_FSPartPDG");
  fPx = bs.require("mc_FSPartPx");
  fPy = bs.require("mc_FSPartPy");
  fPz = bs.require("mc_FSPartPz");
  fE = bs.require("mc_FSPartE");
  fLep = bs.require("mc_primFSLepton");
  fNu = bs.require("mc_incomingPartVec");
}

void TruthSummary::book(TTree* t) {
  auto I = [&](const std::string& n, int* v) { BindOutput(t, n, v, 'I'); };
  auto F = [&](const std::string& n, float* v) { BindOutput(t, n, v, 'F'); };
  I("tfs_np", &n_p); I("tfs_nn", &n_n); I("tfs_npip", &n_pip); I("tfs_npim", &n_pim); I("tfs_npi0", &n_pi0);
  I("tfs_nkpm", &n_kpm); I("tfs_nk0", &n_k0); I("tfs_nomeson", &n_omeson); I("tfs_nhbaryon", &n_hbaryon);
  I("tfs_ngamma", &n_gamma); I("tfs_nlep", &n_lep); I("tfs_np_T120", &n_p_120); I("tfs_nnucfrag", &n_nucfrag);
  F("tfs_sumTp", &sum_Tp); F("tfs_sumTn", &sum_Tn); F("tfs_maxgammaE", &max_gammaE); F("tfs_eavail", &eavail);
  I("tfs_leadnuc_pdg", &lead_nucleon_pdg); F("tfs_leadnuc_E", &lead_nucleon_E);
  for (int k = 0; k < K; ++k) {
    std::string s = std::to_string(k);
    F("tfs_p_T_" + s, &p_T[k]); F("tfs_p_px_" + s, &p_px[k]);
    F("tfs_p_py_" + s, &p_py[k]); F("tfs_p_pz_" + s, &p_pz[k]);
    I("tfs_pi_pdg_" + s, &pi_pdg[k]); F("tfs_pi_T_" + s, &pi_T[k]);
    F("tfs_pi_px_" + s, &pi_px[k]); F("tfs_pi_py_" + s, &pi_py[k]);
    F("tfs_pi_pz_" + s, &pi_pz[k]);
  }
  // true muon kinematics in the beam frame and energy/momentum transfer [MeV]
  F("t_theta_mu", &theta_mu); F("t_p_mu", &p_mu); F("t_pt_mu", &pt_mu); F("t_pz_mu", &pz_mu);
  F("t_q0", &q0); F("t_q3", &q3);
}

void TruthSummary::compute() {
  n_p = n_n = n_pip = n_pim = n_pi0 = n_kpm = n_k0 = n_omeson = n_hbaryon = n_gamma = n_lep = n_p_120 = n_nucfrag = 0;
  sum_Tp = sum_Tn = max_gammaE = eavail = 0;
  lead_nucleon_pdg = 0;
  lead_nucleon_E = 0;
  struct P { double T, px, py, pz; int pdg; };
  std::vector<P> protons, pions;
  const size_t n = std::min<size_t>((size_t)std::lround(fN->get(0)), fPDG->size());
  for (size_t i = 0; i < n; ++i) {
    const int pdg = (int)std::lround(fPDG->get(i));
    const int apdg = std::abs(pdg);
    const double E = fE->get(i), px = fPx->get(i), py = fPy->get(i), pz = fPz->get(i);
    if (pdg == 2212) {
      const double T = std::max(E - kMp, 0.0);  // protons at rest: avoid -O(keV) rounding below 0
      ++n_p; sum_Tp += T; eavail += T;
      if (T > 120) ++n_p_120;
      protons.push_back({T, px, py, pz, pdg});
      if (E > lead_nucleon_E) { lead_nucleon_E = E; lead_nucleon_pdg = pdg; }
    } else if (pdg == 2112) {
      ++n_n; sum_Tn += std::max(E - kMn, 0.0);
      if (E > lead_nucleon_E) { lead_nucleon_E = E; lead_nucleon_pdg = pdg; }
    } else if (apdg == 211) {
      if (pdg > 0) ++n_pip; else ++n_pim;
      eavail += E - kMpi;
      pions.push_back({E - kMpi, px, py, pz, pdg});
    } else if (pdg == 111) {
      ++n_pi0; eavail += E;
    } else if (apdg == 321) {
      ++n_kpm; eavail += E;
    } else if (apdg == 311 || apdg == 130 || apdg == 310) {
      ++n_k0; eavail += E;
    } else if (pdg == 22) {
      ++n_gamma; max_gammaE = std::max<float>(max_gammaE, E); eavail += E;
    } else if (apdg >= 11 && apdg <= 16) {
      ++n_lep;
    } else if (apdg > 1000000000) {
      ++n_nucfrag;
    } else if (apdg >= 100 && apdg < 1000) {
      ++n_omeson; eavail += E;
    } else if (apdg >= 1000 && apdg < 10000) {
      // baryons other than p/n (Lambda, Sigma, charmed, ...)
      ++n_hbaryon; eavail += E - kMp;
    } else {
      eavail += E;
    }
  }
  // primary lepton counted above among leptons; remove it
  if (n_lep > 0) --n_lep;
  auto byT = [](const P& a, const P& b) { return a.T > b.T; };
  std::sort(protons.begin(), protons.end(), byT);
  std::sort(pions.begin(), pions.end(), byT);
  for (int k = 0; k < K; ++k) {
    if ((size_t)k < protons.size()) {
      p_T[k] = protons[k].T; p_px[k] = protons[k].px; p_py[k] = protons[k].py; p_pz[k] = protons[k].pz;
    } else {
      p_T[k] = p_px[k] = p_py[k] = p_pz[k] = NAN;
    }
    if ((size_t)k < pions.size()) {
      pi_pdg[k] = pions[k].pdg; pi_T[k] = pions[k].T;
      pi_px[k] = pions[k].px; pi_py[k] = pions[k].py; pi_pz[k] = pions[k].pz;
    } else {
      pi_pdg[k] = 0; pi_T[k] = pi_px[k] = pi_py[k] = pi_pz[k] = NAN;
    }
  }
  // muon in the beam frame: rotate about x by kBeamAngle
  const double lx = fLep->get(0), ly = fLep->get(1), lz = fLep->get(2), lE = fLep->get(3);
  const double c = std::cos(kBeamAngle), s = std::sin(kBeamAngle);
  const double ry = c * ly - s * lz, rz = s * ly + c * lz;
  p_mu = std::sqrt(lx * lx + ly * ly + lz * lz);
  pz_mu = rz;
  pt_mu = std::sqrt(lx * lx + ry * ry);
  theta_mu = p_mu > 0 ? std::acos(std::max(-1.0, std::min(1.0, rz / p_mu))) : NAN;
  const double nx = fNu->get(0), ny = fNu->get(1), nz = fNu->get(2), nE = fNu->get(3);
  q0 = nE - lE;
  q3 = std::sqrt((nx - lx) * (nx - lx) + (ny - ly) * (ny - ly) + (nz - lz) * (nz - lz));
}

// ---------------------------------------------------------------------------
void Preselection::requestReco(BranchSet& bs) {
  fVtx = bs.require("vtx");
  fMinos = bs.require("isMinosMatchTrack");
  fTx = bs.require("muon_thetaX");
  fTy = bs.require("muon_thetaY");
}

void Preselection::requestTruth(BranchSet& bs) {
  fMcVtx = bs.require("mc_vtx");
  fMcNu = bs.require("mc_incoming");
  fMcCur = bs.require("mc_current");
  fMcLep = bs.require("mc_primFSLepton");
}

bool Preselection::passReco() const {
  if (!fC.enabled) return true;
  const double z = fVtx->get(2);
  if (!(z > fC.zmin && z < fC.zmax)) return false;
  if (!InHexagon(fVtx->get(0), fVtx->get(1), fC.apothem)) return false;
  if (fC.requireMinos && std::lround(fMinos->get(0)) != 1) return false;
  const double th = Theta3D(fTx->get(0), fTy->get(0));
  if (!(th * 180. / M_PI < fC.maxThetaDeg)) return false;
  return true;
}

bool Preselection::passTruth() const {
  if (!fC.enabled) return true;
  const int nu = (int)std::lround(fMcNu->get(0));
  if (fC.truthCCOnly) {
    if (std::abs(nu) != 14) return false;
    if (std::lround(fMcCur->get(0)) != 1) return false;
  }
  const double z = fMcVtx->get(2);
  if (!(z > fC.zmin && z < fC.zmax)) return false;
  if (!InHexagon(fMcVtx->get(0), fMcVtx->get(1), fC.apothem)) return false;
  const double lx = fMcLep->get(0), ly = fMcLep->get(1), lz = fMcLep->get(2);
  const double c = std::cos(kBeamAngle), s = std::sin(kBeamAngle);
  const double rz = s * ly + c * lz;
  const double p = std::sqrt(lx * lx + ly * ly + lz * lz);
  if (!(p > 0 && std::acos(std::max(-1.0, std::min(1.0, rz / p))) * 180. / M_PI < fC.maxThetaDeg)) return false;
  return true;
}

}  // namespace makedf
