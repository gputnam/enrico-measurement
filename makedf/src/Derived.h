// Derived per-event quantities computed in C++ (variable-length inputs that would
// be expensive to store), plus the preselection.
#pragma once

#include <string>
#include <vector>

#include "Columns.h"

class TTree;

namespace makedf {

constexpr double kBeamAngle = -0.05887;  // NuMI beam pitch [rad], rotation about x (MAT convention)
constexpr double kMp = 938.272, kMn = 939.565, kMpi = 139.570, kMmu = 105.658;

struct PreselConfig {
  // Reco: loose superset of the tracker-fiducial analyses
  double zmin = 5800, zmax = 8600, apothem = 900, maxThetaDeg = 25;
  bool requireMinos = true;
  // Truth: CC numu/numubar with true vertex in loose FV and true theta_mu < maxThetaDeg
  bool truthCCOnly = true;
  bool enabled = true;
};

bool InHexagon(double x, double y, double apothem);
// Angle w.r.t. the beam from MINERvA thetaX/thetaY (MAT theta3D)
double Theta3D(double thetaX, double thetaY);

// Final-state truth summary from mc_FSPart* arrays (MC only).
class TruthSummary {
 public:
  void request(BranchSet& bs);
  void book(TTree* out);  // creates or re-binds output branches
  void compute();

 private:
  InBranch *fN = nullptr, *fPDG = nullptr, *fPx = nullptr, *fPy = nullptr, *fPz = nullptr, *fE = nullptr;
  InBranch *fLep = nullptr, *fNu = nullptr;
  static constexpr int K = 2;  // leading particles kept
  int n_p, n_n, n_pip, n_pim, n_pi0, n_kpm, n_k0, n_omeson, n_hbaryon, n_gamma, n_lep, n_p_120, n_nucfrag;
  float sum_Tp, sum_Tn, max_gammaE, eavail, lead_nucleon_E;
  int lead_nucleon_pdg;
  float p_T[K], p_px[K], p_py[K], p_pz[K];
  float pi_T[K], pi_px[K], pi_py[K], pi_pz[K];
  int pi_pdg[K];
  float theta_mu, p_mu, pt_mu, pz_mu, q0, q3;
};

class Preselection {
 public:
  Preselection(const PreselConfig& c) : fC(c) {}
  void requestReco(BranchSet& bs);
  void requestTruth(BranchSet& bs);
  bool passReco() const;
  bool passTruth() const;

 private:
  PreselConfig fC;
  InBranch *fVtx = nullptr, *fMinos = nullptr, *fTx = nullptr, *fTy = nullptr;
  InBranch *fMcVtx = nullptr, *fMcNu = nullptr, *fMcCur = nullptr, *fMcLep = nullptr;
};

}  // namespace makedf
