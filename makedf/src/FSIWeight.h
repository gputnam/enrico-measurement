// GENIE (v2, hA) elastic-FSI / pion-absorption "bug fix" weights.
//
// MAT-independent port of PlotUtils::weight_fsi + PlotUtils::fate_code
// (MAT-MINERvA weighters/weight_fsi.{h,cxx}) as configured by
// PlotUtils::FSIReweighter (weighters/FSIReweighter.h):
//     fCalculator.UseTrackingThreshold();
//     fCalculator.calcWeights(...);
//     w = GetElasticWeight(1) [* GetAbsorptionWeight()]
// This is the "GENIE_FSI_Bug_Fix" reweight used by MINERvA tunes vX.Y.1.
//
// Physics constants and logic are kept identical to the original, including its
// quirks (marked "QUIRK" in FSIWeight.cxx). Only GetElasticWeight(1|2) and
// GetAbsorptionWeight() are provided. The per-entry cache (fLastEntry) and the
// verbose printout of the original are dropped: every call recomputes.
#pragma once

#include <vector>

#include "Columns.h"

class TTree;

namespace makedf {
namespace fsi {

struct IntranukeParticle {
  int pdg;
  int fate;
  double ke;              // kinetic energy [MeV]
  double weight1 = -1.0;  // elastic -> no FSI
  double weight2 = -1.0;  // elastic -> other FSI
};

// GENIE v2 nuclear binding energies (UserPhysicsOptions.xml), MeV.
double GenieBEinMeV(int A);

// Port of fate_code::calcFates: status-14 p/n/pi0/pi+/pi- in the GENIE event record,
// classified by hA fate code, sorted by decreasing KE. Input arrays have nPart entries.
std::vector<IntranukeParticle> CalcFates(int mc_intType, int mc_targetA, int mc_er_nPart, const int* mc_er_ID,
                                         const int* mc_er_status, const int* mc_er_FD, const int* mc_er_LD,
                                         const double* mc_er_Px, const double* mc_er_Py, const double* mc_er_Pz,
                                         const double* mc_er_E);

// Port of weight_fsi. mc_incoming, mc_primaryLepton and mc_er_mother were arguments of the
// original but are never used in the calculation, so they are not taken here.
class WeightFSI {
 public:
  void UseTrackingThreshold() { fUseTrackingThreshold = true; }

  void calcWeights(int mc_charm, int mc_intType, int mc_targetA, int mc_targetZ, int mc_resID, int mc_er_nPart,
                   const int* mc_er_ID, const int* mc_er_status, const int* mc_er_FD, const int* mc_er_LD,
                   const double* mc_er_Px, const double* mc_er_Py, const double* mc_er_Pz, const double* mc_er_E);

  double GetElasticWeight(int config) const;  // config 1: elastic -> no FSI; 2: elastic -> other FSI
  double GetAbsorptionWeight() const { return fAbsorptionWeight; }

  int count_elastic_nucleon = 0, count_elastic_pion = 0;
  double elastic_weight_nucleon = -1.0, elastic_weight_pion = -1.0;

 private:
  void Reset();
  int GetQEElasticWeight(double& w1, double& w2, int fate, double keGeV, int A, int pdg) const;
  int GetQEElasticFSIWeightFromParam(double& w1, double& w2, double keGeV, int A, int pdg) const;

  static constexpr double m_breakpointlow = 0.110;  // GeV
  static constexpr double m_breakpointhigh = 1.01;
  static constexpr double m_highcutoff = 3.8;
  bool fUseTrackingThreshold = false;
  double fElasticWeight1 = 1.0, fElasticWeight2 = 1.0, fAbsorptionWeight = 1.0;
};

}  // namespace fsi

// Per-event FSI bug-fix weights (MC only). Works on rows of both the MasterAnaDev and the
// Truth tree (same branch names/types). Outputs:
//   w_fsi_elastic    = GetElasticWeight(1) with UseTrackingThreshold()
//   w_fsi_absorption = GetAbsorptionWeight()
// FSIReweighter(useElastic, useAbsorption) == product of the enabled ones.
class FSIWeights {
 public:
  void request(BranchSet& bs);
  void book(TTree* out);
  void compute();

 private:
  InBranch *fCharm = nullptr, *fIntType = nullptr, *fA = nullptr, *fZ = nullptr, *fResID = nullptr;
  InBranch *fN = nullptr, *fID = nullptr, *fStatus = nullptr, *fFD = nullptr, *fLD = nullptr;
  InBranch *fPx = nullptr, *fPy = nullptr, *fPz = nullptr, *fE = nullptr;
  fsi::WeightFSI fCalc;
  std::vector<int> vID, vStatus, vFD, vLD;
  std::vector<double> vPx, vPy, vPz, vE;
  float w_elastic = 1, w_absorption = 1;
};

}  // namespace makedf
