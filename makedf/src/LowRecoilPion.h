// MAT-independent port of the MINERvA "LowRecoilPion" (untracked / trackless pion)
// Michel reconstruction, as used by the ME numu CC1pi+ analysis (CC-CH-pip-ana).
//
// Sources ported (behaviour reproduced faithfully, including known quirks):
//   MAT-MINERvA  utilities/LowRecoilPionReco.{h,cxx}   Michel ctor, DoesMichelMatchVtx,
//                                                      DoesMichelMatchClus, GetPionAngle
//   MAT-MINERvA  calculators/LowRecoilPionCuts.h       hasMichel, BestMichelDistance2D,
//                                                      GetClosestMichel
//   MAT-MINERvA  calculators/LowRecoilPionFunctions.h  GetTpiFromRange, GetTrackerECALMuFuzz,
//                                                      NewEavail, GetNMichels
//   MAT-MINERvA  utilities/NSFDefaults.h               Tpi(range) fit constants
// The CC1pi analysis uses the MAT-MINERvA branch "low_recoil_michels_EverChanges"
// (tuple_idx set in the Michel ctor, MichelEvent::m_bestthetaangle); the algorithm is
// otherwise identical to MAT-MINERvA main.
//   CC-CH-pip-ana includes/Michel.h                   endpoint::GetQualityMichels (tracked)
//
// Units: mm, MeV, microseconds (as in MAT), angles in radians.
#pragma once

#include <map>
#include <string>
#include <vector>

#include "Columns.h"

class TTree;

namespace makedf {
namespace lrp {

// --- constants (MAT values) --------------------------------------------------
constexpr double kNuMIBeamAngle = -0.05887;   // MinervaUnits::numi_beam_angle_rad
constexpr double kMpiLRP = 139.57;            // mass used in Michel::GetPionAngle
constexpr double kMpionMAT = 139.5701;        // MinervaUnits::M_pion (W_exp, Ehad)
constexpr double kTpiRangeP0 = 0.128706;      // NSFDefaults::tpi_from_michel_range_fit_p0_cv
constexpr double kTpiRangeP1 = 3.42486;       // NSFDefaults::tpi_from_michel_range_fit_p1_cv
constexpr double kTpiRangeP0Err = 3.28304e-3; // ..._p0_err (TpiFromMichelRangeFit universes)
constexpr double kTpiRangeP1Err = 7.96999e-2; // ..._p1_err
constexpr double kMAX2DDIST = 150.;           // LowRecoilPion::kMAX2DDIST
constexpr double kMAX3DDIST = 2600.;          // GetClosestMichel range cut
constexpr double kEavailScale = 1.17;         // NewEavail / CVUniverse::GetEavail

// Tpi [MeV] from Michel range [mm] (CVUniverse::GetTpiFromRange, CV parameters).
double TpiFromRange(double range, double p0 = kTpiRangeP0, double p1 = kTpiRangeP1);
// CVUniverse::thetaWRTBeam
double ThetaWRTBeam(double x, double y, double z);

// Interaction vertex (MAT GetVertex(): vtx[0..3]; time converted ns -> us).
struct Vertex {
  double x = 0, y = 0, z = 0, t = 0;
};

// Input branches needed by the reconstruction.
struct Inputs {
  // reco (data + MC)
  InBranch *vtx = nullptr;
  InBranch *fitPass = nullptr, *energy = nullptr, *time = nullptr, *datafrac = nullptr;
  InBranch *x1 = nullptr, *y1 = nullptr, *u1 = nullptr, *v1 = nullptr, *z1 = nullptr;
  InBranch *x2 = nullptr, *y2 = nullptr, *u2 = nullptr, *v2 = nullptr, *z2 = nullptr;
  InBranch *clusViewSz = nullptr, *clusView = nullptr, *clusSubdet = nullptr, *clusIsMu = nullptr;
  InBranch *clusEnergy = nullptr, *clusPos = nullptr, *clusTime = nullptr, *clusZ = nullptr;
  // MC only (nullptr on data)
  InBranch *tParentPDG = nullptr, *tParentE = nullptr, *tParentP = nullptr;
  InBranch *tParentPx = nullptr, *tParentPy = nullptr, *tParentPz = nullptr;
  InBranch *tMichelInitZ = nullptr;
  bool isMC = false;

  void request(BranchSet& bs);
  Vertex vertex() const;
  int nMichelCandidates() const;  // GetNMichels(): FittedMichel_michel_fitPass_sz
};

// LowRecoilPion::Michel<T> (only members that feed any output or decision).
struct Michel {
  int tuple_idx = -1;
  double energy = -999., time = -999.;
  int is_fitted = -1;
  double m_x1 = 9999., m_y1 = 9999., m_u1 = 9999., m_v1 = 9999., m_z1 = 9999.;
  double m_x2 = 9999., m_y2 = 9999., m_u2 = 9999., m_v2 = 9999., m_z2 = 9999.;
  double overlay_fraction = -1.;
  int nclusters = -1;

  double up_to_vertex_XZ = 9999., up_to_vertex_UZ = 9999., up_to_vertex_VZ = 9999.;
  double down_to_vertex_XZ = 9999., down_to_vertex_UZ = 9999., down_to_vertex_VZ = 9999.;
  double down_to_clus_XZ = 9999., down_to_clus_UZ = 9999., down_to_clus_VZ = 9999.;
  double up_to_clus_XZ = 9999., up_to_clus_UZ = 9999., up_to_clus_VZ = 9999.;
  double up_to_vertex_dist3D = 9999., down_to_vertex_dist3D = 9999.;
  double up_to_cluster_dist3D = -9999., down_to_cluster_dist3D = -9999.;
  double up_clus_michel_dist3D = -9999., down_clus_michel_dist3D = -9999.;
  double up_clus_michvtx_dist3D = -9999., down_clus_michvtx_dist3D = 9999.;
  double vtx_michel_timediff = 9999.;
  int vtx_endpoint = 0, clus_endpoint = 0;

  int BestMatch = 0, SecondBestMatch = 0;  // 1 UpVtx, 2 DownVtx, 3 UpClus, 4 DownClus
  double Best3Ddist = 9999.;
  double best_XZ = -9999., best_UZ = -9999., best_VZ = -9999.;
  int bestview = -1, recoEndpoint = -1, OrderOfMichel = -1;
  int passable_matchtype[4] = {-1, -1, -1, -1};

  // GetPionAngle
  double best_angle = -9999., best_phi = -9999.;
  double reco_Epi = -9999., reco_KE = -9999., reco_ppi = -9999.;
  double reco_ppix = -9999., reco_ppiy = -9999., reco_ppiz = -9999.;

  // truth (MC only)
  double true_initialz = -9999., true_angle = -9999., pionKE = 99999.;
  int true_parentpdg = -1, trueEndpoint = -1;

  void Init(const Inputs& in, int ci, bool fillTruth);   // ctor body w/o matching
  void DoesMichelMatchVtx(const Vertex& vtx);
  void DoesMichelMatchClus(const Inputs& in, const Vertex& vtx);
  void GetPionAngle(const Vertex& vtx);
};

// LowRecoilPion::MichelEvent<T> (members used downstream).
struct MichelEvent {
  int m_idx = -1;
  double m_bestdist = 9999.;
  double m_bestthetaangle = 9999.;
  double m_best_XZ = 9999., m_best_UZ = 9999., m_best_VZ = 9999.;
  int m_matchtype = 0;
  std::vector<Michel> m_nmichels, m_nmichelspass, m_allmichels;
  double nclusters = -1;
  double best_x = 9999., best_y = 9999., best_z = 9999.;
  double lowTpi = 9999.;
  int selection = 0;
  int ClosestMichelsIsOverlay = -1;
};

// The three MAT cuts (static XXXCut functions); same return values and side effects.
bool HasMichelCut(const Inputs& in, MichelEvent& evt, bool fillTruth);
bool BestMichelDistance2DCut(MichelEvent& evt, double max_dist = kMAX2DDIST);
bool GetClosestMichelCut(const Inputs& in, MichelEvent& evt, double max_dist = kMAX3DDIST);

// MuonFuzz sums (GetTrackerECALMuFuzz) [MeV]
void TrackerECALMuFuzz(const InBranch* planeIDs, const InBranch* energies, double& trk, double& ecal);

// --- tracked-pion (endpoint) Michels: CC-CH-pip-ana includes/Michel.h -------------
namespace endpoint {
enum EMatchCategory { kNoMatch, kFit, kNoFit, kOV, kNMatchCategories };
struct Michel {
  int idx = -105, vtx = -106, had_idx = -107;
  EMatchCategory match_category = kNoMatch;
  double fit_distance = -1.;
};
struct Inputs {
  InBranch *idx = nullptr, *endDist = nullptr, *avgDist = nullptr, *ovDist = nullptr;
  void request(BranchSet& bs);
};
// GetQualityMichels: map michel cluster idx -> Michel (one per hadron track endpoint).
std::map<int, Michel> GetQualityMichels(const Inputs& in);
}  // namespace endpoint

}  // namespace lrp

// Per-event summary columns (prefix lrp_). Follows the TruthSummary pattern.
class LowRecoilPionSummary {
 public:
  void request(BranchSet& bs);  // MC-only columns are booked iff truth branches exist
  void book(TTree* out);
  void compute();
  bool isMC() const { return fIn.isMC; }
  const lrp::MichelEvent& event() const { return fEvt; }

 private:
  lrp::Inputs fIn;
  lrp::endpoint::Inputs fEp;
  InBranch *fFuzzID = nullptr, *fFuzzE = nullptr, *fRecTrk = nullptr, *fRecEcal = nullptr;
  lrp::MichelEvent fEvt;

  // untracked (data + MC)
  int ncand, nfitted, nclusters, nall, npass2d, pass_has, pass_2d, pass_closest, good;
  int idx, tuple_idx, matchtype, bestview;
  float range, tpi, theta, phi, best_x, best_y, best_z, michel_e, michel_dt, ovfrac;
  // event-level Eavail ingredients
  float mufuzz_trk, mufuzz_ecal, eavail;
  // tracked endpoint michels
  int trk_nmich, trk_hadmask;
  // MC-only Michel truth (closest Michel)
  int t_parentpdg, t_endpoint;
  float t_parent_tpi, t_parent_theta;
};

}  // namespace makedf
