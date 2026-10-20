// MAT-independent port of MINERvA's Geant4 hadron-interaction reweighting:
//   MAT-MINERvA weighters/MnvHadronReweight.{h,cxx}            (the physics)
//   MAT-MINERvA universes/GeantHadronSystematics.{h,cxx}        (GEANT_Proton/Neutron/Pion universes)
//   MAT-MINERvA weighters/GeantNeutronCVReweighter.h            (neutron CV weight)
//   MAT-MINERvA universes/MinervaUniverse.cxx SetupMHRWeighter  (default configuration)
//
// Configuration reproduced (MAT defaults as used by MINERvA-101 / GetGeantHadronSystematicsMap):
//   useDefaultParticles()   : pi (pi+ and pi- treated together), p, n
//   setBasicFiducial(NSFDefaults::TrackerFace 5991.29, TrackerBack 8408.91, StandardApothem 850)
//   elastic reweight ON for neutrons, useReweightedNeutronCV() ON (=> HD neutron cross sections)
//   kinematic renormalization from mhrwKineRenorm/
//     Renorm_Kine_Truth_Minerva<PLAYLIST>_v10_HD_withElastics_RwgtNeutronCV_z5991-8408_a850.root
//
// Universe weights are returned as ratios to CV exactly as GeantHadronUniverse::GetWeightRatioToCV()
// (weightUp/Down[pdg] if the event "has" that pdg, else 1). The neutron CV weight is
// MnvHadronReweight::reweightNeutronCV (GeantNeutronCVReweighter::GetWeight).
//
// Truth vs reco rows in MAT: MnvHadronReweight itself only ever reads the truth_hadronReweight*
// and mc_* branches of the row it is given (always through the same "data" branch container), so
// the weights are identical for a Truth-tree row and a MasterAnaDev row of the same interaction.
// The only place MAT distinguishes (BaseUniverse::m_is_truth, set by SetTruth()) is the legacy
// MinervaUniverse::GetGeantHadronWeight() CV weight, which returns 1 for truth rows. Since
// 10/2024 the Reweighter-based GeantNeutronCVReweighter and the GEANT universes do NOT check
// IsTruth(). GeantWeights reproduces either convention through its constructor flags.
#pragma once

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Columns.h"

class TGraph;
class TF1;
class TH2D;
class TTree;

namespace makedf {

// Directory holding MAT's $PLOTUTILSROOT/data contents (hadronReweight/, mhrwKineRenorm/).
// auxDir is resolved as: argument if non-empty, else $MAKEDF_AUX, else <configDir>/../aux.
// If <auxDir>/MATFluxAndReweightFiles/hadronReweight exists that subdirectory is used,
// otherwise <auxDir> itself must contain hadronReweight/ and mhrwKineRenorm/.
std::string ResolveHadronReweightDataDir(const std::string& auxDir, const std::string& configDir = "");

struct HadronReweightConfig {
  std::string dataDir;                // see ResolveHadronReweightDataDir
  double minZ = 5991.29;              // NSFDefaults::TrackerFace
  double maxZ = 8408.91;              // NSFDefaults::TrackerBack
  double apothem = 850.;              // NSFDefaults::StandardApothem
  bool useNeutronCVReweight = true;   // MinervaUniverse::m_use_mhrw_neutronCV_reweight default
  bool useElastics = true;            // MinervaUniverse::m_use_mhrw_elastics default
  std::string projectName = "Renorm_Kine_Truth";
};

// One event's inputs (what MAT's HadronBranchContainer::FillUnivEntry copies).
struct HadronEvent {
  bool hasBranches = true;  // false -> MAT "default mode": all weights 1
  int nPoints = 0;          // truth_hadronReweightNPoints
  std::vector<int> intCodePerSegment, nuke, trackID, pdg;  // pdg is per track (indexed by trackID)
  std::vector<double> density, finalE, initialE, posX, posY, posZ;
  int targetZ = 0;
  double vtx[3] = {0, 0, 0};
  std::vector<int> fsPDG;
  std::vector<double> fsE, fsPx, fsPy, fsPz;
};

struct HadronWeights {
  // Same semantics as MAT's InelXSecWeights (std::map keyed by signed pdg, -211 included).
  std::map<int, double> weightUp, weightDown;
  std::map<int, bool> eventHas;
  double up(int pdg) const;
  double down(int pdg) const;
  bool has(int pdg) const;
  // GeantHadronUniverse::GetWeightRatioToCV for nsigma = +-1
  double ratio(int pdg, double nsigma) const {
    return has(pdg) ? (nsigma < 0 ? down(pdg) : up(pdg)) : 1.;
  }
};

class HadronReweight {
 public:
  explicit HadronReweight(const HadronReweightConfig& c);
  ~HadronReweight();
  HadronReweight(const HadronReweight&) = delete;
  HadronReweight& operator=(const HadronReweight&) = delete;

  // Select the playlist (MAT string e.g. "minervame1M", "me1M" or "ME1M"); loads (and caches) the
  // kinematic renormalization file for it. Must be called before getWeights/reweightNeutronCV.
  void setPlaylist(const std::string& playlist);
  const std::string& playlist() const { return fPlaylist; }

  HadronWeights getWeights(const HadronEvent& b) const;      // MnvHadronReweight::getWeights
  double reweightNeutronCV(const HadronEvent& b) const;       // MnvHadronReweight::reweightNeutronCV

  // MINERvA-101 util::GetPlaylist mapping from mc_run (or ev_run) to playlist name.
  static std::string PlaylistFromRun(int run);
  // MAT setPlaylist normalisation: strip leading "minerva", upper case.
  static std::string NormalizePlaylist(const std::string& playlist);

  std::string renormFileName(const std::string& normPlaylist) const;  // MnvHadronReweight::makekinefilename
  // Every aux file this configuration reads, for the given playlists (normalised or not).
  std::vector<std::string> requiredFiles(const std::vector<std::string>& playlists) const;

  // Diagnostics
  long nSegmentsOutsideArrays() const { return fBadIndex; }

 private:
  enum XSecType { kInelastic, kElastic };
  struct Renorm;  // per-playlist kinematic renorm histograms

  HadronReweightConfig fC;
  std::string fPlaylist;
  std::map<std::string, std::unique_ptr<Renorm>> fRenorm;
  const Renorm* fCurRenorm = nullptr;
  std::map<int, bool> fParticles;  // pdg -> treat antiparticle same

  std::vector<std::unique_ptr<TGraph>> fGraphs;
  std::map<std::string, std::unique_ptr<TF1>> fFuncs;
  mutable long fBadIndex = 0;

  TGraph* loadGraph(const std::string& file, const char* name);
  void buildFunctions();
  const TF1* getXSecFunc(int pdg, int nuke, XSecType t) const;
  double safeEval(const TF1* f, double Ti, double Tf) const;
  double getAverageCrossSection2(int pdg, int nuke, double Ti, double Tf, XSecType t, bool isCVCalc = false) const;
  double getAverageNewNeutronCrossSection(int nuke, double Ti, double Tf, XSecType t) const;

  bool inFiducial(double x, double y, double z) const;
  double getDelta(int pdg, int target, double ke) const;
  double getDeltaElastic(int pdg, int nuke, double ke) const;
  bool getWeight(int nuke, int pdg, int intCode, double density, double Ei, double Ef, double delta,
                 double& up, double& down) const;
  bool getWeightWithElastic(int nuke, int pdg, int intCode, double density, double Ei, double Ef, double dInel,
                            double dEl, double& up, double& down) const;
  bool getWeightToNewNeutronXSec(int nuke, int intCode, double density, double Ei, double Ef, double& wgt) const;

  struct Leading { double p, theta; };
  std::map<int, Leading> getLeadingFSParticles(const HadronEvent& b) const;
  std::pair<int, int> getTargetNuclei(const HadronEvent& b) const;
  const TH2D* renormHist(int tar, int nuc, int pdg, int shift) const;
};

// Per-row GEANT hadron weights (follows the TruthSummary request/book/compute pattern).
// Columns (Float_t): geant_p_dn, geant_p_up, geant_n_dn, geant_n_up, geant_pi_dn, geant_pi_up
// (universe weight ratios to CV) and geant_neutron_cv (neutron CV weight).
//
//   isTruth      : row comes from the Truth tree (MAT BaseUniverse::SetTruth(true))
//   legacyTruthCV: if true, reproduce the legacy MinervaUniverse::GetGeantHadronWeight convention
//                  (neutron CV weight = 1 on truth rows). Default false = current
//                  GeantNeutronCVReweighter (weight applied to truth rows too, i.e. to the
//                  efficiency denominator). The universe ratios never depend on isTruth in MAT.
//   playlist     : if empty, inferred per row from mc_run (as MINERvA-101 util::GetPlaylist).
class GeantWeights {
 public:
  GeantWeights(HadronReweight& hr, bool isTruth, bool legacyTruthCV = false, const std::string& playlist = "");
  void request(BranchSet& bs);
  void book(TTree* out);
  void compute();

 private:
  HadronReweight& fHR;
  bool fIsTruth, fLegacyTruthCV;
  std::string fFixedPlaylist;
  int fLastRun = -1;
  HadronEvent fEv;
  InBranch *fNPoints = nullptr, *fIntCode = nullptr, *fNuke = nullptr, *fTrackID = nullptr, *fPDG = nullptr;
  InBranch *fDens = nullptr, *fFinalE = nullptr, *fInitE = nullptr, *fX = nullptr, *fY = nullptr, *fZ = nullptr;
  InBranch *fTargetZ = nullptr, *fVtx = nullptr, *fRun = nullptr;
  InBranch *fNFS = nullptr, *fFSPDG = nullptr, *fFSE = nullptr, *fFSPx = nullptr, *fFSPy = nullptr, *fFSPz = nullptr;
  float p_dn, p_up, n_dn, n_up, pi_dn, pi_up, neutron_cv;
};

}  // namespace makedf
