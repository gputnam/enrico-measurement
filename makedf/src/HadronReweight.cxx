// Port of MAT-MINERvA MnvHadronReweight (+ GeantHadronSystematics / GeantNeutronCVReweighter).
// Function names and control flow follow the original so the two can be compared line by line;
// comments marked "MAT:" document intentional reproductions of quirks of the original.
#include "HadronReweight.h"

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <tuple>

#include "Math/IntegratorOptions.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"

namespace makedf {

namespace {

constexpr double pion_mass = 139.570, neutron_mass = 939.565, proton_mass = 938.272, kaon_mass = 493.677;
constexpr int PDGPION = 211, PDGPROTON = 2212, PDGKAON = 321, PDGNEUTRON = 2112, PDGNEUTRONCV = 2112 + 10000;
constexpr int UP = 0, DOWN = 1, CV = 2;
constexpr int minKineEntries = 3;
const char* const mhrw_shift_names[3] = {"up", "down", "cv"};
constexpr int mhrw_nTars = 12;
constexpr int mhrw_targets[12] = {1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7};
constexpr int mhrw_nuclei[12] = {26, 82, 26, 82, 6, 26, 82, 82, 26, 82, 6, 8};
constexpr int mhrw_nNTR_Scint = 7;
constexpr double numi_beam_angle_rad = -0.05887;  // MinervaUnits::numi_beam_angle_rad
constexpr int AL = 4096;                            // MAT HadronBranchContainer array length

// XSecFunctions defaults
constexpr double DEFAULT_MIN = 100, DEFAULT_MAX = 5000;
constexpr double DEFAULT_MIN_ALT = 40, DEFAULT_MAX_ALT = 2000;

// TargetUtils / TargetProp (MC geometry)
namespace TP {
constexpr double offset_pb_fe = 205.0;
constexpr double TrackerFace = 5991.29, TrackerBack = 8408.91;
constexpr double WaterFace = 5200.0, WaterBack = 5420.0;
constexpr double T1Fe = 25.75, T1Pb = 25.75, T2Fe = 25.75, T2Pb = 25.75, T3C = 76.3, T3Fe = 25.75, T3Pb = 25.75,
                 T4Pb = 8., T5Fe = 13., T5Pb = 13.;
// useCenterZMCOffset_ = true
constexpr double C1 = 4478.016 + 4.3, C2 = 4698.824 + 3.5, C3 = 4915.595 + 5.1, C3C = 4940.82 + 5.1,
                 C4 = 5641.68 + 0.0, C5 = 5774.32 + 2.0;
}  // namespace TP

// MAT: TargetUtils::IsInHexagon / MnvHadronReweight::isInHexagon, including the
// "yp*yp < lenOfSide/2" condition of the original.
bool IsInHexagonMAT(double x, double y, double apothem) {
  double lenOfSide = apothem * (2 / sqrt(3));
  double slope = (lenOfSide / 2.0) / apothem;
  double xp = fabs(x);
  double yp = fabs(y);
  if ((xp * xp + yp * yp) < apothem * apothem) return true;
  else if (xp <= apothem && yp * yp < lenOfSide / 2.0) return true;
  else if (xp <= apothem && yp < lenOfSide - xp * slope) return true;
  return false;
}

double CoordU(double x, double y) { return -x * cos(TMath::Pi() / 6) + y * sin(TMath::Pi() / 6); }
double CoordD(double x, double y) { return x * cos(TMath::Pi() / 6) + y * sin(TMath::Pi() / 6); }
double CoordC(double x, double y) { return x * sin(TMath::Pi() / 6) + y * cos(TMath::Pi() / 6); }
bool InZ(double z, double center, double thick) { return fabs(z - center) <= thick / 2; }

bool InIron1VolMC(double x, double y, double z) {
  if (!InZ(z, TP::C1, TP::T1Fe) || !IsInHexagonMAT(x, y, 850.)) return false;
  return CoordU(x, y) - TP::offset_pb_fe < 0.0;
}
bool InLead1VolMC(double x, double y, double z) {
  if (!InZ(z, TP::C1, TP::T1Pb) || !IsInHexagonMAT(x, y, 850.)) return false;
  return CoordU(x, y) - TP::offset_pb_fe >= 0.0;
}
bool InIron2VolMC(double x, double y, double z) {
  if (!InZ(z, TP::C2, TP::T2Fe) || !IsInHexagonMAT(x, y, 850.)) return false;
  return CoordD(x, y) - TP::offset_pb_fe < 0.0;
}
bool InLead2VolMC(double x, double y, double z) {
  if (!InZ(z, TP::C2, TP::T2Pb) || !IsInHexagonMAT(x, y, 850.)) return false;
  return CoordD(x, y) - TP::offset_pb_fe >= 0.0;
}
bool InCarbon3VolMC(double x, double y, double z) {
  if (!InZ(z, TP::C3C, TP::T3C) || !IsInHexagonMAT(x, y, 850.)) return false;
  return CoordC(x, y) >= 0.0;
}
bool InIron3VolMC(double x, double y, double z) {
  if (!InZ(z, TP::C3, TP::T3Fe) || !IsInHexagonMAT(x, y, 850.)) return false;
  return CoordC(x, y) < 0.0 && x < 0.0;
}
bool InLead3VolMC(double x, double y, double z) {
  if (!InZ(z, TP::C3, TP::T3Pb) || !IsInHexagonMAT(x, y, 850.)) return false;
  return CoordC(x, y) < 0.0 && x >= 0.0;
}
bool InLead4VolMC(double x, double y, double z) { return InZ(z, TP::C4, TP::T4Pb) && IsInHexagonMAT(x, y, 850.); }
bool InIron5VolMC(double x, double y, double z) {
  if (!InZ(z, TP::C5, TP::T5Fe) || !IsInHexagonMAT(x, y, 850.)) return false;
  return CoordU(x, y) - TP::offset_pb_fe < 0.0;
}
bool InLead5VolMC(double x, double y, double z) {
  if (!InZ(z, TP::C5, TP::T5Pb) || !IsInHexagonMAT(x, y, 850.)) return false;
  return CoordU(x, y) - TP::offset_pb_fe >= 0.0;
}
bool InWaterTargetVolMC(double x, double y, double z) {
  if (z < TP::WaterFace || TP::WaterBack < z) return false;
  return IsInHexagonMAT(x, y, 850.);
}
bool InTrackerMAT(double x, double y, double z) {
  if (z < TP::TrackerFace || TP::TrackerBack < z) return false;
  return IsInHexagonMAT(x, y, 850.);
}

double getMass(int pdg) {
  switch (std::abs(pdg)) {
    case PDGPION: return pion_mass;
    case PDGPROTON: return proton_mass;
    case PDGKAON: return kaon_mass;
    case PDGNEUTRON: return neutron_mass;
  }
  throw std::runtime_error("HadronReweight: getMass: unknown pdg " + std::to_string(pdg));
}

double getCorrectedDensity(int nuke, double density) {
  double c = 0;
  switch (nuke) {
    case -6: c = 4.626e22; break;
    case 6: c = 5.01398e22; break;
    case 26: c = 1.0784e22; break;
    case 82: c = 2.906e21; break;
    case 13: c = 2.2319486e22; break;
    case 8: c = 3.3428e22; break;
    default: throw std::runtime_error("HadronReweight: molecules per g unknown for nuke=" + std::to_string(nuke));
  }
  return density * c;
}

double doWeightCalcInteracted(double rhoX, double sp_data, double sp_geant, double st_data, double st_geant) {
  double denom = (1.0 - exp(-1.0 * rhoX * st_geant));
  if (denom <= 0) return 0.0;
  double nom = (1.0 - exp(-1.0 * rhoX * st_data));
  double a = sp_data / st_data;
  double b = sp_geant / st_geant;
  return (nom / denom) * a / b;
}
double doWeightCalcInelasticOnly(double rhoX, double si_data, double si_geant) {
  double denom = (1.0 - exp(-1.0 * rhoX * si_geant));
  if (denom <= 0) return 0.0;
  double nom = (1.0 - exp(-1.0 * rhoX * si_data));
  return (nom / denom);
}
double doWeightCalcNoninteracted(double rhoX, double st_data, double st_geant) {
  return exp(-1.0 * rhoX * (st_data - st_geant));
}

double NeutronReweightAmount(double neutronKE) {
  return std::max(0.0, std::min(150.0, -150 + neutronKE) / 1500.) +
         std::min(0.25, 0.1 + (3 * std::max(0.0, 50 - neutronKE)) / 500.);
}

// XSecFunctions::tgrapheval
double tgrapheval(const TGraph* graph, double x, double min, double max, bool lerpzero) {
  if (x > max) x = max;
  if (x < min) {
    if (lerpzero) {
      double foo = graph->Eval(min) * 1.0e-27;
      if (x < 0.5) x = 0.5;
      return foo * x / min;
    }
    return graph->Eval(min) * 1.0e-27;  // changeLowEElasticXSec is never set (-1)
  }
  return graph->Eval(x) * 1.0e-27;
}
// XSecFunctions::tgrapheval_withHD with useHDXSec = true (always the case after
// MnvHadronReweight::tryLoadingKineFile, which calls useHDXSec(true) unconditionally)
double tgrapheval_withHD(const TGraph* graph, double x, double /*min*/, double max) {
  if (x > max) x = max;
  return graph->Eval(x) * 1.0e-27;
}

bool FileExists(const std::string& p) {
  struct stat st;
  return stat(p.c_str(), &st) == 0;
}

}  // namespace

// ---------------------------------------------------------------------------
std::string ResolveHadronReweightDataDir(const std::string& auxDir, const std::string& configDir) {
  std::string aux = auxDir;
  if (aux.empty()) {
    if (const char* e = std::getenv("MAKEDF_AUX")) aux = e;
  }
  if (aux.empty()) {
    std::string cd = configDir;
    while (cd.size() > 1 && cd.back() == '/') cd.pop_back();
    auto sl = cd.find_last_of('/');
    aux = (sl == std::string::npos ? std::string(".") : cd.substr(0, sl)) + "/aux";
  }
  if (FileExists(aux + "/MATFluxAndReweightFiles/hadronReweight")) return aux + "/MATFluxAndReweightFiles";
  if (FileExists(aux + "/hadronReweight")) return aux;
  throw std::runtime_error("HadronReweight: no hadronReweight/ under " + aux + " (set MAKEDF_AUX)");
}

double HadronWeights::up(int pdg) const {
  auto it = weightUp.find(pdg);
  return it == weightUp.end() ? 0. : it->second;
}
double HadronWeights::down(int pdg) const {
  auto it = weightDown.find(pdg);
  return it == weightDown.end() ? 0. : it->second;
}
bool HadronWeights::has(int pdg) const {
  auto it = eventHas.find(pdg);
  return it != eventHas.end() && it->second;
}

struct HadronReweight::Renorm {
  std::map<std::tuple<int, int, int, int>, std::unique_ptr<TH2D>> h;  // (tar, nuc, pdg, shift)
};

// ---------------------------------------------------------------------------
HadronReweight::HadronReweight(const HadronReweightConfig& c) : fC(c) {
  if (fC.dataDir.empty()) fC.dataDir = ResolveHadronReweightDataDir("");
  // useDefaultParticles()
  fParticles[PDGPION] = true;
  fParticles[PDGPROTON] = false;
  fParticles[PDGNEUTRON] = false;
  // MAT's safeEval uses TF1::Integral, i.e. ROOT's default 1-D integrator. With MathMore available
  // that is GSL QAGS ("AdaptiveSingular", epsabs = epsrel = 1e-12). Since the integrands are
  // O(1e-25) the absolute tolerance dominates, so the integrator configuration matters at the
  // per-mille level: a different default integrator changes the results.
  if (ROOT::Math::IntegratorOneDimOptions::DefaultIntegrator() != "AdaptiveSingular")
    std::cerr << "HadronReweight: WARNING default 1D integrator is "
              << ROOT::Math::IntegratorOneDimOptions::DefaultIntegrator()
              << " (MAT used AdaptiveSingular/GSL); weights will differ slightly\n";
  buildFunctions();
}

HadronReweight::~HadronReweight() = default;

TGraph* HadronReweight::loadGraph(const std::string& file, const char* name) {
  const std::string path = fC.dataDir + "/hadronReweight/" + file;
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("HadronReweight: cannot open " + path);
  auto* g = dynamic_cast<TGraph*>(f->Get(name));
  if (!g) throw std::runtime_error("HadronReweight: TGraph " + std::string(name) + " missing in " + path);
  auto* out = static_cast<TGraph*>(g->Clone());
  out->Sort();  // as XSecFunctions::loadTGraph (sets kIsSortedX -> binary search in Eval)
  fGraphs.emplace_back(out);
  return out;
}

void HadronReweight::buildFunctions() {
  // XSecFunctions::loadTGraphs(hdXSec = true) (+ useReweightedNeutronCV's new_geant4_xsec graphs).
  TGraph* C_n_in = loadGraph("cross_section_neutron_carbon_hd.root", "inelXS");
  TGraph* C_n_el = loadGraph("cross_section_neutron_carbon_hd.root", "elXS");
  TGraph* Pb_n_in = loadGraph("cross_section_neutron_lead_hd.root", "inelXS");
  TGraph* Pb_n_el = loadGraph("cross_section_neutron_lead_hd.root", "elXS");
  TGraph* H_n_in = loadGraph("cross_section_neutron_hydrogen_hd.root", "inelXS");
  TGraph* H_n_el = loadGraph("cross_section_neutron_hydrogen_hd.root", "elXS");
  TGraph* C_p_in = loadGraph("cross_section_proton_carbon.root", "inelXS");
  TGraph* C_p_el = loadGraph("cross_section_proton_carbon.root", "elXS");
  TGraph* Pb_p_in = loadGraph("cross_section_proton_lead.root", "inelXS");
  TGraph* Pb_p_el = loadGraph("cross_section_proton_lead.root", "elXS");
  TGraph* C_pi_in = loadGraph("cross_section_pion_carbon.root", "inelXS");
  TGraph* C_pi_el = loadGraph("cross_section_pion_carbon.root", "elXS");
  TGraph* Fe_pi_in = loadGraph("cross_section_pion_iron.root", "inelXS");
  TGraph* Fe_pi_el = loadGraph("cross_section_pion_iron.root", "elXS");
  TGraph* Pb_pi_in = loadGraph("cross_section_pion_lead.root", "inelXS");
  TGraph* Pb_pi_el = loadGraph("cross_section_pion_lead.root", "elXS");
  TGraph *New_C_n_el = nullptr, *New_C_n_in = nullptr;
  if (fC.useNeutronCVReweight) {
    New_C_n_el = loadGraph("new_geant4_xsec.root", "elastic_xsec");
    New_C_n_in = loadGraph("new_geant4_xsec.root", "inelastic_xsec");
  }
  // MAT: the iron proton/neutron functions use the *lead* graphs, and cross_section_*_iron(_hd)
  // for protons/neutrons is therefore never used. Reproduced as is.

  auto mk = [&](const std::string& name, std::function<double(double)> fn) {
    auto f = std::make_unique<TF1>(("makedf_mhrw_" + name).c_str(), [fn](double* x, double*) { return fn(x[0]); },
                                   0.0, 5000.0, 0, 1, TF1::EAddToList::kNo);
    fFuncs[name] = std::move(f);
  };
  const double m0 = DEFAULT_MIN, m1 = DEFAULT_MAX, a0 = DEFAULT_MIN_ALT, a1 = DEFAULT_MAX_ALT;
  // Carbon
  mk("carbon_pion_inelastic", [=](double x) { return tgrapheval(C_pi_in, x, m0, m1, true); });
  mk("carbon_proton_inelastic", [=](double x) { return tgrapheval(C_p_in, x, m0, m1, true); });
  mk("carbon_neutron_inelastic", [=](double x) { return tgrapheval_withHD(C_n_in, x, m0, m1); });
  mk("carbon_pion_elastic", [=](double x) { return tgrapheval(C_pi_el, x, m0, m1, false); });
  mk("carbon_proton_elastic", [=](double x) { return tgrapheval(C_p_el, x, m0, m1, false); });
  mk("carbon_neutron_elastic", [=](double x) { return tgrapheval_withHD(C_n_el, x, m0, m1); });
  // Scintillator neutron (useHDXSec: carbon + hydrogen)
  mk("scint_neutron_inelastic",
     [=](double x) { return tgrapheval_withHD(C_n_in, x, m0, m1) + tgrapheval_withHD(H_n_in, x, m0, m1); });
  mk("scint_neutron_elastic",
     [=](double x) { return tgrapheval_withHD(C_n_el, x, m0, m1) + tgrapheval_withHD(H_n_el, x, m0, m1); });
  // Lead
  mk("lead_pion_inelastic", [=](double x) { return tgrapheval(Pb_pi_in, x, a0, a1, true); });
  mk("lead_proton_inelastic", [=](double x) { return tgrapheval(Pb_p_in, x, a0, a1, true); });
  mk("lead_neutron_inelastic", [=](double x) { return tgrapheval_withHD(Pb_n_in, x, a0, a1); });
  mk("lead_pion_elastic", [=](double x) { return tgrapheval(Pb_pi_el, x, a0, a1, false); });
  mk("lead_proton_elastic", [=](double x) { return tgrapheval(Pb_p_el, x, a0, a1, false); });
  mk("lead_neutron_elastic", [=](double x) { return tgrapheval_withHD(Pb_n_el, x, a0, a1); });
  // Iron (MAT: proton and neutron use the lead graphs with the iron range)
  mk("iron_pion_inelastic", [=](double x) { return tgrapheval(Fe_pi_in, x, m0, m1, true); });
  mk("iron_proton_inelastic", [=](double x) { return tgrapheval(Pb_p_in, x, m0, m1, true); });
  mk("iron_neutron_inelastic", [=](double x) { return tgrapheval_withHD(Pb_n_in, x, m0, m1); });
  mk("iron_pion_elastic", [=](double x) { return tgrapheval(Fe_pi_el, x, m0, m1, false); });
  mk("iron_proton_elastic", [=](double x) { return tgrapheval(Pb_p_el, x, m0, m1, false); });
  mk("iron_neutron_elastic", [=](double x) { return tgrapheval_withHD(Pb_n_el, x, m0, m1); });
  // New (reweighted CV) scintillator neutron
  if (fC.useNeutronCVReweight) {
    mk("new_carbon_neutron_elastic",
       [=](double x) { return tgrapheval_withHD(New_C_n_el, x, m0, m1) + tgrapheval_withHD(H_n_el, x, m0, m1); });
    mk("new_carbon_neutron_inelastic",
       [=](double x) { return tgrapheval_withHD(New_C_n_in, x, m0, m1) + tgrapheval_withHD(H_n_in, x, m0, m1); });
  }
}

const TF1* HadronReweight::getXSecFunc(int pdg, int nuke, XSecType t) const {
  if (pdg == -211) pdg = 211;
  const char* pname = pdg == PDGPION ? "pion" : pdg == PDGPROTON ? "proton" : pdg == PDGNEUTRON ? "neutron" : nullptr;
  if (!pname) return nullptr;  // kaons and anything else: no function
  std::string mat;
  if (nuke == 6 || nuke == -6) mat = (pdg == PDGNEUTRON && nuke == -6) ? "scint" : "carbon";
  else if (nuke == 26) mat = "iron";
  else if (nuke == 82) mat = "lead";
  else return nullptr;  // TODO in MAT: water, aluminium, ...
  return fFuncs.at(mat + "_" + pname + (t == kInelastic ? "_inelastic" : "_elastic")).get();
}

// MAT safeEval(TF1, Ti, Tf): same TF1::Integral call (ROOT default 1-D integrator, epsrel 1e-12).
double HadronReweight::safeEval(const TF1* fc, double Ti, double Tf) const {
  TF1* f = const_cast<TF1*>(fc);
  if (Ti == Tf) return f->Eval(Ti);
  return f->Integral(Ti, Tf) / (Tf - Ti);
}

double HadronReweight::getAverageCrossSection2(int pdg, int nuke, double Ti, double Tf, XSecType t,
                                               bool isCVCalc) const {
  if (fC.useNeutronCVReweight && pdg == PDGNEUTRON && nuke == -6 && !isCVCalc)
    return getAverageNewNeutronCrossSection(nuke, Ti, Tf, t);
  double integral = 0.0;
  if (const TF1* f = getXSecFunc(pdg, nuke, t)) integral = safeEval(f, Ti, Tf);
  return integral;
}

double HadronReweight::getAverageNewNeutronCrossSection(int nuke, double Ti, double Tf, XSecType t) const {
  if (nuke != -6) return 0.0;
  return safeEval(fFuncs.at(t == kElastic ? "new_carbon_neutron_elastic" : "new_carbon_neutron_inelastic").get(), Ti,
                  Tf);
}

bool HadronReweight::inFiducial(double X, double Y, double Z) const {
  if (Z < fC.minZ || Z > fC.maxZ) return false;
  return IsInHexagonMAT(X, Y, fC.apothem);
}

double HadronReweight::getDelta(int pdg, int target, double ke) const {
  if (fC.useNeutronCVReweight && pdg == PDGNEUTRON && target == -6) {
    if (ke > 200) return 0.02;
    if (ke > 100) return 0.03;
    if (ke > 25) return 0.15;
    if (ke > 10) return 0.20;
    return 0.25;
  }
  if (target == -6 || target == 6) {
    if (pdg == PDGPION || pdg == PDGPROTON) return 0.1;
    else if (pdg == PDGNEUTRON) return NeutronReweightAmount(ke);
  } else if (target == 26) {
    return 0.15;
  } else if (target == 82) {
    return 0.20;
  } else if (target == 13) {
    return 0;
  } else if (target == 8) {
    return 0.15;
  }
  // MAT: "This shouldn't happen!! :(" + throw 66
  throw std::runtime_error("HadronReweight::getDelta: unexpected pdg/target " + std::to_string(pdg) + "/" +
                           std::to_string(target));
}

double HadronReweight::getDeltaElastic(int pdg, int nuke, double ke) const {
  if (fC.useNeutronCVReweight && pdg == PDGNEUTRON && nuke == -6) {
    if (ke > 25) return 0.03;
    if (ke > 10) return 0.10;
    return 0.25;
  }
  if (pdg == PDGNEUTRON) {
    if (nuke == 82) {
      if (ke > 110) return 0.15;
      else if (ke > 30) return 0.25;
      else return 0.10;
    } else if (nuke == 26) {
      if (ke > 400) return 0.20;
      else if (ke > 160) return 0.25;
      else return 0.15;
    } else {
      if (ke > 80) return 0.20;
      else if (ke > 30) return 0.40;
      else if (ke > 10) return 0.25;
      else if (ke > 8) return 0.40;
      else if (ke > 5) return 0.70;
      else return 0.40;
    }
  }
  return 0.8;
}

bool HadronReweight::getWeightWithElastic(int nuke, int pdg, int intCode, double density, double Ei, double Ef,
                                          double delta_inelastic, double delta_elastic, double& up,
                                          double& down) const {
  if (density == 0) return false;
  const double mass = getMass(pdg);
  double Ti = Ei - mass;
  double Tf = Ef - mass;
  if (Ti < 0 || Tf < 0) {
    if (Ef == 0) return false;
    if (Ti <= 0 && Tf <= 0) return false;
    if (Ti > 0 && Tf < 0) Tf = 0.0;
  }
  const double rhoX = getCorrectedDensity(nuke, density);
  double sigma_elastic_geant = getAverageCrossSection2(pdg, nuke, Ti, Tf, kElastic);
  double sigma_inelastic_geant = getAverageCrossSection2(pdg, nuke, Ti, Tf, kInelastic);
  double sigma_total_geant = sigma_elastic_geant + sigma_inelastic_geant;
  if (sigma_total_geant <= 0 || sigma_elastic_geant <= 0 || sigma_inelastic_geant <= 0) return false;

  double sigma_total_data_up = sigma_elastic_geant * (1 + delta_elastic) + sigma_inelastic_geant * (1 + delta_inelastic);
  double sigma_total_data_down =
      sigma_elastic_geant * (1 - delta_elastic) + sigma_inelastic_geant * (1 - delta_inelastic);

  if (intCode != 0 && intCode != 2 && intCode != 3) {  // inelastic interacted (1, 4)
    double sp_geant = sigma_inelastic_geant;
    double sp_up = sp_geant * (1 + delta_inelastic);
    if (sp_up < 0) sp_up = 0;
    double sp_down = sp_geant * (1 - delta_inelastic);
    if (sp_down < 0) sp_down = 0;
    up = sigma_total_data_up <= 0 ? 0 : doWeightCalcInteracted(rhoX, sp_up, sp_geant, sigma_total_data_up, sigma_total_geant);
    down = sigma_total_data_down <= 0
               ? 0
               : doWeightCalcInteracted(rhoX, sp_down, sp_geant, sigma_total_data_down, sigma_total_geant);
  } else if (intCode == 3) {  // elastic
    double sp_geant = sigma_elastic_geant;
    double sp_up = sp_geant * (1 + delta_elastic);
    if (sp_up < 0) sp_up = 0;
    double sp_down = sp_geant * (1 - delta_elastic);
    if (sp_down < 0) sp_down = 0;
    up = sigma_total_data_up <= 0 ? 0 : doWeightCalcInteracted(rhoX, sp_up, sp_geant, sigma_total_data_up, sigma_total_geant);
    down = sigma_total_data_down <= 0
               ? 0
               : doWeightCalcInteracted(rhoX, sp_down, sp_geant, sigma_total_data_down, sigma_total_geant);
  } else {  // not interacted (0, 2)
    up = doWeightCalcNoninteracted(rhoX, sigma_total_data_up, sigma_total_geant);
    down = doWeightCalcNoninteracted(rhoX, sigma_total_data_down, sigma_total_geant);
  }
  if (std::isnan(up) || std::isnan(down)) throw std::runtime_error("HadronReweight::getWeightWithElastic: NaN");
  return true;
}

bool HadronReweight::getWeightToNewNeutronXSec(int nuke, int intCode, double density, double Ei, double Ef,
                                               double& wgt) const {
  if (density == 0) return false;
  if (nuke != -6) return false;
  const double mass = getMass(PDGNEUTRON);
  double Ti = Ei - mass;
  double Tf = Ef - mass;
  if (Ti < 0 || Tf < 0) {
    if (Ef == 0) return false;
    if (Ti <= 0 && Tf <= 0) return false;
    if (Ti > 0 && Tf < 0) Tf = 0.0;
  }
  const double rhoX = getCorrectedDensity(nuke, density);
  double sigma_elastic_geant = getAverageCrossSection2(PDGNEUTRON, nuke, Ti, Tf, kElastic, true);
  double sigma_inelastic_geant = getAverageCrossSection2(PDGNEUTRON, nuke, Ti, Tf, kInelastic, true);
  double sigma_total_geant = sigma_elastic_geant + sigma_inelastic_geant;
  if (sigma_total_geant <= 0 || sigma_elastic_geant <= 0 || sigma_inelastic_geant <= 0) return false;

  double sigma_elastic_data_new = getAverageNewNeutronCrossSection(nuke, Ti, Tf, kElastic);
  double sigma_inelastic_data_new = getAverageNewNeutronCrossSection(nuke, Ti, Tf, kInelastic);
  double sigma_total_data_new = sigma_elastic_data_new + sigma_inelastic_data_new;
  if (sigma_total_data_new <= 0 || sigma_elastic_data_new <= 0 || sigma_inelastic_data_new <= 0) return false;

  if (intCode != 0 && intCode != 2 && intCode != 3)
    wgt = doWeightCalcInteracted(rhoX, sigma_inelastic_data_new, sigma_inelastic_geant, sigma_total_data_new,
                                 sigma_total_geant);
  else if (intCode == 3)
    wgt = doWeightCalcInteracted(rhoX, sigma_elastic_data_new, sigma_elastic_geant, sigma_total_data_new,
                                 sigma_total_geant);
  else
    wgt = doWeightCalcNoninteracted(rhoX, sigma_total_data_new, sigma_total_geant);
  if (std::isnan(wgt)) throw std::runtime_error("HadronReweight::getWeightToNewNeutronXSec: NaN");
  return true;
}

bool HadronReweight::getWeight(int nuke, int pdg, int intCode, double density, double Ei, double Ef, double delta,
                               double& up, double& down) const {
  if (density == 0) return false;
  const double mass = getMass(pdg);
  double Ti = Ei - mass;
  double Tf = Ef - mass;
  if (Ti < 0 || Tf < 0) {
    if (Ef == 0) return false;
    if (Ti <= 0 && Tf <= 0) return false;
    if (Ti > 0 && Tf < 0) Tf = 0.0;
  }
  const double rhoX = getCorrectedDensity(nuke, density);
  double sigma_inelastic_geant = getAverageCrossSection2(pdg, nuke, Ti, Tf, kInelastic);
  if (sigma_inelastic_geant <= 0) return false;
  double si_up = sigma_inelastic_geant * (1 + delta);
  double si_down = sigma_inelastic_geant * (1 - delta);
  if (intCode != 0 && intCode != 2 && intCode != 3) {
    up = doWeightCalcInelasticOnly(rhoX, si_up, sigma_inelastic_geant);
    down = doWeightCalcInelasticOnly(rhoX, si_down, sigma_inelastic_geant);
  } else {
    up = doWeightCalcNoninteracted(rhoX, si_up, sigma_inelastic_geant);
    down = doWeightCalcNoninteracted(rhoX, si_down, sigma_inelastic_geant);
  }
  if (std::isnan(up) || std::isnan(down)) throw std::runtime_error("HadronReweight::getWeight: NaN");
  return true;
}

// ---------------------------------------------------------------------------
std::map<int, HadronReweight::Leading> HadronReweight::getLeadingFSParticles(const HadronEvent& b) const {
  std::map<int, std::pair<int, double>> leadingPartIdxE;
  int nFS = (int)b.fsPDG.size();
  if (nFS >= AL || nFS < 0) nFS = 0;  // MAT: "mc_nFSPart out of bounds"
  for (int iPart = 0; iPart < nFS; ++iPart) {
    const int pdg = b.fsPDG[iPart];
    const bool bPart = fParticles.count(pdg) > 0;
    const bool bAntiPart = fParticles.count(-pdg) > 0;
    if (!bPart && !bAntiPart) continue;
    if (bAntiPart && !fParticles.at(-pdg)) continue;
    const int part_pdg = bPart ? pdg : -pdg;
    auto it = leadingPartIdxE.find(part_pdg);
    if (it == leadingPartIdxE.end()) leadingPartIdxE[part_pdg] = {iPart, b.fsE[iPart]};
    else if (b.fsE[iPart] > it->second.second) it->second = {iPart, b.fsE[iPart]};
  }
  std::map<int, Leading> out;
  for (auto& l : leadingPartIdxE) {
    const int i = l.second.first;
    // TVector3::RotateX(numi_beam_angle_rad), Mag(), Theta()
    const double s = std::sin(numi_beam_angle_rad), c = std::cos(numi_beam_angle_rad);
    const double x = b.fsPx[i], yy = b.fsPy[i], zz = b.fsPz[i];
    const double y = c * yy - s * zz, z = s * yy + c * zz;
    const double perp = std::sqrt(x * x + y * y);
    const double theta = (x == 0.0 && y == 0.0 && z == 0.0) ? 0.0 : std::atan2(perp, z);
    out[l.first] = {std::sqrt(x * x + y * y + z * z), theta * 180 / 3.14159};
  }
  return out;
}

std::pair<int, int> HadronReweight::getTargetNuclei(const HadronEvent& b) const {
  int nuc = b.targetZ;
  const double x = b.vtx[0], y = b.vtx[1], z = b.vtx[2];
  int tar = -1;
  if (!IsInHexagonMAT(x, y, 850.)) return {200, 0};
  if (z < TP::TrackerFace) {
    if (InIron1VolMC(x, y, z)) { tar = 1; nuc = 26; }
    else if (InLead1VolMC(x, y, z)) { tar = 1; nuc = 82; }
    else if (InIron2VolMC(x, y, z)) { tar = 2; nuc = 26; }
    else if (InLead2VolMC(x, y, z)) { tar = 2; nuc = 82; }
    else if (InCarbon3VolMC(x, y, z)) { tar = 3; nuc = 6; }
    else if (InIron3VolMC(x, y, z)) { tar = 3; nuc = 26; }
    else if (InLead3VolMC(x, y, z)) { tar = 3; nuc = 82; }
    else if (InLead4VolMC(x, y, z)) { tar = 4; nuc = 82; }
    else if (InIron5VolMC(x, y, z)) { tar = 5; nuc = 26; }
    else if (InLead5VolMC(x, y, z)) { tar = 5; nuc = 82; }
    else if (InWaterTargetVolMC(x, y, z)) { tar = 7; nuc = 8; }
    if (tar > 0) return {tar, nuc};
    if (z < (TP::C1 - TP::T1Pb / 2)) { tar = 101; nuc = 6; }
    else if ((TP::C1 + TP::T1Pb / 2) < z && z < (TP::C2 - TP::T2Pb / 2)) { tar = 102; nuc = 6; }
    else if ((TP::C2 + TP::T2Pb / 2) < z && z < (TP::C3 - TP::T3Pb / 2)) { tar = 103; nuc = 6; }
    else if ((TP::C3 + TP::T3Pb / 2) < z && z < TP::WaterFace) { tar = 104; nuc = 6; }
    else if (TP::WaterBack < z && z < (TP::C4 - TP::T4Pb / 2)) { tar = 105; nuc = 6; }
    else if ((TP::C4 + TP::T4Pb / 2) < z && z < (TP::C5 - TP::T5Pb / 2)) { tar = 106; nuc = 6; }
    else if ((TP::C5 + TP::T5Pb / 2) < z && z < TP::TrackerFace) { tar = 107; nuc = 6; }
  } else if (TP::TrackerFace < z && z < TP::TrackerBack) {
    if (InTrackerMAT(x, y, z)) { tar = 6; nuc = 6; }
  } else {
    tar = 300; nuc = 0;
  }
  return {tar, nuc};
}

const TH2D* HadronReweight::renormHist(int tar, int nuc, int pdg, int shift) const {
  if (!fCurRenorm) throw std::runtime_error("HadronReweight: setPlaylist() not called");
  auto it = fCurRenorm->h.find({tar, nuc, pdg, shift});
  if (it == fCurRenorm->h.end())  // MAT would dereference a null TH2D* here
    throw std::runtime_error("HadronReweight: no kine renorm histogram for tar " + std::to_string(tar) + " nuc " +
                             std::to_string(nuc) + " pdg " + std::to_string(pdg));
  return it->second.get();
}

// ---------------------------------------------------------------------------
std::string HadronReweight::NormalizePlaylist(const std::string& playlist) {
  std::string t = playlist;
  if (playlist.find("minerva") != std::string::npos || playlist.find("Minerva") != std::string::npos)
    t = playlist.substr(7, playlist.size() - 7);
  std::transform(t.begin(), t.end(), t.begin(), ::toupper);
  return t;
}

std::string HadronReweight::renormFileName(const std::string& normPlaylist) const {
  std::string suffix = !normPlaylist.empty() ? "_Minerva" + normPlaylist : "";
  suffix += "_v10";
  suffix += "_HD";  // useHDXSec is always on once the kine file is loaded
  if (fC.useElastics) suffix += "_withElastics";
  if (fC.useNeutronCVReweight) suffix += "_RwgtNeutronCV";
  suffix += "_z" + std::to_string((int)fC.minZ) + "-" + std::to_string((int)fC.maxZ) + "_a" +
            std::to_string((int)fC.apothem);
  return fC.dataDir + "/mhrwKineRenorm/" + fC.projectName + suffix + ".root";
}

std::vector<std::string> HadronReweight::requiredFiles(const std::vector<std::string>& playlists) const {
  std::vector<std::string> out;
  for (const char* f : {"cross_section_neutron_carbon_hd.root", "cross_section_neutron_lead_hd.root",
                        "cross_section_neutron_hydrogen_hd.root", "cross_section_proton_carbon.root",
                        "cross_section_proton_lead.root", "cross_section_pion_carbon.root",
                        "cross_section_pion_iron.root", "cross_section_pion_lead.root"})
    out.push_back(fC.dataDir + "/hadronReweight/" + f);
  if (fC.useNeutronCVReweight) out.push_back(fC.dataDir + "/hadronReweight/new_geant4_xsec.root");
  for (auto& p : playlists) out.push_back(renormFileName(NormalizePlaylist(p)));
  return out;
}

void HadronReweight::setPlaylist(const std::string& playlist) {
  const std::string pl = NormalizePlaylist(playlist);
  if (pl == fPlaylist && fCurRenorm) return;
  fPlaylist = pl;
  auto& slot = fRenorm[pl];
  if (!slot) {
    // MnvHadronReweight::tryLoadingKineFile
    const std::string fn = renormFileName(pl);
    std::unique_ptr<TFile> f(TFile::Open(fn.c_str(), "READ"));
    if (!f || f->IsZombie())
      throw std::runtime_error("HadronReweight: cannot open kinematic renorm file " + fn +
                               " (MAT would regenerate it from the Truth tree; not supported)");
    auto r = std::make_unique<Renorm>();
    std::vector<int> pdgs;
    for (auto& p : fParticles) pdgs.push_back(p.first);
    if (fC.useNeutronCVReweight) pdgs.push_back(PDGNEUTRONCV);
    auto load = [&](const std::string& name, int tar, int nuc, int pdg, int shift) {
      auto* h = dynamic_cast<TH2D*>(f->Get(name.c_str()));
      if (!h) throw std::runtime_error("HadronReweight: cannot load " + name + " from " + fn);
      h = static_cast<TH2D*>(h->Clone());
      h->SetDirectory(nullptr);
      r->h[{tar, nuc, pdg, shift}].reset(h);
    };
    for (int iTar = 0; iTar < mhrw_nTars; ++iTar) {
      const int tar = mhrw_targets[iTar], nuc = mhrw_nuclei[iTar];
      const std::string tn = tar == 6 ? "tracker" : tar == 7 ? "water" : Form("tar%d_nuc%d", tar, nuc);
      for (int s = 0; s < 3; ++s)
        for (int pdg : pdgs) load(Form("kine_rw_%s_p_theta_%d_%s", tn.c_str(), pdg, mhrw_shift_names[s]), tar, nuc, pdg, s);
    }
    for (int iTar = 1; iTar <= mhrw_nNTR_Scint; ++iTar)
      for (int s = 0; s < 3; ++s)
        for (int pdg : pdgs) load(Form("kine_rw_NTR_Scint%d_p_theta_%d_%s", iTar, pdg, mhrw_shift_names[s]), 100 + iTar, 6, pdg, s);
    for (int s = 0; s < 3; ++s)
      for (int pdg : pdgs) {
        load(Form("kine_rw_OD_p_theta_%d_%s", pdg, mhrw_shift_names[s]), 200, 0, pdg, s);
        load(Form("kine_rw_EHCAL_p_theta_%d_%s", pdg, mhrw_shift_names[s]), 300, 0, pdg, s);
      }
    std::cout << "HadronReweight: loaded kinematic renorm " << fn << std::endl;
    slot = std::move(r);
  }
  fCurRenorm = slot.get();
}

std::string HadronReweight::PlaylistFromRun(int run) {
  // MINERvA-101 util/GetPlaylist.cpp, verbatim (duplicate keys: the first entry wins, as in std::map
  // initialisation, e.g. 126240 -> minervame6B and 127100 -> minervame6G).
  static const std::map<int, std::string> m = {
      {110000, "minervame1A"}, {114000, "minervame1A"}, {111000, "minervame1B"}, {115000, "minervame1B"},
      {111030, "minervame1C"}, {115030, "minervame1C"}, {111100, "minervame1D"}, {115100, "minervame1D"},
      {111325, "minervame1E"}, {115325, "minervame1E"}, {111490, "minervame1F"}, {115490, "minervame1F"},
      {110150, "minervame1G"}, {114150, "minervame1G"}, {113000, "minervame1L"}, {117000, "minervame1L"},
      {113020, "minervame1M"}, {117020, "minervame1M"}, {113270, "minervame1N"}, {117270, "minervame1N"},
      {113375, "minervame1O"}, {117375, "minervame1O"}, {112000, "minervame1P"}, {116000, "minervame1P"},
      {123000, "minervame5A"}, {127000, "minervame5A"}, {122000, "minervame6A"}, {126000, "minervame6A"},
      {122240, "minervame6B"}, {126240, "minervame6B"}, {122390, "minervame6C"}, {126390, "minervame6C"},
      {122565, "minervame6D"}, {126565, "minervame6D"}, {122730, "minervame6E"}, {126730, "minervame6E"},
      {122880, "minervame6F"}, {126880, "minervame6F"}, {123100, "minervame6G"}, {127100, "minervame6G"},
      {123250, "minervame6H"}, {126240, "minervame6H"}, {123425, "minervame6I"}, {127425, "minervame6I"},
      {123565, "minervame6J"}, {127100, "minervame6J"}, {131000, "minervame3A"}, {135000, "minervame3A"},
      {6038, "minervame1A"},   {7000, "minervame2"},    {9000, "minervame3"},    {16854, "minervame3A"},
      {16838, "minervame4A"},  {16842, "minervame4B"},  {16845, "minervame4C"},  {16848, "minervame4D"},
      {16849, "minervame4E"},  {10068, "minervame1B"},  {10129, "minervame1C"},  {10255, "minervame1D"},
      {16003, "minervame1E"},  {16523, "minervame1F"},  {18668, "minervame1G"},  {19167, "minervame1H"},
      {19168, "minervame1L"},  {19186, "minervame1M"},  {19500, "minervame5A"},  {22000, "minervame6A"},
      {22398, "minervame6B"},  {20328, "minervame1N"},  {20619, "minervame1O"},  {20800, "minervame1P"},
      {25001, "minervame6C"},  {25174, "minervame6D"},  {25283, "minervame6E"},  {25374, "minervame6F"},
      {25612, "minervame6G"},  {26000, "minervame6H"},  {26201, "minervame6I"},  {26277, "minervame6J"},
      {20300, "minervame7"}};
  auto found = m.upper_bound(run);
  if (found == m.begin()) return "NoSuchPlaylist";
  return std::prev(found)->second;
}

// ---------------------------------------------------------------------------
HadronWeights HadronReweight::getWeights(const HadronEvent& b) const {
  HadronWeights weights;
  // initializeWeights(weights, 1.0)
  for (auto& p : fParticles) {
    weights.weightUp[p.first] = 1.0;
    weights.weightDown[p.first] = 1.0;
    weights.eventHas[p.first] = false;
    if (p.second) {
      weights.weightUp[-p.first] = 1.0;
      weights.weightDown[-p.first] = 1.0;
      weights.eventHas[-p.first] = false;
    }
  }
  if (fC.useNeutronCVReweight) {
    weights.weightUp[PDGNEUTRONCV] = 1.0;
    weights.weightDown[PDGNEUTRONCV] = 1.0;
  }
  if (!b.hasBranches) return weights;
  // MAT throws if NPoints > AL (=4096); we have no fixed arrays and simply proceed.

  const size_t nPdg = b.pdg.size();
  int prevElasticTrackIndex = -1;
  for (int i = 0; i < b.nPoints; ++i) {
    const size_t ui = (size_t)i;
    if (ui >= b.density.size() || ui >= b.trackID.size() || ui >= b.intCodePerSegment.size() ||
        ui >= b.nuke.size() || ui >= b.finalE.size() || ui >= b.initialE.size() || ui >= b.posX.size() ||
        ui >= b.posY.size() || ui >= b.posZ.size()) {
      ++fBadIndex;  // MAT would read stale array contents here
      continue;
    }
    const double density = b.density[ui];
    if (density == 0) continue;
    const int index = b.trackID[ui];
    if (index < 0 || (size_t)index >= nPdg) continue;
    const int intCode = b.intCodePerSegment[ui];
    const int pdg = b.pdg[index];
    const int antiPDG = -pdg;
    bool treatAnti = false;
    auto itP = fParticles.find(pdg);
    if (itP == fParticles.end()) {
      auto itA = fParticles.find(antiPDG);
      if (itA == fParticles.end()) continue;
      else if (!itA->second) continue;
      else treatAnti = itA->second;
    } else treatAnti = itP->second;

    const int nuke = b.nuke[ui];
    const double Ef = b.finalE[ui];
    const double Ei = b.initialE[ui];
    double minE = 0.0;
    if (std::abs(pdg) == 211) minE = 141.8;
    if (std::abs(pdg) == 2212) minE = 943.6;
    if (std::abs(pdg) == 2112) minE = 951.5;
    if (Ef < minE && Ei < minE) continue;

    weights.eventHas[pdg] = true;
    if (treatAnti) weights.eventHas[antiPDG] = true;

    if (!inFiducial(b.posX[ui], b.posY[ui], b.posZ[ui])) continue;

    double up = 1.0, down = 1.0;
    // getDelta(b, iPrimTraj): |pdg| of the track, nuke of the segment, KE from the *final* energy
    const int apdg = std::abs(pdg);
    const double delta_inelastic = getDelta(apdg, nuke, Ef - getMass(apdg)) * 1.0 /*deltascale*/;
    bool success;
    if (pdg == PDGNEUTRON && fC.useElastics) {  // || elasticReweightForAllParticles (false)
      if (intCode == 3 && index == prevElasticTrackIndex) continue;
      const double ke = Ei - getMass(pdg);
      const double delta_elastic = getDeltaElastic(pdg, nuke, ke);
      success = getWeightWithElastic(nuke, pdg, intCode, density, Ei, Ef, delta_inelastic, delta_elastic, up, down);
    } else {
      success = getWeight(nuke, pdg, intCode, density, Ei, Ef, delta_inelastic, up, down);
    }
    if (success) {
      weights.weightUp[pdg] *= up;
      if (treatAnti) weights.weightUp[antiPDG] *= up;
      weights.weightDown[pdg] *= down;
      if (treatAnti) weights.weightDown[antiPDG] *= down;
      if (intCode == 3) prevElasticTrackIndex = index;
    }
  }

  // Kinematic renormalization (renormKine == true once the file is loaded; "renorm" is false)
  const auto tarnuc = getTargetNuclei(b);
  if (tarnuc.first > 0) {
    const auto fs = getLeadingFSParticles(b);
    for (auto& p : fs) {
      const int pdg = p.first;
      if (!weights.has(pdg)) continue;
      const TH2D* hcv = renormHist(tarnuc.first, tarnuc.second, pdg, CV);
      const int bin = hcv->FindFixBin(p.second.p, p.second.theta);
      if (hcv->GetBinContent(bin) > minKineEntries) {
        weights.weightUp[pdg] /= renormHist(tarnuc.first, tarnuc.second, pdg, UP)->GetBinContent(bin);
        weights.weightDown[pdg] /= renormHist(tarnuc.first, tarnuc.second, pdg, DOWN)->GetBinContent(bin);
      }
    }
  }
  return weights;
}

double HadronReweight::reweightNeutronCV(const HadronEvent& b) const {
  if (!fC.useNeutronCVReweight) return 1.0;
  double out = 1.0;
  if (!b.hasBranches) return 1.0;
  const size_t nPdg = b.pdg.size();
  int prevElasticTrackIndex = -1;
  for (int i = 0; i < b.nPoints; ++i) {
    const size_t ui = (size_t)i;
    if (ui >= b.density.size() || ui >= b.trackID.size() || ui >= b.intCodePerSegment.size() ||
        ui >= b.nuke.size() || ui >= b.finalE.size() || ui >= b.initialE.size() || ui >= b.posX.size() ||
        ui >= b.posY.size() || ui >= b.posZ.size()) {
      ++fBadIndex;
      continue;
    }
    const double density = b.density[ui];
    if (density == 0) continue;
    const int index = b.trackID[ui];
    const int intCode = b.intCodePerSegment[ui];
    if (index < 0 || (size_t)index >= nPdg) {  // MAT reads PDG[index] unchecked here
      ++fBadIndex;
      continue;
    }
    const int pdg = b.pdg[index];
    if (pdg != PDGNEUTRON) continue;
    const int nuke = b.nuke[ui];
    if (nuke != -6) continue;
    const double Ef = b.finalE[ui];
    const double Ei = b.initialE[ui];
    const double minE = 951.5;
    if (Ef < minE && Ei < minE) continue;
    if (!inFiducial(b.posX[ui], b.posY[ui], b.posZ[ui])) continue;
    double wgt = 1.0;
    if (intCode == 3 && index == prevElasticTrackIndex) continue;
    if (getWeightToNewNeutronXSec(nuke, intCode, density, Ei, Ef, wgt)) {
      out *= wgt;
      if (intCode == 3) prevElasticTrackIndex = index;
    }
  }
  const auto fs = getLeadingFSParticles(b);
  auto itn = fs.find(PDGNEUTRON);
  if (itn != fs.end()) {
    const auto tarnuc = getTargetNuclei(b);
    if (tarnuc.first > 0) {
      const TH2D* hcv = renormHist(tarnuc.first, tarnuc.second, PDGNEUTRONCV, CV);
      const int bin = hcv->FindFixBin(itn->second.p, itn->second.theta);
      if (hcv->GetBinContent(bin) > minKineEntries)
        out /= renormHist(tarnuc.first, tarnuc.second, PDGNEUTRONCV, UP)->GetBinContent(bin);
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
GeantWeights::GeantWeights(HadronReweight& hr, bool isTruth, bool legacyTruthCV, const std::string& playlist)
    : fHR(hr), fIsTruth(isTruth), fLegacyTruthCV(legacyTruthCV), fFixedPlaylist(playlist) {
  if (!fFixedPlaylist.empty()) fHR.setPlaylist(fFixedPlaylist);
}

void GeantWeights::request(BranchSet& bs) {
  // HadronBranchContainer: missing truth_hadronReweightNPaths -> "default mode" (all weights 1)
  fEv.hasBranches = bs.has("truth_hadronReweightNPaths");
  if (fEv.hasBranches) {
    fNPoints = bs.require("truth_hadronReweightNPoints");
    fIntCode = bs.require("truth_hadronReweightIntCodePerSegment");
    fNuke = bs.require("truth_hadronReweightNuke");
    fTrackID = bs.require("truth_hadronReweightTrackID");
    fPDG = bs.require("truth_hadronReweightPDG");
    fDens = bs.require("truth_hadronReweightColumnarDensity");
    fFinalE = bs.require("truth_hadronReweightFinalE");
    fInitE = bs.require("truth_hadronReweightInitialE");
    fX = bs.require("truth_hadronReweightPosX");
    fY = bs.require("truth_hadronReweightPosY");
    fZ = bs.require("truth_hadronReweightPosZ");
  } else {
    std::cerr << "GeantWeights: WARNING truth_hadronReweight branches missing; GEANT weights set to 1\n";
  }
  fTargetZ = bs.require("mc_targetZ");
  fVtx = bs.require("mc_vtx");
  fNFS = bs.require("mc_nFSPart");
  fFSPDG = bs.require("mc_FSPartPDG");
  fFSE = bs.require("mc_FSPartE");
  fFSPx = bs.require("mc_FSPartPx");
  fFSPy = bs.require("mc_FSPartPy");
  fFSPz = bs.require("mc_FSPartPz");
  if (fFixedPlaylist.empty()) fRun = bs.require("mc_run");
}

void GeantWeights::book(TTree* t) {
  BindOutput(t, "geant_p_dn", &p_dn, 'F');
  BindOutput(t, "geant_p_up", &p_up, 'F');
  BindOutput(t, "geant_n_dn", &n_dn, 'F');
  BindOutput(t, "geant_n_up", &n_up, 'F');
  BindOutput(t, "geant_pi_dn", &pi_dn, 'F');
  BindOutput(t, "geant_pi_up", &pi_up, 'F');
  BindOutput(t, "geant_neutron_cv", &neutron_cv, 'F');
}

namespace {
template <class T>
void FillVec(std::vector<T>& v, const InBranch* b) {
  const size_t n = b->size();
  v.resize(n);
  for (size_t i = 0; i < n; ++i) v[i] = (T)b->get(i);
}
}  // namespace

void GeantWeights::compute() {
  if (fRun) {
    const int run = (int)std::lround(fRun->get(0));
    if (run != fLastRun) {
      fHR.setPlaylist(HadronReweight::PlaylistFromRun(run));
      fLastRun = run;
    }
  }
  HadronEvent& e = fEv;
  if (e.hasBranches) {
    e.nPoints = (int)std::lround(fNPoints->get(0));
    FillVec(e.intCodePerSegment, fIntCode);
    FillVec(e.nuke, fNuke);
    FillVec(e.trackID, fTrackID);
    FillVec(e.pdg, fPDG);
    FillVec(e.density, fDens);
    FillVec(e.finalE, fFinalE);
    FillVec(e.initialE, fInitE);
    FillVec(e.posX, fX);
    FillVec(e.posY, fY);
    FillVec(e.posZ, fZ);
  }
  e.targetZ = (int)std::lround(fTargetZ->get(0));
  for (int k = 0; k < 3; ++k) e.vtx[k] = fVtx->get(k);
  const long nfs = std::lround(fNFS->get(0));
  FillVec(e.fsPDG, fFSPDG);
  FillVec(e.fsE, fFSE);
  FillVec(e.fsPx, fFSPx);
  FillVec(e.fsPy, fFSPy);
  FillVec(e.fsPz, fFSPz);
  if (nfs >= 0 && (size_t)nfs < e.fsPDG.size()) {  // mc_nFSPart is the count branch; keep consistent
    e.fsPDG.resize(nfs); e.fsE.resize(nfs); e.fsPx.resize(nfs); e.fsPy.resize(nfs); e.fsPz.resize(nfs);
  }

  const HadronWeights w = fHR.getWeights(e);
  p_dn = (float)w.ratio(2212, -1);
  p_up = (float)w.ratio(2212, +1);
  n_dn = (float)w.ratio(2112, -1);
  n_up = (float)w.ratio(2112, +1);
  pi_dn = (float)w.ratio(211, -1);
  pi_up = (float)w.ratio(211, +1);
  neutron_cv = (fIsTruth && fLegacyTruthCV) ? 1.f : (float)fHR.reweightNeutronCV(e);
}

}  // namespace makedf
