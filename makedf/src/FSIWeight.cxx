// Port of MAT-MINERvA weighters/weight_fsi.cxx (PlotUtils::fate_code / weight_fsi).
// See FSIWeight.h. Lines marked QUIRK reproduce behaviour of the original that looks
// unintended but is kept for exact agreement with MINERvA tunes.
#include "FSIWeight.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "TTree.h"

namespace makedf {
namespace fsi {

double GenieBEinMeV(int A) {
  // From UserPhysicsOptions.xml (GENIE v2 as used by MINERvA); hard coded for an exact match.
  if (A == 1) return 0;                    // hydrogen
  if (A == 6) return 17.0;                 // lithium
  if (A == 12) return 25.0;                // carbon
  if (A == 16) return 27.0;                // oxygen
  if (A == 24) return 32.0;                // magnesium
  if (A == 40) return 29.5;                // argon
  if (A == 48) return 30.0;                // Ti48
  if (A == 56) return 36.0;                // 56 iron
  if (A == 58) return 36.0;                // 58 nickel
  if (A >= 206 && A <= 208) return 44.0;   // 208 lead
  // MINERvA specials, picked off by hand
  if (A == 28) return 8.219751;  // silicon
  if (A == 27) return 8.115287;  // aluminum
  if (A == 14) return 7.185166;  // nitrogen
  if (A == 55) return 8.653063;  // iron55 or manganese55
  if (A == 35) return 8.347164;  // chlorine
  if (A == 4) return 5.0;        // rough
  return 8.0;                    // GENIE semi-empirical default
}

std::vector<IntranukeParticle> CalcFates(int mc_intType, int mc_targetA, int mc_er_nPart, const int* mc_er_ID,
                                         const int* mc_er_status, const int* mc_er_FD, const int* mc_er_LD,
                                         const double* mc_er_Px, const double* mc_er_Py, const double* mc_er_Pz,
                                         const double* mc_er_E) {
  // Fate codes are those of GENIE v2 hA (not hN / GENIE v3).
  std::vector<IntranukeParticle> out;
  for (int i = 0; i < mc_er_nPart; ++i) {
    if (mc_er_status[i] != 14) continue;
    const int id = mc_er_ID[i];
    // Only track fates for proton, neutron, pizero, piplus, piminus.
    if (!(id == 2212 || id == 2112 || id == 111 || std::abs(id) == 211)) continue;
    int tempfate = -1;
    const int fd = mc_er_FD[i];
    const int ld = mc_er_LD[i];

    // DEVIATION (safety only): the original indexes mc_er_ID[fd..ld] unconditionally, which is
    // undefined behaviour if the daughter range is invalid. Such particles get fate -1 here,
    // which (as in the original for fate <= 0) never contributes to either weight.
    const bool validDaughters = fd >= 0 && ld >= fd && ld < mc_er_nPart;
    if (validDaughters) {
      // QUIRK: original is `TMath::Abs(mc_er_ID[i] == 211)`, so pi- is NOT flagged as a pion here.
      const int isaPion = (id == 111 || id == 211) ? 1 : 0;

      int FSpions = 0;
      for (int jj = fd; jj <= ld; ++jj)
        if (std::abs(mc_er_ID[jj]) == 211 || mc_er_ID[jj] == 111) FSpions++;

      if (fd == ld) {
        // one daughter: 1 = no interaction, 3 = elastic, 5 = absorption into a 2000000xxx hadblob
        double offset = GenieBEinMeV(mc_targetA);
        double tolerance = 0.0001;
        if (std::fabs(offset - 8.0) < 0.1) tolerance = 1.2;  // unknown nuclei
        if (offset < 6.0) tolerance = 1.2;
        if (mc_intType != 1) offset = 0.0;
        if (mc_er_ID[fd] > 2000000000) {
          tempfate = 5;
        } else if (std::fabs(mc_er_E[i] - mc_er_E[fd] - offset) < tolerance) {
          tempfate = 1;  // no scattering
        } else if (mc_er_ID[fd] != mc_er_ID[i]) {
          tempfate = 2;  // charge exchange
        } else {
          tempfate = 3;  // elastic
        }
      } else if (FSpions > 0 && !isaPion) {
        tempfate = 8;  // nucleon (or, by the QUIRK above, pi-) goes to pions
      } else if (isaPion && FSpions == 0) {
        tempfate = 5;
        if (ld - fd == 1) {  // two-body absorption: nn / pn / pp
          int p = 0, n = 0;
          if (mc_er_ID[fd] == 2212) p++;
          if (mc_er_ID[fd] == 2112) n++;
          if (mc_er_ID[ld] == 2212) p++;
          if (mc_er_ID[ld] == 2112) n++;
          if (p == 1 && n == 1) tempfate = 51;
          if (p == 2) tempfate = 52;
          if (n == 2) tempfate = 50;
        }
      } else if (isaPion && FSpions >= 2) {
        tempfate = 8;
      } else if (isaPion && FSpions == 1) {
        tempfate = 4;
        for (int jj = fd; jj <= ld; ++jj)
          if ((mc_er_ID[jj] == 111 || std::abs(mc_er_ID[jj]) == 211) && mc_er_ID[jj] != id) tempfate = 2;
      } else {
        tempfate = 4;  // nucleon: can't tell 2 (CEX) from 4
      }
    }

    const double px = mc_er_Px[i], py = mc_er_Py[i], pz = mc_er_Pz[i], e = mc_er_E[i];
    const double p = std::sqrt(px * px + py * py + pz * pz);
    const double m = std::sqrt(e * e - p * p);
    const double k = e - m;
    IntranukeParticle part;
    part.pdg = id;
    part.fate = tempfate;
    part.ke = k;
    out.push_back(part);
  }
  // same comparator and algorithm as the original (std::sort, decreasing KE)
  std::sort(out.begin(), out.end(), [](const IntranukeParticle& a, const IntranukeParticle& b) { return a.ke > b.ke; });
  return out;
}

//-------------------------------------------------------------------------------------------
void WeightFSI::Reset() {
  fElasticWeight1 = 1.0;
  fElasticWeight2 = 1.0;
  fAbsorptionWeight = 1.0;
  count_elastic_nucleon = 0;
  count_elastic_pion = 0;
  elastic_weight_nucleon = -1.0;
  elastic_weight_pion = -1.0;
}

double WeightFSI::GetElasticWeight(int config) const {
  if (config == 1) return fElasticWeight1;
  if (config == 2) return fElasticWeight2;
  std::cerr << "Invalid config: " << config << std::endl;
  return 1.0;
}

void WeightFSI::calcWeights(int mc_charm, int mc_intType, int mc_targetA, int mc_targetZ, int mc_resID,
                            int mc_er_nPart, const int* mc_er_ID, const int* mc_er_status, const int* mc_er_FD,
                            const int* mc_er_LD, const double* mc_er_Px, const double* mc_er_Py,
                            const double* mc_er_Pz, const double* mc_er_E) {
  Reset();
  std::vector<IntranukeParticle> particles = CalcFates(mc_intType, mc_targetA, mc_er_nPart, mc_er_ID, mc_er_status,
                                                       mc_er_FD, mc_er_LD, mc_er_Px, mc_er_Py, mc_er_Pz, mc_er_E);

  if ((!mc_charm) && mc_targetA >= 4) {
    const bool is_delta = (mc_intType == 2) && (mc_resID == 0);
    for (auto& p : particles) {
      const double ke = p.ke / 1.e3;  // GeV
      double elastic_weight1 = -1.0;
      double elastic_weight2 = -1.0;
      GetQEElasticWeight(elastic_weight1, elastic_weight2, p.fate, ke, mc_targetA, p.pdg);
      // "Rik's magic"
      if (is_delta && ke > 0.1 && ke < 0.3 && elastic_weight1 > 1.05)
        elastic_weight1 = 1.0 + (elastic_weight1 - 1.0) * 0.90;
      p.weight1 = elastic_weight1;
      p.weight2 = elastic_weight2;
    }

    double pionElasticWeight[4] = {-1.0, -1.0, -1.0, -1.0};
    double nucleonElasticWeight[4] = {-1.0, -1.0, -1.0, -1.0};
    int index_pion = 0, index_nucleon = 0;
    for (const auto& p : particles) {
      const int pdg = p.pdg;
      const int fate = p.fate;
      if (fate != 1 && fate != 3) continue;
      if (fUseTrackingThreshold) {
        if (pdg == 2212 && p.ke < 125.0) continue;  // proton tracking threshold 125 MeV
      }
      // QUIRK: fate-1/3 particles below the binding-energy threshold keep weight1 = -1 but
      // still occupy one of the 4 slots (they are then not counted below).
      if ((index_nucleon < 4) && (pdg == 2212 || pdg == 2112)) {
        nucleonElasticWeight[index_nucleon] = p.weight1;
        ++index_nucleon;
      }
      if ((index_pion < 4) && (pdg == 211 || pdg == 111 || pdg == -211)) {
        pionElasticWeight[index_pion] = p.weight1;
        ++index_pion;
      }
    }

    double weightElastic = 1.0;
    int npi = 0, nnuc = 0, ntot = 0;
    // kludge that approximately preserves probability
    const double factor = 0.70;
    for (int i = 0; i < 4; i++) {
      if (pionElasticWeight[i] >= -0.001) npi++;
      if (nucleonElasticWeight[i] >= -0.001) nnuc++;
    }
    count_elastic_nucleon = nnuc;
    count_elastic_pion = npi;
    if (nnuc > 0) elastic_weight_nucleon = nucleonElasticWeight[0];
    if (npi > 0) elastic_weight_pion = pionElasticWeight[0];
    ntot = npi + nnuc;

    if (ntot == 0) {
      // no weight
    } else if (ntot == 1) {
      // Uses slot [0]; since particles are KE-sorted, a below-threshold (-1) slot 0 implies all
      // later slots of that species are -1 too, so slot 0 holds the single counted weight.
      if (nucleonElasticWeight[0] >= -0.001) weightElastic = nucleonElasticWeight[0];
      else weightElastic = pionElasticWeight[0];
    } else if (ntot == 2) {
      const bool is_res = mc_intType == 2;
      for (int i = 0; i < 4; i++) {
        if (pionElasticWeight[i] >= -0.001) weightElastic *= pionElasticWeight[i];
        if (is_res && std::fabs(pionElasticWeight[i] - 1.0) <= 0.0001) weightElastic *= 1.1;
        if (std::fabs(nucleonElasticWeight[i]) < 0.001) weightElastic *= 0.0;
        else if (nucleonElasticWeight[i] >= 0.999) weightElastic *= 1.0 + (nucleonElasticWeight[i] - 1.0) * factor;
      }
    } else if (ntot == 3) {
      for (int i = 0; i < 4; i++) {
        // QUIRK: float, so the factor is 1.0299999713897705 in double arithmetic
        const float factorenhance = 1.03;
        if (std::fabs(pionElasticWeight[i]) < 0.001) weightElastic *= 0.0;
        else if (pionElasticWeight[i] >= 1.005) weightElastic *= pionElasticWeight[i];
        else if (pionElasticWeight[i] >= 0.999) weightElastic *= factorenhance;
        if (std::fabs(nucleonElasticWeight[i]) < 0.001) weightElastic *= 0.0;
        else if (nucleonElasticWeight[i] >= 1.005) weightElastic *= 1.0 + (nucleonElasticWeight[i] - 1.0) * factor;
        else if (nucleonElasticWeight[i] >= 0.999) weightElastic *= factorenhance;
      }
    } else {
      // more than three weights: leave at 1 (original: "skip this event")
    }
    fElasticWeight1 = weightElastic;
    // NB: fElasticWeight2 is never set by the original (stays 1.0).
  }

  // ---- pion absorption weight (applied for any target / charm) ----
  double weightPionAbsorptionBug = 1.0;
  for (const auto& p : particles) {
    if (p.ke <= 0.0) break;
    if (p.pdg != 211) continue;
    const int piplus14fate = p.fate;
    double thisZoverA = (double)mc_targetZ / (double)mc_targetA;
    double wrongpn = 0.88 * (1 - thisZoverA) * (thisZoverA);
    double wrongpp = 0.14 * (thisZoverA) * (thisZoverA);
    double wrongnn = 0.14 * (1 - thisZoverA) * (1 - thisZoverA);
    double wrongdenom = wrongpn + wrongpp + wrongnn;
    wrongpp = (1.0 - wrongpn / wrongdenom) * (wrongpn + wrongpp) / wrongdenom;
    wrongpn = wrongpn / wrongdenom;
    wrongnn = (1.0 - wrongpn - wrongpp);
    (void)wrongnn;
    double rightpn = 0.083 * (1 - thisZoverA) * (1 - thisZoverA);
    double rightpp = 2.0 * (1 - thisZoverA) * (thisZoverA);
    double rightnn = 0.0;
    double rightdenom = rightpn + rightpp + rightnn;
    rightpn = rightpn / rightdenom;
    rightpp = rightpp / rightdenom;
    if (piplus14fate == 51) weightPionAbsorptionBug *= rightpn / wrongpn;
    else if (piplus14fate == 50) weightPionAbsorptionBug *= 0.0;
    else if (piplus14fate == 52) weightPionAbsorptionBug *= rightpp / wrongpp;
  }
  for (const auto& p : particles) {
    if (p.ke <= 0.0) break;
    if (p.pdg != 111) continue;
    const int pizero14fate = p.fate;
    double thisZoverA = (double)mc_targetZ / (double)mc_targetA;
    double rightpn = 0.88 * (1 - thisZoverA) * (thisZoverA);
    double rightpp = 0.14 * (thisZoverA) * (thisZoverA);
    double rightnn = 0.14 * (1 - thisZoverA) * (1 - thisZoverA);
    double rightdenom = rightpn + rightpp + rightnn;
    double wrongpp = (1.0 - rightpn / rightdenom) * (rightpn + rightpp) / rightdenom;
    double wrongpn = rightpn / rightdenom;
    double wrongnn = (1.0 - wrongpn - wrongpp);
    rightnn = rightnn / rightdenom;
    rightpp = rightpp / rightdenom;
    if (pizero14fate == 51) weightPionAbsorptionBug *= 1.0;
    else if (pizero14fate == 50) weightPionAbsorptionBug *= rightnn / wrongnn;
    else if (pizero14fate == 52) weightPionAbsorptionBug *= rightpp / wrongpp;
  }
  // QUIRK: pi- absorption is never reweighted (and, by the isaPion quirk, pi- never gets 5x fates).
  fAbsorptionWeight = weightPionAbsorptionBug;
}

int WeightFSI::GetQEElasticWeight(double& weightElastic1, double& weightElastic2, int fatecodePN, double inputKEgev,
                                  int inputA, int pdg) const {
  const double GeV = 1.e3;
  const double smallKElargeWeightBuffer = 0.005;
  const double lowthreshold = GenieBEinMeV(inputA) / GeV + smallKElargeWeightBuffer;
  if (inputA <= 4) return 1;
  if (inputA >= 220) return 1;
  if (fatecodePN <= 0) return 1;
  double tempweightElastic1 = 1.0;
  double tempweightElastic2 = 1.0;
  int errorcode = GetQEElasticFSIWeightFromParam(tempweightElastic1, tempweightElastic2, inputKEgev, inputA, pdg);
  if (inputKEgev > lowthreshold) {
    if (fatecodePN == 1) { weightElastic1 = tempweightElastic1; }
    else if (fatecodePN != 1 && fatecodePN != 3) { weightElastic2 = tempweightElastic2; }
    else if (fatecodePN == 3) { weightElastic1 = 0.0; weightElastic2 = 0.0; }
  }
  return errorcode;
}

int WeightFSI::GetQEElasticFSIWeightFromParam(double& weight1, double& weight2, double inputKEgev, int inputA,
                                              int pdg) const {
  // Polynomial parameterisations, copied verbatim from weight_fsi.cxx.
  static const double protonE1polyHe4[10] = { 39.32856672, -3930.858433, 198879.6178, -5992910.025, 114791528.5, -1431038143, 1.155973558e+10, -5.830506854e+10, 1.667805605e+11, -2.064812952e+11};
  static const double protonE2polyHe4[10] = { 2.448600795, -21.2815073, 175.9270834, -837.1818884, 2452.148013, -4549.790084, 5362.000843, -3893.423638, 1589.436743, -279.1563134};
  static const double protonE3polyHe4[10] = { 1.156044571, -0.03764619617, 0.005627980527, 0., 0., 0., 0., 0., 0., 0.};
  
  static const double protonE1polyC12[10] = { 17.82282668, 1551.604933, -170931.3066, 6480983.701, -134101728.6, 1696097963, -1.350258893e+10, 6.620396447e+10, -1.829733656e+11, 2.184157212e+11};
  static const double protonE2polyC12[10] = { 4.174891296, -45.54402543, 345.7329055, -1484.999656, 3875.448189, -6356.522435, 6588.903725, -4194.602797, 1499.087099, -230.4838552};
  static const double protonE3polyC12[10] = { 1.231471752, -0.05270919298, 0.008157195446, 0., 0., 0., 0., 0., 0., 0.};
  
  static const double protonE1polyO16[10] = { -11.42825559, 4987.758032, -340193.7326, 11199009.42, -216572955.5, 2635669472, -2.048933232e+10, 9.894487269e+10, -2.708119885e+11, 3.213193282e+11};
  static const double protonE2polyO16[10] = { 4.858781299, -55.02831481, 418.1617951, -1798.993092, 4700.672222, -7716.976536, 8001.460326, -5089.561265, 1814.087529, -277.4558181};
  static const double protonE3polyO16[10] = { 1.2676478, -0.0583832634, 0.008531849866, 0., 0., 0., 0., 0., 0., 0.};

  static const double protonE1polySi28[10] = { 323.0786729, -3756.15957, -732485.8957, 37802444.09, -882957457.5, 1.192714772e+10, -9.873018114e+10, 4.952931632e+11, -1.385578606e+12, 1.661480791e+12};
  static const double protonE2polySi28[10] = { 24.61820776, -403.5263055, 3189.486319, -14362.55838, 40075.08404, -71634.16657, 82127.95804, -58397.4893, 23440.19114, -4058.274728};
  static const double protonE3polySi28[10] = { 1.40181183, -0.09856004905, 0.01487881423, 0., 0., 0., 0., 0., 0., 0.};

  
  static const double protonE1polyAr40[10] = { 1193.624967, -55946.07016, 1113419.834, -11908293.37, 71837495.09, -231106211.9, 309219073.5, 0, 0, 0};
  static const double protonE2polyAr40[10] = { 31.09732866, -517.1986901, 4089.582571, -18404.71037, 51320.2566, -91701.47518, 105134.2984, -74779.77238, 30032.94532, -5203.640836};
  static const double protonE3polyAr40[10] = { 1.474620561, -0.1118058107, 0.01667045434, 0., 0., 0., 0., 0., 0., 0.};

  static const double protonE1polyFe56[10] = { 2040.917669, -106698.8941, 2353861.651, -27766888.52, 183882774.5, -646544109.9, 941601815.8, 0, 0, 0};
  static const double protonE2polyFe56[10] = { 31.90738462, -525.4002497, 4131.334337, -18520.37588, 51485.77356, -91755.43134, 104939.7732, -74465.92015, 29837.60618, -5157.844509};
  static const double protonE3polyFe56[10] = { 1.526848651, -0.1274368016, 0.01986318525, 0., 0., 0., 0., 0., 0., 0.};

  static const double protonE1polyPb208[10] = { 5049.811026, -264931.1183, 5844125.51, -68809380.99, 454319447.8, -1591535516, 2308315357, 0, 0, 0};
  static const double protonE2polyPb208[10] = { 61.61613899, -1022.962821, 7948.199277, -35212.29804, 96827.23696, -170874.4408, 193710.7751, -136368.3507, 54248.0508, -9316.087101};
  static const double protonE3polyPb208[10] = { 1.901054456, -0.1953860148, 0.02861291685, 0., 0., 0., 0., 0., 0., 0.};

  ///

  static const double neutronE1polyHe4[10] = { 76.21450444, -7284.080281, 350409.6144, -10144233.29, 188295757, -2289715752, 1.813028138e+10, -8.996810885e+10, 2.539078252e+11, -3.108221822e+11};
  static const double neutronE2polyHe4[10] = { 3.819189505, -41.52923432, 328.7530687, -1528.585143, 4443.389389, -8283.145704, 9883.531893, -7291.837115, 3026.956436, -540.240589};
  static const double neutronE3polyHe4[10] = { 1.1301232, -0.02361999224, 0.003182730455, 0., 0., 0., 0., 0., 0., 0.};
  
  
  static const double neutronE1polyC12[10] = { 91.8661909, -3841.369929, 17803.80266, 2466426.335, -77777993.43, 1160008348, -1.006039598e+10, 5.189316382e+10, -1.481598665e+11, 1.807900697e+11};
  static const double neutronE2polyC12[10] = { 6.87059066, -87.0000531, 658.0906236, -2865.934219, 7713.872241, -13257.88765, 14592.16131, -9963.08823, 3845.829084, -641.7424001};
  static const double neutronE3polyC12[10] = { 1.191527119, -0.02739756607, 0.003201116244, 0., 0., 0., 0., 0., 0., 0.};

  // this looks radically different.  Huh?
  static const double neutronE1polyO16[10] = { -2.15186349, 6545.892859, -480405.8834, 16120949.72, -313111105.5, 3803439779, -2.941740429e+10, 1.410667677e+11, -3.829205886e+11, 4.502057086e+11};
  static const double neutronE2polyO16[10] = { 7.669365414, -96.49518275, 722.5640045, -3118.543517, 8308.536397, -14117.09555, 15347.33815, -10345.77359, 3942.650035, -649.6469222};
  static const double neutronE3polyO16[10] = { 1.2325431, -0.04131824065, 0.005846054663, 0., 0., 0., 0., 0., 0., 0.};

  static const double neutronE1polySi28[10] = { 1942.270031, -145015.1637, 4619879.234, -74613366.29, 537963557.4, 938586048, -4.898302179e+10, 3.844432317e+11, -1.353278306e+12, 1.872616614e+12};
  static const double neutronE2polySi28[10] = { 45.15947021, -752.4355742, 5833.867552, -25789.37417, 70848.97048, -125102.6277, 142078.2779, -100275.8391, 40006.37836, -6891.096095};
  static const double neutronE3polySi28[10] = { 1.342152187, -0.06350811765, 0.008383670794, 0., 0., 0., 0., 0., 0., 0.};
  
  static const double neutronE1polyAr40[10] = { 1675.010584, -72524.21889, 1329865.841, -13064326.4, 72102012.51, -211141714.4, 255512775, 0., 0., 0.};
  static const double neutronE2polyAr40[10] = { 53.50240758, -895.3341982, 6935.402331, -30619.68803, 84012.55091, -148171.6227, 168093.1796, -118513.9921, 47235.85732, -8128.517134};
  static const double neutronE3polyAr40[10] = { 1.419983761, -0.08493487878, 0.01244785896, 0., 0., 0., 0., 0., 0., 0.};
  
  static const double neutronE1polyFe56[10] = { 3241.107116, -165590.9076, 3582114.609, -41549932.77, 271131806.6, -940872093, 1354028645, 0., 0., 0.};
  static const double neutronE2polyFe56[10] = { 54.64574278, -907.6702308, 7006.764962, -30871.98759, 84591.8376, -149056.734, 168992.8464, -119102.8865, 47461.72133, -8167.162769};
  static const double neutronE3polyFe56[10] = { 1.443437845, -0.07152361183, 0.00858548564, 0., 0., 0., 0., 0., 0., 0.};

  static const double neutronE1polyPb208[10] = { 6635.553396, -337195.0226, 7238352.564, -83226870.84, 538080687.7, -1849719333, 2637236834, 0., 0., 0.};
  static const double neutronE2polyPb208[10] = { 96.25604035, -1597.334355, 12200.22151, -53178.48721, 144204.542, -251602.1285, 282621.0265, -197471.025, 78062.22019, -13333.6168};
  static const double neutronE3polyPb208[10] = { 1.783187249, -0.121317479, 0.01488332601, 0., 0., 0., 0., 0., 0., 0.};

  
  // 
  // for weights to otherFSI.   These are same neutrino and anti-neutrino
  //

  static const double lowEpolyHe4w2[10] = { 0.8525499509, 106.5613564, -1349.864879, -129938.4005, 5706909.418, -106576730.4, 1102103401, -6559583996, 2.108872155e+10, -2.843225144e+10};
  static const double midEpolyHe4w2[10] ={ 1.401592311, 0.5340768238, 35.9410343, -351.4320019, 1465.967144, -3424.247666, 4784.59239, -3969.621213, 1803.437515, -345.4102984};

  static const double highEpolyHe4w2[10] = { 1.158515008, 0.0001009351013, -0.0001265995844, 0., 0., 0., 0., 0., 0., 0.};

  
  static const double lowEpolyC12w2[10] = { -1.09785639, 353.8971638, -14438.15969, 249803.2065, -961138.7354, -33007787.3, 592230359.6, -4425023560, 1.622582081e+10, -2.391266333e+10};
  static const double midEpolyC12w2[10] = { 1.405735843, 0.5803084824, 33.33315803, -323.6899774, 1330.142548, -3060.151485, 4217.656957, -3458.103359, 1555.362503, -295.3732419};

  static const double highEpolyC12w2[10] = { 1.16047678, -0.002525807226, 0.0006168334479, 0., 0., 0., 0., 0., 0., 0.};

  
  static const double lowEpolyO16w2[10] = { -1.600476996, 435.9526759, -19971.39691, 452148.723, -5407872.958, 28338610.09, 57010988.68, -1558249122, 7621399671, -1.283154539e+10};
  static const double midEpolyO16w2[10] = { 1.448614827, -0.4064702399, 42.84475735, -374.3659837, 1495.641126, -3402.083242, 4662.460835, -3809.164928, 1708.297187, -323.5100567};
  static const double highEpolyO16w2[10] = { 1.154381145, 0.004104892683, -0.0008410639001, 0., 0., 0., 0., 0., 0., 0.};

  static const double lowEpolySi28w2[10] = { -11.14679259, 1619.43455, -82842.44543, 2334771.427, -40502433.64, 451579395.5, -3251861347, 1.464296207e+10, -3.753948905e+10, 4.186264037e+10};
  static const double midEpolySi28w2[10] = { 1.402604068, 0.6760161318, 32.45166596, -320.1710152, 1323.602905, -3055.847889, 4221.183186, -3465.956677, 1560.301404, -296.4792756};
  static const double highEpolySi28w2[10] = { 1.15626498, 0.001690036606, -0.0002929508535, 0., 0., 0., 0., 0., 0., 0.};
  
  static const double lowEpolyAr40w2[10] = { -15.91880593, 2230.471033, -116235.7268, 3358158.198, -59926068.7, 688940043.9, -5124126962, 2.385706678e+10, -6.326926353e+10, 7.298762608e+10};
  static const double midEpolyAr40w2[10] =  { 1.370127894, 1.510992183, 23.88831308, -274.2049331, 1180.476078, -2784.735678, 3904.85482, -3244.201353, 1474.709714, -282.5053889};
  static const double highEpolyAr40w2[10] = { 1.157107574, 0.001164884679, -0.0003047637304, 0., 0., 0., 0., 0., 0., 0.};

  static const double lowEpolyFe56w2[10] = { -17.55855486, 2416.340033, -125628.1801, 3636120.214, -65207995.22, 755278024.5, -5671137415, 2.669951459e+10, -7.169405909e+10, 8.382943525e+10};
  static const double midEpolyFe56w2[10] = { 1.364364871, 1.617991401, 22.78343352, -266.8142579, 1149.241143, -2702.282161, 3770.631447, -3114.17639, 1406.301429, -267.5046458};
  static const double highEpolyFe56w2[10] = { 1.158295471, -0.0006771545028, 0.0002835267823, 0., 0., 0., 0., 0., 0., 0.};

  static const double lowEpolyPb208w2[10] = { -22.77635755, 3056.796849, -158957.7289, 4605296.37, -82631008.84, 956843795.2, -7177117819, 3.372961824e+10, -9.035554312e+10, 1.053493431e+11};
  static const double midEpolyPb208w2[10] = { 1.380837506, 1.340129355, 24.57882428, -272.4572318, 1158.247826, -2708.591884, 3770.535156, -3112.007007, 1405.812079, -267.6767955};
  static const double highEpolyPb208w2[10] = { 1.156387959, 0.001577782994, -0.0003703094127, 0., 0., 0., 0., 0., 0., 0.};



      // for pions elastic
  // For weights to pions, there are four ranges not three.
  static const double piplusE1polyC12[10] = {9.94936346e-01, 7.74247776e+00, -4.23721592e+02, 1.32366985e+04, -2.33889397e+05, 2.79519162e+06, -2.21397423e+07, 1.18115508e+08, -3.83775844e+08, 5.45834201e+08 };
  static const double piplusE2polyC12[10] = {-5.39724541e+01, 1.55691418e+03, -1.67801045e+04, 9.93671606e+04, -3.67977823e+05, 8.96122625e+05, -1.44363781e+06, 1.48618119e+06, -8.86776372e+05, 2.33312818e+05 };
  static const double piplusE3polyC12[10] = {-3.78237012e+02, 2.79989352e+03, -8.26927857e+03, 1.20837207e+04, -8.14884538e+03, 6.46957629e+02, 1.74115927e+03, 4.69547890e+01, -7.63025197e+02, 2.41820393e+02 };
  static const double piplusE4polyC12[10] = {1.22087378e+00, -1.00805455e-01, 1.84681810e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const double pizeroE1polyC12[10] = {9.72689495e-01, 1.49689756e+01, -9.90177911e+02, 3.55450768e+04, -7.21870395e+05, 9.12147217e+06, -7.14915776e+07, 3.45439500e+08, -9.53865861e+08, 1.14643003e+09 };
  static const double pizeroE2polyC12[10] = {-1.46650952e+02, 3.93088339e+03, -4.28398943e+04, 2.61209013e+05, -9.96804542e+05, 2.48538114e+06, -4.06181960e+06, 4.20231531e+06, -2.49950732e+06, 6.51500752e+05 };
  static const double pizeroE3polyC12[10] = {6.63546229e+02, -6.86281963e+03, 3.05525328e+04, -7.60833778e+04, 1.15413063e+05, -1.08313669e+05, 6.01468700e+04, -1.67258496e+04, 7.73760912e+02, 4.37044419e+02 };
  static const double pizeroE4polyC12[10] = {1.19889836e+00, -8.92934830e-02, 1.65346281e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const double piminusE1polyC12[10] = {1.00219141e+00, 6.20103564e+00, -2.93890871e+02, 7.52521760e+03, -9.05524938e+04, 6.47415899e+05, -2.68078401e+06, 1.39463349e+07, -8.15626192e+07, 1.80637117e+08 };
  static const double piminusE2polyC12[10] = {-1.94103153e+02, 5.09910032e+03, -5.56735222e+04, 3.42809070e+05, -1.32569322e+06, 3.35387690e+06, -5.56192442e+06, 5.83553454e+06, -3.51657942e+06, 9.27643678e+05 };
  static const double piminusE3polyC12[10] = {5.41948092e+01, -1.01098776e+03, 5.84798413e+03, -1.59287600e+04, 2.23801446e+04, -1.36642006e+04, -3.07297926e+03, 9.92854217e+03, -5.63490320e+03, 1.10208433e+03 };
  static const double piminusE4polyC12[10] = {1.21189831e+00, -9.20444945e-02, 1.67066575e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  if (pdg == 2212 || pdg == 2112) {
    const double *low, *mid, *high;
    const double *low2, *mid2, *high2;
    if (inputA > 3 && inputA <= 10) {
      low2 = lowEpolyHe4w2; mid2 = midEpolyHe4w2; high2 = highEpolyHe4w2;
      low = protonE1polyHe4; mid = protonE2polyHe4; high = protonE3polyHe4;
      if (pdg == 2112) { low = neutronE1polyHe4; mid = neutronE2polyHe4; high = neutronE3polyHe4; }
    } else if (inputA > 10 && inputA <= 15) {
      low2 = lowEpolyC12w2; mid2 = midEpolyC12w2; high2 = highEpolyC12w2;
      low = protonE1polyC12; mid = protonE2polyC12; high = protonE3polyC12;
      if (pdg == 2112) { low = neutronE1polyC12; mid = neutronE2polyC12; high = neutronE3polyC12; }
    } else if (inputA > 15 && inputA <= 23) {
      low2 = lowEpolyO16w2; mid2 = midEpolyO16w2; high2 = highEpolyO16w2;
      low = protonE1polyO16; mid = protonE2polyO16; high = protonE3polyO16;
      if (pdg == 2112) { low = neutronE1polyO16; mid = neutronE2polyO16; high = neutronE3polyO16; }
    } else if (inputA > 23 && inputA < 32) {
      low2 = lowEpolySi28w2; mid2 = midEpolySi28w2; high2 = highEpolySi28w2;
      low = protonE1polySi28; mid = protonE2polySi28; high = protonE3polySi28;
      if (pdg == 2112) { low = neutronE1polySi28; mid = neutronE2polySi28; high = neutronE3polySi28; }
    } else if (inputA > 32 && inputA <= 48) {  // QUIRK: A == 32 falls through to "not tested"
      low2 = lowEpolyAr40w2; mid2 = midEpolyAr40w2; high2 = highEpolyAr40w2;
      low = protonE1polyAr40; mid = protonE2polyAr40; high = protonE3polyAr40;
      if (pdg == 2112) { low = neutronE1polyAr40; mid = neutronE2polyAr40; high = neutronE3polyAr40; }
    } else if (inputA > 48 && inputA <= 70) {
      low2 = lowEpolyFe56w2; mid2 = midEpolyFe56w2; high2 = highEpolyFe56w2;
      low = protonE1polyFe56; mid = protonE2polyFe56; high = protonE3polyFe56;
      if (pdg == 2112) { low = neutronE1polyFe56; mid = neutronE2polyFe56; high = neutronE3polyFe56; }
    } else if (inputA > 195 && inputA <= 240) {
      low2 = lowEpolyPb208w2; mid2 = midEpolyPb208w2; high2 = highEpolyPb208w2;
      low = protonE1polyPb208; mid = protonE2polyPb208; high = protonE3polyPb208;
      if (pdg == 2112) { low = neutronE1polyPb208; mid = neutronE2polyPb208; high = neutronE3polyPb208; }
    } else {
      // QUIRK: untested nuclei return errorcode 1 with weights untouched (caller's 1.0)
      return 1.0;
    }
    if (inputKEgev > m_highcutoff) inputKEgev = m_highcutoff;
    const double *poly = high, *poly2 = high2;
    if (inputKEgev < m_breakpointlow) { poly = low; poly2 = low2; }
    else if (inputKEgev < m_breakpointhigh) { poly = mid; poly2 = mid2; }
    double powke = 1.0, tempweight1 = 0.0, tempweight2 = 0.0;
    for (int i = 0; i < 10; i++) {
      tempweight1 += powke * poly[i];
      tempweight2 += powke * poly2[i];
      powke *= inputKEgev;
    }
    weight1 = tempweight1;
    weight2 = tempweight2;
  }

  if (pdg == 211 || pdg == 111 || pdg == -211) {
    // Pions: carbon parameterisation for all nuclei, four KE ranges.
    const double *vlow, *low, *mid, *high;
    if (pdg == 211) {
      vlow = piplusE1polyC12; low = piplusE2polyC12; mid = piplusE3polyC12; high = piplusE4polyC12;
    } else if (pdg == 111) {
      vlow = pizeroE1polyC12; low = pizeroE2polyC12; mid = pizeroE3polyC12; high = pizeroE4polyC12;
    } else {
      vlow = piminusE1polyC12; low = piminusE2polyC12; mid = piminusE3polyC12; high = piminusE4polyC12;
    }
    const double breakpointvlow = 0.195, breakpointlow = 0.625, breakpointhigh = 1.200, highcutoff = 3.8;
    if (inputKEgev > highcutoff) inputKEgev = highcutoff;
    const double* poly = high;
    if (inputKEgev < breakpointvlow) poly = vlow;
    else if (inputKEgev < breakpointlow) poly = low;
    else if (inputKEgev < breakpointhigh) poly = mid;
    double powke = 1.0, tempweight = 0.0;
    for (int i = 0; i < 10; i++) {
      tempweight += powke * poly[i];
      powke *= inputKEgev;
    }
    weight1 = tempweight;
    weight2 = 1.0;
  }
  return 0;
}

}  // namespace fsi

// ---------------------------------------------------------------------------
void FSIWeights::request(BranchSet& bs) {
  fCharm = bs.require("mc_charm");
  fIntType = bs.require("mc_intType");
  fA = bs.require("mc_targetA");
  fZ = bs.require("mc_targetZ");
  fResID = bs.require("mc_resID");
  fN = bs.require("mc_er_nPart");
  fID = bs.require("mc_er_ID");
  fStatus = bs.require("mc_er_status");
  fFD = bs.require("mc_er_FD");
  fLD = bs.require("mc_er_LD");
  fPx = bs.require("mc_er_Px");
  fPy = bs.require("mc_er_Py");
  fPz = bs.require("mc_er_Pz");
  fE = bs.require("mc_er_E");
  fCalc.UseTrackingThreshold();  // as in FSIReweighter
}

void FSIWeights::book(TTree* t) {
  BindOutput(t, "w_fsi_elastic", &w_elastic, 'F');
  BindOutput(t, "w_fsi_absorption", &w_absorption, 'F');
}

void FSIWeights::compute() {
  auto I = [](const InBranch* b, size_t i) { return (int)std::lround(b->get(i)); };
  // all mc_er_* arrays share the mc_er_nPart count; clamp to the read buffers
  size_t n = (size_t)std::max<long>(0, std::lround(fN->get(0)));
  for (const InBranch* b : {fID, fStatus, fFD, fLD, fPx, fPy, fPz, fE}) n = std::min(n, b->size());
  vID.resize(n); vStatus.resize(n); vFD.resize(n); vLD.resize(n);
  vPx.resize(n); vPy.resize(n); vPz.resize(n); vE.resize(n);
  for (size_t i = 0; i < n; ++i) {
    vID[i] = I(fID, i); vStatus[i] = I(fStatus, i); vFD[i] = I(fFD, i); vLD[i] = I(fLD, i);
    vPx[i] = fPx->get(i); vPy[i] = fPy->get(i); vPz[i] = fPz->get(i); vE[i] = fE->get(i);
  }
  fCalc.calcWeights(I(fCharm, 0), I(fIntType, 0), I(fA, 0), I(fZ, 0), I(fResID, 0), (int)n, vID.data(),
                    vStatus.data(), vFD.data(), vLD.data(), vPx.data(), vPy.data(), vPz.data(), vE.data());
  w_elastic = (float)fCalc.GetElasticWeight(1);
  w_absorption = (float)fCalc.GetAbsorptionWeight();
}

}  // namespace makedf
