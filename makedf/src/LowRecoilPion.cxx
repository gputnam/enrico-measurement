#include "LowRecoilPion.h"

#include <algorithm>
#include <cmath>
#include <utility>

#include "TTree.h"

namespace makedf {
namespace lrp {

namespace {
inline double at(const InBranch* b, size_t i) { return (b && i < b->size()) ? b->get(i) : NAN; }
inline size_t nElem(const InBranch* b) { return b ? b->size() : 0; }
}  // namespace

double TpiFromRange(double range, double p0, double p1) { return p0 * range + p1 * std::sqrt(range); }

double ThetaWRTBeam(double x, double y, double z) {
  const double pyp = -1.0 * std::sin(kNuMIBeamAngle) * z + std::cos(kNuMIBeamAngle) * y;
  const double pzp = std::cos(kNuMIBeamAngle) * z + std::sin(kNuMIBeamAngle) * y;
  const double denom2 = x * x + pyp * pyp + pzp * pzp;
  if (0. == denom2) return -9999.;
  return std::acos(pzp / std::sqrt(denom2));
}

// ---------------------------------------------------------------------------
void Inputs::request(BranchSet& bs) {
  vtx = bs.require("vtx");
  fitPass = bs.require("FittedMichel_michel_fitPass");
  energy = bs.require("FittedMichel_michel_energy");
  time = bs.require("FittedMichel_michel_time");
  datafrac = bs.require("FittedMichel_michel_datafraction");
  x1 = bs.require("FittedMichel_michel_x1"); y1 = bs.require("FittedMichel_michel_y1");
  u1 = bs.require("FittedMichel_michel_u1"); v1 = bs.require("FittedMichel_michel_v1");
  z1 = bs.require("FittedMichel_michel_z1");
  x2 = bs.require("FittedMichel_michel_x2"); y2 = bs.require("FittedMichel_michel_y2");
  u2 = bs.require("FittedMichel_michel_u2"); v2 = bs.require("FittedMichel_michel_v2");
  z2 = bs.require("FittedMichel_michel_z2");
  clusViewSz = bs.require("cluster_view_sz");
  clusView = bs.require("cluster_view");
  clusSubdet = bs.require("cluster_subdet");
  clusIsMu = bs.require("cluster_isMuontrack");
  clusEnergy = bs.require("cluster_energy");
  clusPos = bs.require("cluster_pos");
  clusTime = bs.require("cluster_time");
  clusZ = bs.require("cluster_z");
  // MC-only Michel truth
  tParentPDG = bs.request("truth_FittedMichel_true_primaryparent_pdg");
  tParentE = bs.request("truth_FittedMichel_true_primaryparent_energy");
  tParentP = bs.request("truth_FittedMichel_true_primaryparent_momentum");
  tParentPx = bs.request("truth_FittedMichel_true_primaryparent_momentumx");
  tParentPy = bs.request("truth_FittedMichel_true_primaryparent_momentumy");
  tParentPz = bs.request("truth_FittedMichel_true_primaryparent_momentumz");
  tMichelInitZ = bs.request("truth_FittedMichel_reco_micheltrajectory_initialz");
  isMC = tParentPDG && tParentE && tParentP && tParentPx && tParentPy && tParentPz && tMichelInitZ;
}

Vertex Inputs::vertex() const {
  Vertex v;
  v.x = vtx->get(0); v.y = vtx->get(1); v.z = vtx->get(2);
  v.t = vtx->get(3) / std::pow(10, 3);  // ns -> us
  return v;
}

int Inputs::nMichelCandidates() const { return (int)std::lround(fitPass->count->get(0)); }

// ---------------------------------------------------------------------------
// LowRecoilPion::Michel<T>::Michel(univ, ci), without the matching calls
void Michel::Init(const Inputs& in, int ci, bool fillTruth) {
  tuple_idx = ci;  // EverChanges branch: Michel(...) : tuple_idx(ci)
  energy = at(in.energy, ci);
  time = at(in.time, ci) / std::pow(10, 3);
  is_fitted = (int)at(in.fitPass, ci);
  m_x1 = at(in.x1, ci); m_y1 = at(in.y1, ci); m_u1 = at(in.u1, ci); m_v1 = at(in.v1, ci); m_z1 = at(in.z1, ci);
  m_x2 = at(in.x2, ci); m_y2 = at(in.y2, ci); m_u2 = at(in.u2, ci); m_v2 = at(in.v2, ci); m_z2 = at(in.z2, ci);
  nclusters = (int)std::lround(in.clusViewSz->get(0));
  overlay_fraction = at(in.datafrac, ci);

  if (fillTruth) {
    true_initialz = at(in.tMichelInitZ, ci);
    true_parentpdg = (int)at(in.tParentPDG, ci);
    const double pe = at(in.tParentE, ci), pp = at(in.tParentP, ci);
    const double mass = std::sqrt(pe * pe - pp * pp);
    pionKE = pe - mass;
    true_angle = ThetaWRTBeam(at(in.tParentPx, ci), at(in.tParentPy, ci), at(in.tParentPz, ci));
    const double end1diff = std::fabs(true_initialz - m_z1);
    const double end2diff = std::fabs(true_initialz - m_z2);
    if (overlay_fraction > 0.5) trueEndpoint = 0;
    else if (true_parentpdg == 211 && end1diff < end2diff) trueEndpoint = 1;
    else if (true_parentpdg == 211 && end2diff < end1diff) trueEndpoint = 2;
  }
}

void Michel::DoesMichelMatchVtx(const Vertex& v) {
  const double vtx_u = (0.5 * (v.x - std::sqrt(3.) * v.y));
  const double vtx_v = (0.5 * (v.x + std::sqrt(3.) * v.y));
  const double zdiff1 = v.z - m_z1;
  const double zdiff2 = v.z - m_z2;
  vtx_michel_timediff = time - v.t;

  double xdiff = std::fabs(v.x - m_x1), udiff = std::fabs(vtx_u - m_u1), vdiff = std::fabs(vtx_v - m_v1);
  up_to_vertex_XZ = std::sqrt(xdiff * xdiff + zdiff1 * zdiff1);
  up_to_vertex_UZ = std::sqrt(udiff * udiff + zdiff1 * zdiff1);
  up_to_vertex_VZ = std::sqrt(vdiff * vdiff + zdiff1 * zdiff1);

  xdiff = std::fabs(v.x - m_x2); udiff = std::fabs(vtx_u - m_u2); vdiff = std::fabs(vtx_v - m_v2);
  down_to_vertex_XZ = std::sqrt(xdiff * xdiff + zdiff2 * zdiff2);
  down_to_vertex_UZ = std::sqrt(udiff * udiff + zdiff2 * zdiff2);
  down_to_vertex_VZ = std::sqrt(vdiff * vdiff + zdiff2 * zdiff2);

  const double xdiff1 = std::fabs(v.x - m_x1), xdiff2 = std::fabs(v.x - m_x2);
  const double ydiff1 = std::fabs(v.y - m_y1), ydiff2 = std::fabs(v.y - m_y2);
  const double dist1 = std::sqrt(zdiff1 * zdiff1 + xdiff1 * xdiff1 + ydiff1 * ydiff1);
  const double dist2 = std::sqrt(zdiff2 * zdiff2 + xdiff2 * xdiff2 + ydiff2 * ydiff2);
  up_to_vertex_dist3D = dist1;
  down_to_vertex_dist3D = dist2;
  if (dist1 < dist2) vtx_endpoint = 1;
  else if (dist2 < dist1) vtx_endpoint = 2;
}

void Michel::DoesMichelMatchClus(const Inputs& in, const Vertex& v) {
  int ncl = (int)std::lround(in.clusViewSz->get(0));
  // guard against inconsistent array lengths (never expected)
  for (const InBranch* b : {in.clusView, in.clusSubdet, in.clusIsMu, in.clusEnergy, in.clusPos, in.clusTime, in.clusZ})
    ncl = std::min<int>(ncl, (int)nElem(b));

  const double michelx1 = m_x1, michelx2 = m_x2, michelu1 = m_u1, michelu2 = m_u2;
  const double michelv1 = m_v1, michelv2 = m_v2, michelz1 = m_z1, michelz2 = m_z2;
  const double michely1 = m_y1, michely2 = m_y2;
  const double micheltime = time;

  double closestdistance1x = 9999., closestdistance1u = 9999., closestdistance1v = 9999.;
  double closestdistance2x = 9999., closestdistance2u = 9999., closestdistance2v = 9999.;
  int x1_idx = -1, u1_idx = -1, v1_idx = -1, x2_idx = -1, u2_idx = -1, v2_idx = -1;

  // First loop over clusters: closest cluster per view per endpoint
  for (int i = 0; i < ncl; i++) {
    const int subdet = (int)in.clusSubdet->get(i);
    if (subdet != 2 and subdet != 3) continue;
    const int ismuon = (int)in.clusIsMu->get(i);
    if (ismuon != 0) continue;
    const double energy_ = in.clusEnergy->get(i);
    if (energy_ < 2.) continue;
    const double pos = in.clusPos->get(i);
    if (std::fabs(pos) > 1000) continue;
    const double ctime = in.clusTime->get(i) / std::pow(10, 3);
    const double timediff = micheltime - ctime;
    if (timediff < 0.) continue;
    const double zpos = in.clusZ->get(i);
    const int view = (int)in.clusView->get(i);
    const double zdiff1 = std::fabs(zpos - michelz1);
    const double zdiff2 = std::fabs(zpos - michelz2);
    if (view == 1) {
      const double xdiff1 = std::fabs(pos - michelx1);
      const double d1 = std::sqrt(xdiff1 * xdiff1 + zdiff1 * zdiff1);
      const double xdiff2 = std::fabs(pos - michelx2);
      const double d2 = std::sqrt(xdiff2 * xdiff2 + zdiff2 * zdiff2);
      if (d1 <= closestdistance1x) { closestdistance1x = d1; x1_idx = i; }
      if (d2 <= closestdistance2x) { closestdistance2x = d2; x2_idx = i; }
    } else if (view == 2) {
      const double udiff1 = std::fabs(pos - michelu1);
      const double d1 = std::sqrt(udiff1 * udiff1 + zdiff1 * zdiff1);
      const double udiff2 = std::fabs(pos - michelu2);
      const double d2 = std::sqrt(udiff2 * udiff2 + zdiff2 * zdiff2);
      if (d1 < closestdistance1u) { closestdistance1u = d1; u1_idx = i; }
      if (d2 < closestdistance2u) { closestdistance2u = d2; u2_idx = i; }
    } else if (view == 3) {
      const double vdiff1 = std::fabs(pos - michelv1);
      const double d1 = std::sqrt(vdiff1 * vdiff1 + zdiff1 * zdiff1);
      const double vdiff2 = std::fabs(pos - michelv2);
      const double d2 = std::sqrt(vdiff2 * vdiff2 + zdiff2 * zdiff2);
      if (d1 < closestdistance1v) { closestdistance1v = d1; v1_idx = i; }
      if (d2 <= closestdistance2v) { closestdistance2v = d2; v2_idx = i; }
    }
  }

  // Second "loop": only the closest clusters, with an extra z (tracker+ECAL) requirement.
  // Iterating over the (at most 6) indices in increasing order is equivalent to the MAT
  // loop over all clusters with a find() in closestidx.
  std::vector<double> clusx1, clusx2, clusu1, clusu2, clusv1, clusv2;
  int closestidx[6] = {u1_idx, u2_idx, x1_idx, x2_idx, v1_idx, v2_idx};
  std::sort(closestidx, closestidx + 6);
  int prev = -2;
  for (int k = 0; k < 6; ++k) {
    const int i = closestidx[k];
    if (i < 0 || i == prev) continue;
    prev = i;
    const double energy_ = in.clusEnergy->get(i);
    const double ctime = in.clusTime->get(i) / std::pow(10, 3);
    const double pos = in.clusPos->get(i);
    const double zpos = in.clusZ->get(i);
    const int view = (int)in.clusView->get(i);
    const double timediff = micheltime - ctime;
    const int ismuon = (int)in.clusIsMu->get(i);
    if (ismuon != 0) continue;
    if (energy_ < 2.) continue;
    if (timediff < 0.) continue;
    if (zpos < 5980 or zpos > 9038) continue;
    if (view == 1) {
      if (i == x1_idx) { clusx1.push_back(pos); clusx1.push_back(zpos); }
      if (i == x2_idx) { clusx2.push_back(pos); clusx2.push_back(zpos); }
    } else if (view == 2) {
      if (i == u1_idx) { clusu1.push_back(pos); clusu1.push_back(zpos); }
      if (i == u2_idx) { clusu2.push_back(pos); clusu2.push_back(zpos); }
    } else if (view == 3) {
      if (i == v1_idx) { clusv1.push_back(pos); clusv1.push_back(zpos); }
      if (i == v2_idx) { clusv2.push_back(pos); clusv2.push_back(zpos); }
    }
  }

  std::vector<double> matchclus1, matchclus2;  // [0]=x [1]=y [2]=z
  double XZdist1 = 9999., UZdist1 = 9999., VZdist1 = 9999.;
  if (!clusx1.empty()) { double a = std::fabs(m_x1 - clusx1[0]), b = std::fabs(m_z1 - clusx1[1]); XZdist1 = std::sqrt(a * a + b * b); }
  if (!clusu1.empty()) { double a = std::fabs(m_u1 - clusu1[0]), b = std::fabs(m_z1 - clusu1[1]); UZdist1 = std::sqrt(a * a + b * b); }
  if (!clusv1.empty()) { double a = std::fabs(m_v1 - clusv1[0]), b = std::fabs(m_z1 - clusv1[1]); VZdist1 = std::sqrt(a * a + b * b); }
  up_to_clus_XZ = XZdist1; up_to_clus_UZ = UZdist1; up_to_clus_VZ = VZdist1;

  double XZdist2 = 9999., UZdist2 = 9999., VZdist2 = 9999.;
  if (!clusx2.empty()) { double a = std::fabs(m_x2 - clusx2[0]), b = std::fabs(m_z2 - clusx2[1]); XZdist2 = std::sqrt(a * a + b * b); }
  if (!clusu2.empty()) { double a = std::fabs(m_u2 - clusu2[0]), b = std::fabs(m_z2 - clusu2[1]); UZdist2 = std::sqrt(a * a + b * b); }
  if (!clusv2.empty()) { double a = std::fabs(m_v2 - clusv2[0]), b = std::fabs(m_z2 - clusv2[1]); VZdist2 = std::sqrt(a * a + b * b); }
  down_to_clus_XZ = XZdist2; down_to_clus_UZ = UZdist2; down_to_clus_VZ = VZdist2;

  const double s3 = 1. / std::sqrt(3.);
  // 3D cluster point from the two closest views (MAT logic reproduced verbatim,
  // including the x/y swap in the "VU" branch).
  if (XZdist1 < UZdist1 && XZdist1 < VZdist1 && UZdist1 < VZdist1) {  // XU
    if (!clusu1.empty() && !clusx1.empty()) matchclus1 = {clusx1[0], s3 * (clusx1[0] - 2 * clusu1[0]), clusx1[1]};
  } else if (XZdist1 < UZdist1 && XZdist1 < VZdist1 && UZdist1 > VZdist1) {  // XV
    if (!clusv1.empty() && !clusx1.empty()) matchclus1 = {clusx1[0], s3 * (2 * clusv1[0] - clusx1[0]), clusx1[1]};
  } else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 && VZdist1 < XZdist1) {  // UV
    if (!clusv1.empty() && !clusu1.empty()) matchclus1 = {clusu1[0] + clusv1[0], s3 * (clusv1[0] - clusu1[0]), clusu1[1]};
  } else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 && VZdist1 > XZdist1) {  // UX
    if (!clusu1.empty() && !clusx1.empty()) matchclus1 = {clusx1[0], s3 * (clusx1[0] - 2 * clusu1[0]), clusu1[1]};
  } else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 && XZdist1 < UZdist1) {  // VX
    if (!clusv1.empty() && !clusx1.empty()) matchclus1 = {clusx1[0], s3 * (2 * clusv1[0] - clusx1[0]), clusv1[1]};
  } else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 && XZdist1 > UZdist1) {  // VU
    if (!clusu1.empty() && !clusv1.empty()) matchclus1 = {s3 * (clusv1[0] - clusu1[0]), clusu1[0] + clusv1[0], clusv1[1]};
  }
  if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 < VZdist2) {
    if (!clusu2.empty() && !clusx2.empty()) matchclus2 = {clusx2[0], s3 * (clusx2[0] - 2 * clusu2[0]), clusx2[1]};
  } else if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 > VZdist2) {
    if (!clusv2.empty() && !clusx2.empty()) matchclus2 = {clusx2[0], s3 * (2 * clusv2[0] - clusx2[0]), clusx2[1]};
  } else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 < XZdist2) {
    if (!clusu2.empty() && !clusv2.empty()) matchclus2 = {clusu2[0] + clusv2[0], s3 * (clusv2[0] - clusu2[0]), clusu2[1]};
  } else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 > XZdist2) {
    if (!clusu2.empty() && !clusx2.empty()) matchclus2 = {clusx2[0], s3 * (clusx2[0] - 2 * clusu2[0]), clusu2[1]};
  } else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 < UZdist2) {
    if (!clusv2.empty() && !clusx2.empty()) matchclus2 = {clusx2[0], s3 * (2 * clusv2[0] - clusx2[0]), clusv2[1]};
  } else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 > UZdist2) {
    if (!clusu2.empty() && !clusv2.empty()) matchclus2 = {s3 * (clusv2[0] - clusu2[0]), clusu2[0] + clusv2[0], clusv2[1]};
  }

  double clusx1diff = 9999., clusy1diff = 9999., clusz1diff = 9999.;
  double mclusx1diff = 9999., mclusy1diff = 9999., mclusz1diff = 9999.;
  double michvtx_x1diff = 9999., michvtx_x2diff = 9999., michvtx_y1diff = 9999.;
  double michvtx_y2diff = 9999., michvtx_z1diff = 9999., michvtx_z2diff = 9999.;
  if (!matchclus1.empty()) {
    clusx1diff = v.x - matchclus1[0]; clusy1diff = v.y - matchclus1[1]; clusz1diff = v.z - matchclus1[2];
    mclusx1diff = michelx1 - matchclus1[0]; mclusy1diff = michely1 - matchclus1[1]; mclusz1diff = michelz1 - matchclus1[2];
    michvtx_x1diff = michelx1 - v.x; michvtx_y1diff = michely1 - v.y; michvtx_z1diff = michelz1 - v.z;
  }
  double clusx2diff = 9999., clusy2diff = 9999., clusz2diff = 9999.;
  double mclusx2diff = 9999., mclusy2diff = 9999., mclusz2diff = 9999.;
  if (!matchclus2.empty()) {
    clusx2diff = v.x - matchclus2[0]; clusy2diff = v.y - matchclus2[1]; clusz2diff = v.z - matchclus2[2];
    mclusx2diff = michelx2 - matchclus2[0]; mclusy2diff = michely2 - matchclus2[1]; mclusz2diff = michelz2 - matchclus2[2];
    michvtx_x2diff = michelx2 - v.x; michvtx_y2diff = michely2 - v.y; michvtx_z2diff = michelz2 - v.z;
  }
  auto norm = [](double a, double b, double c) { return std::sqrt(std::pow(a, 2) + std::pow(b, 2) + std::pow(c, 2)); };
  const double dist1 = norm(clusx1diff, clusy1diff, clusz1diff);   // cluster - vertex
  const double dist2 = norm(clusx2diff, clusy2diff, clusz2diff);
  down_clus_michel_dist3D = norm(mclusx2diff, mclusy2diff, mclusz2diff);  // cluster - Michel
  up_clus_michel_dist3D = norm(mclusx1diff, mclusy1diff, mclusz1diff);
  up_to_cluster_dist3D = dist1;
  down_to_cluster_dist3D = dist2;
  up_clus_michvtx_dist3D = norm(michvtx_x1diff, michvtx_y1diff, michvtx_z1diff);  // Michel - vertex (pion range)
  down_clus_michvtx_dist3D = norm(michvtx_x2diff, michvtx_y2diff, michvtx_z2diff);
  if (dist1 < dist2) clus_endpoint = 1;
  else if (dist1 > dist2) clus_endpoint = 2;
}

void Michel::GetPionAngle(const Vertex& v) {
  double ex, ey, ez;
  if (BestMatch == 1 || BestMatch == 3) { ex = m_x1; ey = m_y1; ez = m_z1; }
  else if (BestMatch == 2 || BestMatch == 4) { ex = m_x2; ey = m_y2; ez = m_z2; }
  else { ex = ey = ez = 9999.; }
  const double rx = ex - v.x, ry = ey - v.y, rz = ez - v.z;
  const double angle = ThetaWRTBeam(rx, ry, rz);
  best_angle = angle;
  const double xyp = -1.0 * std::sin(kNuMIBeamAngle + 0.00000001) * rz + std::cos(kNuMIBeamAngle + 0.00000001) * ry;
  const double phi = std::atan2(xyp, rx);
  best_phi = phi;
  const double tpi = TpiFromRange(Best3Ddist);
  const double ppi = std::sqrt(2 * tpi * kMpiLRP + tpi * tpi);
  reco_Epi = tpi + kMpiLRP;
  reco_KE = tpi;
  reco_ppi = ppi;
  reco_ppiz = ppi * std::cos(angle);
  reco_ppix = ppi * std::sin(angle) * std::cos(phi);
  reco_ppiy = ppi * std::sin(angle) * std::sin(phi);
}

// ---------------------------------------------------------------------------
bool HasMichelCut(const Inputs& in, MichelEvent& evt, bool fillTruth) {
  const int nmichels = in.nMichelCandidates();
  if (nmichels < 1) return false;
  int nfitted = 0;
  for (size_t i = 0; i < in.fitPass->size(); ++i) nfitted += (std::lround(in.fitPass->get(i)) == 1);
  if (nfitted < 1) return false;
  const int nclusters = (int)std::lround(in.clusViewSz->get(0));
  const Vertex vtx = in.vertex();
  for (int i = 0; i < nmichels; ++i) {
    Michel m;
    m.Init(in, i, fillTruth);
    if (m.is_fitted != 1) continue;
    // MAT matches every fitted Michel to vtx+clusters in the ctor before these cuts; the
    // cluster matching is deferred here until after the cuts (identical result, faster).
    m.DoesMichelMatchVtx(vtx);
    if (std::fabs(m.vtx_michel_timediff) < 0.400) continue;  // dead-time Michels
    const double z1 = m.m_z1, z2 = m.m_z2;
    if (z1 < 5980. || z2 < 5980.) continue;  // tracker + ECAL
    else if (z1 > 9038. || z2 > 9038.) continue;
    if (m.nclusters > 0) m.DoesMichelMatchClus(in, vtx);
    evt.m_allmichels.push_back(m);
  }
  evt.nclusters = nclusters;
  return !evt.m_allmichels.empty();
}

bool BestMichelDistance2DCut(MichelEvent& evt, const double max_dist) {
  std::vector<Michel> nmichelspass;
  if (evt.m_allmichels.empty()) return false;
  for (unsigned int i = 0; i < evt.m_allmichels.size(); i++) {
    Michel& m = evt.m_allmichels[i];
    const double upvtxXZ = m.up_to_vertex_XZ, downvtxXZ = m.down_to_vertex_XZ;
    const double upvtxUZ = m.up_to_vertex_UZ, downvtxUZ = m.down_to_vertex_UZ;
    const double upvtxVZ = m.up_to_vertex_VZ, downvtxVZ = m.down_to_vertex_VZ;
    // "two of three views closer than max_dist"
    auto two = [&](double a, double b, double c) {
      return (a < max_dist && (b < max_dist || c < max_dist)) || (b < max_dist && (a < max_dist || c < max_dist)) ||
             (c < max_dist && (a < max_dist || b < max_dist));
    };
    m.passable_matchtype[0] = two(upvtxXZ, upvtxUZ, upvtxVZ) ? 1 : -1;
    m.passable_matchtype[1] = two(downvtxXZ, downvtxUZ, downvtxVZ) ? 2 : -1;
    m.passable_matchtype[2] = two(m.up_to_clus_XZ, m.up_to_clus_UZ, m.up_to_clus_VZ) ? 3 : -1;
    m.passable_matchtype[3] = two(m.down_to_clus_XZ, m.down_to_clus_UZ, m.down_to_clus_VZ) ? 4 : -1;
    if (m.passable_matchtype[1] == -1 && m.passable_matchtype[2] == -1 && m.passable_matchtype[3] == -1 &&
        m.passable_matchtype[0] == -1)
      continue;

    std::vector<double> distances3D;
    if (m.passable_matchtype[0] == 1) distances3D.push_back(m.up_to_vertex_dist3D);
    if (m.passable_matchtype[1] == 2) distances3D.push_back(m.down_to_vertex_dist3D);
    if (m.passable_matchtype[2] == 3) distances3D.push_back(m.up_clus_michel_dist3D);
    if (m.passable_matchtype[3] == 4) distances3D.push_back(m.down_clus_michel_dist3D);
    if (distances3D.empty()) distances3D = {9999., 9999., 9999., 9999.};
    std::sort(distances3D.begin(), distances3D.end());

    if (distances3D[0] == 9999.) continue;
    else if (distances3D[0] == m.up_to_vertex_dist3D) {
      m.best_XZ = m.up_to_vertex_XZ; m.best_UZ = m.up_to_vertex_UZ; m.best_VZ = m.up_to_vertex_VZ;
      m.BestMatch = 1; m.Best3Ddist = m.up_to_vertex_dist3D;
    } else if (distances3D[0] == m.down_to_vertex_dist3D) {
      m.BestMatch = 2;
      m.best_XZ = m.down_to_vertex_XZ; m.best_UZ = m.down_to_vertex_UZ; m.best_VZ = m.down_to_vertex_VZ;
      m.Best3Ddist = m.down_to_vertex_dist3D;
    } else if (distances3D[0] == m.up_clus_michel_dist3D) {
      m.BestMatch = 3;
      m.best_XZ = m.up_to_clus_XZ;
      m.best_UZ = m.up_to_clus_VZ;  // (sic) U/V swapped in MAT
      m.best_VZ = m.up_to_clus_UZ;
      m.Best3Ddist = m.up_clus_michvtx_dist3D;  // range = Michel endpoint to vertex
    } else if (distances3D[0] == m.down_clus_michel_dist3D) {
      m.BestMatch = 4;
      m.Best3Ddist = m.down_clus_michvtx_dist3D;
      m.best_XZ = m.down_to_clus_XZ; m.best_UZ = m.down_to_clus_UZ; m.best_VZ = m.down_to_clus_VZ;
    } else {
      m.BestMatch = 0; m.Best3Ddist = 9999.; m.best_XZ = m.best_UZ = m.best_VZ = 9999.;
      continue;
    }
    // MAT reads distances3D[1] even when size()==1 (out of bounds); guarded here.
    // SecondBestMatch is not used by any downstream decision.
    if (distances3D.size() < 2) m.SecondBestMatch = 0;
    else if (distances3D[1] == m.up_to_vertex_dist3D) m.SecondBestMatch = 1;
    else if (distances3D[1] == m.down_to_vertex_dist3D) m.SecondBestMatch = 2;
    else if (distances3D[1] == m.up_clus_michel_dist3D) m.SecondBestMatch = 3;
    else if (distances3D[1] == m.down_clus_michel_dist3D) m.SecondBestMatch = 4;
    else m.SecondBestMatch = 0;

    if (m.best_XZ < m.best_UZ and m.best_XZ < m.best_VZ) m.bestview = 1;
    else if (m.best_UZ < m.best_XZ and m.best_UZ < m.best_VZ) m.bestview = 2;
    else if (m.best_VZ < m.best_UZ and m.best_VZ < m.best_XZ) m.bestview = 3;
    if (m.BestMatch == 1 || m.BestMatch == 3) m.recoEndpoint = 1;
    else if (m.BestMatch == 2 || m.BestMatch == 4) m.recoEndpoint = 2;
    nmichelspass.push_back(m);
    evt.m_nmichelspass.push_back(m);
  }
  evt.m_nmichels.clear();
  if (!nmichelspass.empty()) {
    evt.selection = 1;
    evt.m_nmichels = nmichelspass;
    return true;
  }
  return false;
}

bool GetClosestMichelCut(const Inputs& in, MichelEvent& evt, const double /*max_dist (unused in MAT)*/) {
  const int nmichels = evt.m_nmichels.size();
  std::vector<Michel> closestMichel;
  if (nmichels == 0) return false;
  evt.m_bestdist = 9999.;
  std::vector<double> allmichel3Ddist;
  for (int i = 0; i < nmichels; ++i) allmichel3Ddist.push_back(evt.m_nmichels[i].Best3Ddist);
  std::sort(allmichel3Ddist.begin(), allmichel3Ddist.end());
  // (sic) MAT: order is only assigned when Michel i sits at sorted position i, so the
  // "closest" Michel is found only if it is the first one in the list.
  for (int i = 0; i < nmichels; ++i)
    if (evt.m_nmichels[i].Best3Ddist == allmichel3Ddist[i]) evt.m_nmichels[i].OrderOfMichel = i + 1;
  const Vertex vtx = in.vertex();
  for (int i = 0; i < nmichels; ++i) {
    Michel& m = evt.m_nmichels[i];
    m.GetPionAngle(vtx);
    const double dist = m.Best3Ddist;
    evt.m_nmichelspass.clear();
    if (m.OrderOfMichel == 1) {
      evt.m_bestdist = dist;
      evt.m_idx = i;
      if (m.overlay_fraction > 0.5) evt.ClosestMichelsIsOverlay = 1;
      evt.m_best_XZ = m.best_XZ; evt.m_best_UZ = m.best_UZ; evt.m_best_VZ = m.best_VZ;
      evt.m_matchtype = m.BestMatch;
      if (m.BestMatch == 1 || m.BestMatch == 3) { evt.best_x = m.m_x1; evt.best_y = m.m_y1; evt.best_z = m.m_z1; }
      else if (m.BestMatch == 2 || m.BestMatch == 4) { evt.best_x = m.m_x2; evt.best_y = m.m_y2; evt.best_z = m.m_z2; }
      closestMichel.push_back(m);
    } else {
      evt.m_idx = -1;
    }
  }
  if (closestMichel.empty()) return false;
  evt.lowTpi = closestMichel[0].pionKE;
  evt.m_nmichels = closestMichel;
  // EverChanges: m_bestthetaangle = m_nmichels[m_idx].best_angle. m_idx is 0 only when a
  // single Michel passed the 2D cut; for m_idx == -1 MAT reads out of bounds (undefined),
  // here NaN.
  evt.m_bestthetaangle = (evt.m_idx >= 0 && evt.m_idx < (int)evt.m_nmichels.size())
                             ? evt.m_nmichels[evt.m_idx].best_angle : NAN;
  return !(closestMichel[0].Best3Ddist > 2600.);
}

void TrackerECALMuFuzz(const InBranch* ids, const InBranch* es, double& trk, double& ecal) {
  trk = ecal = 0.;
  const int nfuzz = (int)std::min(ids->size(), es->size());
  for (int i = 0; i < nfuzz; i++) {
    const int planeID = (int)ids->get(i);
    if (planeID < 1504968704 || planeID > 1709703168) continue;
    const double fuzze = es->get(i);
    if (planeID > 1504968704 and planeID < 1560805376) trk += fuzze;
    else if (planeID > 1700003840 and planeID < 1709703168) ecal += fuzze;
  }
}

// ---------------------------------------------------------------------------
namespace endpoint {

void Inputs::request(BranchSet& bs) {
  idx = bs.require("matched_michel_idx");
  endDist = bs.require("matched_michel_end_dist");
  avgDist = bs.require("matched_michel_avg_dist");
  ovDist = bs.require("matched_michel_ov_dist");
}

static double GetDistMichel(const Inputs& in, EMatchCategory s, unsigned vtx) {
  const InBranch* b = s == kFit ? in.endDist : s == kNoFit ? in.avgDist : s == kOV ? in.ovDist : nullptr;
  if (!b) return -1.;
  double d = at(b, vtx) / 10.;  // mm -> cm (NaN if out of range -> invalid below)
  const bool valid = !std::isnan(d) && (d == 0. || std::fabs(d) > 0.0001) && std::fabs(d) < 500;
  return valid ? d : -1.;
}

static Michel MakeMichel(const Inputs& in, int i, int v) {
  Michel m;
  m.idx = i; m.vtx = v; m.had_idx = v - 1;
  const bool isIntVtx = (v == 0);
  const double fit = GetDistMichel(in, kFit, v), nofit = GetDistMichel(in, kNoFit, v), ov = GetDistMichel(in, kOV, v);
  const double FIT_CUT = isIntVtx ? 9.0 : 7.5;     // cm
  const double NOFIT_CUT = isIntVtx ? 10.0 : 50.;  // cm
  if (fit > 0 && fit < FIT_CUT) { m.match_category = kFit; m.fit_distance = fit; }
  else if (nofit > 0 && nofit < NOFIT_CUT) { m.match_category = kNoFit; m.fit_distance = nofit; }
  else if (ov > 0 && ov * (2. / 3.) < NOFIT_CUT) { m.match_category = kOV; m.fit_distance = ov; }
  else { m.match_category = kNoMatch; m.fit_distance = -2.; }
  return m;
}

static Michel CompareMichels(const Michel& r, const Michel& c) {
  if (r.match_category > c.match_category) return r;  // (sic) larger enum wins
  if (r.match_category < c.match_category) return c;
  if (r.fit_distance < c.fit_distance) return r;
  if (r.fit_distance > c.fit_distance) return c;
  return r;
}

static void AddOrReplaceMichel(std::map<int, Michel>& mm, const Michel& m) {
  auto it = mm.find(m.idx);
  if (it == mm.end()) mm.emplace(m.idx, m);
  else it->second = CompareMichels(it->second, m);
}

std::map<int, Michel> GetQualityMichels(const Inputs& in) {
  std::map<int, Michel> ret;
  const size_t n = in.idx->size();
  for (unsigned vtx = 0; vtx < n; ++vtx) {
    const int mm_idx = (int)std::lround(in.idx->get(vtx));
    if (mm_idx < 0) continue;
    Michel cur = MakeMichel(in, mm_idx, vtx);
    if (cur.match_category == kNoMatch) continue;
    if (vtx == 0) continue;  // interaction-vertex Michels: trackless territory
    if (ret.size() == 0) {
    } else if (ret.size() == 1) {
      const Michel reigning = ret.begin()->second;
      if (reigning.match_category == kOV || cur.match_category == kOV) ret.clear();
      cur = CompareMichels(reigning, cur);
    } else {
      if (cur.match_category == kOV) continue;
    }
    AddOrReplaceMichel(ret, cur);
  }
  return ret;
}

}  // namespace endpoint
}  // namespace lrp

// ---------------------------------------------------------------------------
void LowRecoilPionSummary::request(BranchSet& bs) {
  fIn.request(bs);
  fEp.request(bs);
  fFuzzID = bs.require("muon_fuzz_per_plane_r80_planeIDs");
  fFuzzE = bs.require("muon_fuzz_per_plane_r80_energies");
  fRecTrk = bs.require("blob_recoil_E_tracker");
  fRecEcal = bs.require("blob_recoil_E_ecal");
}

void LowRecoilPionSummary::book(TTree* t) {
  auto I = [&](const std::string& n, int* v) { BindOutput(t, n, v, 'I'); };
  auto F = [&](const std::string& n, float* v) { BindOutput(t, n, v, 'F'); };
  I("lrp_ncand", &ncand); I("lrp_nfitted", &nfitted); I("lrp_nclusters", &nclusters);
  I("lrp_nall", &nall); I("lrp_npass2d", &npass2d); I("lrp_good", &good);
  I("lrp_idx", &idx); I("lrp_tuple_idx", &tuple_idx); I("lrp_matchtype", &matchtype); I("lrp_bestview", &bestview);
  F("lrp_range", &range); F("lrp_tpi", &tpi); F("lrp_theta", &theta); F("lrp_phi", &phi);
  F("lrp_best_x", &best_x); F("lrp_best_y", &best_y); F("lrp_best_z", &best_z);
  F("lrp_michel_e", &michel_e); F("lrp_michel_dt", &michel_dt); F("lrp_ovfrac", &ovfrac);
  F("lrp_mufuzz_trk", &mufuzz_trk); F("lrp_mufuzz_ecal", &mufuzz_ecal); F("lrp_eavail", &eavail);
  I("lrp_trk_nmich", &trk_nmich); I("lrp_trk_hadmask", &trk_hadmask);
  if (fIn.isMC) {
    I("lrp_t_parentpdg", &t_parentpdg); I("lrp_t_endpoint", &t_endpoint);
    F("lrp_t_parent_tpi", &t_parent_tpi); F("lrp_t_parent_theta", &t_parent_theta);
  }
}

void LowRecoilPionSummary::compute() {
  fEvt = lrp::MichelEvent();
  ncand = fIn.nMichelCandidates();
  nfitted = 0;
  for (size_t i = 0; i < fIn.fitPass->size(); ++i) nfitted += (std::lround(fIn.fitPass->get(i)) == 1);
  nclusters = (int)std::lround(fIn.clusViewSz->get(0));

  // Same call sequence as CC-CH-pip-ana xsec/makeCrossSectionMCInputs.C (short-circuit &&)
  bool ok = lrp::HasMichelCut(fIn, fEvt, fIn.isMC);
  nall = fEvt.m_allmichels.size();
  ok = ok && lrp::BestMichelDistance2DCut(fEvt);
  npass2d = ok ? (int)fEvt.m_nmichels.size() : 0;
  ok = ok && lrp::GetClosestMichelCut(fIn, fEvt);
  good = ok;

  const bool hasClosest = npass2d > 0 && !fEvt.m_nmichels.empty() && fEvt.m_nmichels[0].OrderOfMichel == 1;
  if (hasClosest) {
    const lrp::Michel& m = fEvt.m_nmichels[0];
    idx = fEvt.m_idx; tuple_idx = m.tuple_idx; matchtype = m.BestMatch; bestview = m.bestview;
    range = fEvt.m_bestdist; tpi = lrp::TpiFromRange(fEvt.m_bestdist);
    theta = m.best_angle; phi = m.best_phi;
    best_x = fEvt.best_x; best_y = fEvt.best_y; best_z = fEvt.best_z;
    michel_e = m.energy; michel_dt = m.vtx_michel_timediff; ovfrac = m.overlay_fraction;
    t_parentpdg = m.true_parentpdg; t_endpoint = m.trueEndpoint;
    t_parent_tpi = m.pionKE; t_parent_theta = m.true_angle;
  } else {
    idx = tuple_idx = matchtype = bestview = -1;
    range = tpi = theta = phi = best_x = best_y = best_z = michel_e = michel_dt = ovfrac = NAN;
    t_parentpdg = t_endpoint = -1;
    t_parent_tpi = t_parent_theta = NAN;
  }

  double ftrk, fecal;
  lrp::TrackerECALMuFuzz(fFuzzID, fFuzzE, ftrk, fecal);
  mufuzz_trk = ftrk; mufuzz_ecal = fecal;
  eavail = ((fRecTrk->get(0) - ftrk) + (fRecEcal->get(0) - fecal)) * lrp::kEavailScale;

  const auto em = lrp::endpoint::GetQualityMichels(fEp);
  trk_nmich = em.size();
  trk_hadmask = 0;
  for (const auto& kv : em)
    if (kv.second.had_idx >= 0 && kv.second.had_idx < 31) trk_hadmask |= (1 << kv.second.had_idx);
}

}  // namespace makedf
