// makedf: flatten MINERvA MasterAnaDev AnaTuples into preselected dataframe trees.
//
// Output trees:
//   reco  : one row per MasterAnaDev entry passing the reco preselection
//           (columns from reco.cols; MC adds reco_mc.cols, truth.cols and the truth summary)
//   truth : (MC only) one row per Truth-tree entry passing the truth preselection
//           (columns from truth.cols + truth summary)
//   meta  : one row per input file (POT, entries, file name)
//
// Usage: makedf -o out.root [-c configdir] [-l filelist] [-n maxentries] [--no-presel] [--no-truth] files...
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "Columns.h"
#include "Derived.h"
#include "FSIWeight.h"
#include "HadronReweight.h"
#include "LowRecoilPion.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

using namespace makedf;

namespace {

struct Options {
  std::string out, configDir, tag;
  std::vector<std::string> inputs;
  long maxEntries = -1;
  bool noTruth = false, strict = false;
  PreselConfig presel;
};

void Usage() {
  std::cerr << "Usage: makedf -o out.root [-c configdir] [-l filelist] [-n maxentries] [-t tag]\n"
               "              [--no-presel] [--no-truth] [--strict] input.root ...\n";
}

Options Parse(int argc, char** argv) {
  Options o;
  const char* env = std::getenv("MAKEDF_CONFIG");
  o.configDir = env ? env : "";
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) { Usage(); std::exit(2); }
      return argv[++i];
    };
    if (a == "-o") o.out = next();
    else if (a == "-c") o.configDir = next();
    else if (a == "-n") o.maxEntries = std::stol(next());
    else if (a == "-t") o.tag = next();
    else if (a == "-l") {
      std::ifstream in(next());
      std::string f;
      while (in >> f) if (!f.empty() && f[0] != '#') o.inputs.push_back(f);
    } else if (a == "--no-presel") o.presel.enabled = false;
    else if (a == "--no-truth") o.noTruth = true;
    else if (a == "--strict") o.strict = true;
    else if (a == "-h" || a == "--help") { Usage(); std::exit(0); }
    else o.inputs.push_back(a);
  }
  if (o.out.empty() || o.inputs.empty() || o.configDir.empty()) { Usage(); std::exit(2); }
  return o;
}

TFile* OpenWithRetry(const std::string& path, int tries = 4) {
  for (int t = 0; t < tries; ++t) {
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (f && !f->IsZombie()) return f;
    delete f;
    std::cerr << "makedf: failed to open " << path << " (attempt " << t + 1 << ")\n";
    std::this_thread::sleep_for(std::chrono::seconds(10 * (t + 1)));
  }
  return nullptr;
}

// Processes one tree into `out`; returns number of rows written.
long ProcessTree(TTree* in, TTree* out, const std::vector<std::vector<ColumnSpec>*>& specs, bool isTruthTree,
                 bool withSummary, const Options& o, int fileId, int& fileIdOut, HadronReweight* hr) {
  BranchSet bs(in);
  Preselection presel(o.presel);
  if (isTruthTree) presel.requestTruth(bs);
  else presel.requestReco(bs);

  // Output branches are created for the first input file and re-bound for later ones.
  std::vector<std::unique_ptr<OutColumn>> cols;
  for (auto* s : specs) {
    auto c = MakeColumns(*s, bs, out, o.strict);
    for (auto& x : c) cols.push_back(std::move(x));
  }
  TruthSummary summary;
  FSIWeights fsi;
  std::unique_ptr<GeantWeights> geant;
  if (withSummary) {  // MC: truth summary and truth-level reweights
    summary.request(bs);
    summary.book(out);
    fsi.request(bs);
    fsi.book(out);
    if (hr) {
      geant = std::make_unique<GeantWeights>(*hr, isTruthTree);
      geant->request(bs);
      geant->book(out);
    }
  }
  LowRecoilPionSummary lrp;  // untracked-pion (Michel) reconstruction summary, reco rows only
  if (!isTruthTree) {
    lrp.request(bs);
    lrp.book(out);
  }
  BindOutput(out, "file_id", &fileIdOut, 'I');
  fileIdOut = fileId;
  bs.bind();

  long n = in->GetEntries();
  if (o.maxEntries >= 0 && n > o.maxEntries) n = o.maxEntries;
  long written = 0;
  for (long i = 0; i < n; ++i) {
    if (in->GetEntry(i) <= 0) {
      std::cerr << "makedf: ERROR reading entry " << i << " of " << in->GetName() << "\n";
      std::exit(4);
    }
    if (isTruthTree ? !presel.passTruth() : !presel.passReco()) continue;
    for (auto& c : cols) c->fill();
    if (withSummary) {
      summary.compute();
      fsi.compute();
      if (geant) geant->compute();
    }
    if (!isTruthTree) lrp.compute();
    out->Fill();
    ++written;
  }
  // detach output addresses from soon-to-be-destroyed buffers
  out->ResetBranchAddresses();
  return written;
}

}  // namespace

int main(int argc, char** argv) {
  Options o = Parse(argc, argv);
  const std::string cd = o.configDir + "/";
  auto recoSpec = ReadColumnFile(cd + "reco.cols");
  auto recoMCSpec = ReadColumnFile(cd + "reco_mc.cols");
  auto truthSpec = ReadColumnFile(cd + "truth.cols");

  TFile fout(o.out.c_str(), "RECREATE", "", 505);  // ZSTD level 5
  TTree* treco = new TTree("reco", "preselected reco candidates");
  TTree* ttruth = new TTree("truth", "preselected truth interactions");
  TTree* tmeta = new TTree("meta", "per-input-file metadata");
  treco->SetAutoSave(0);
  ttruth->SetAutoSave(0);

  int fileId = 0, isMC = 0, fileIdOut = 0;
  double potUsed = 0, potTotal = 0;
  long long nReco = 0, nTruth = 0, nRecoPass = 0, nTruthPass = 0;
  char fname[2048];
  char tag[256];
  tmeta->Branch("file_id", &fileId, "file_id/I");
  tmeta->Branch("is_mc", &isMC, "is_mc/I");
  tmeta->Branch("pot_used", &potUsed, "pot_used/D");
  tmeta->Branch("pot_total", &potTotal, "pot_total/D");
  tmeta->Branch("n_reco", &nReco, "n_reco/L");
  tmeta->Branch("n_truth", &nTruth, "n_truth/L");
  tmeta->Branch("n_reco_pass", &nRecoPass, "n_reco_pass/L");
  tmeta->Branch("n_truth_pass", &nTruthPass, "n_truth_pass/L");
  tmeta->Branch("file", fname, "file/C");
  tmeta->Branch("tag", tag, "tag/C");
  snprintf(tag, sizeof(tag), "%s", o.tag.c_str());

  std::unique_ptr<HadronReweight> hadronRW;  // created on the first MC file
  int firstMC = -1;
  int nfail = 0;
  for (size_t ifile = 0; ifile < o.inputs.size(); ++ifile) {
    const auto& path = o.inputs[ifile];
    auto t0 = std::chrono::steady_clock::now();
    std::unique_ptr<TFile> fin(OpenWithRetry(path));
    if (!fin) { ++nfail; continue; }
    TTree* mad = (TTree*)fin->Get("MasterAnaDev");
    TTree* truth = (TTree*)fin->Get("Truth");
    TTree* meta = (TTree*)fin->Get("Meta");
    if (!mad || !meta) {
      std::cerr << "makedf: ERROR missing MasterAnaDev/Meta in " << path << "\n";
      ++nfail;
      continue;
    }
    isMC = truth != nullptr;
    if (firstMC < 0) firstMC = isMC;
    if (firstMC != isMC) { std::cerr << "makedf: ERROR mixing data and MC inputs\n"; return 5; }

    potUsed = potTotal = 0;
    if (meta->GetBranch("POT_Used")) { meta->SetBranchAddress("POT_Used", &potUsed); }
    if (meta->GetBranch("POT_Total")) { meta->SetBranchAddress("POT_Total", &potTotal); }
    double sumUsed = 0, sumTotal = 0;
    for (long i = 0; i < meta->GetEntries(); ++i) { meta->GetEntry(i); sumUsed += potUsed; sumTotal += potTotal; }
    potUsed = sumUsed;
    potTotal = sumTotal;

    if (isMC && !hadronRW) {
      HadronReweightConfig hc;
      hc.dataDir = ResolveHadronReweightDataDir("", o.configDir);
      hadronRW = std::make_unique<HadronReweight>(hc);
    }
    std::vector<std::vector<ColumnSpec>*> recoSpecs{&recoSpec};
    if (isMC) { recoSpecs.push_back(&recoMCSpec); recoSpecs.push_back(&truthSpec); }
    nReco = mad->GetEntries();
    nRecoPass = ProcessTree(mad, treco, recoSpecs, false, isMC, o, fileId, fileIdOut, hadronRW.get());
    nTruth = nTruthPass = 0;
    if (isMC && !o.noTruth) {
      nTruth = truth->GetEntries();
      nTruthPass = ProcessTree(truth, ttruth, {&truthSpec}, true, true, o, fileId, fileIdOut, hadronRW.get());
    }
    snprintf(fname, sizeof(fname), "%s", path.c_str());
    tmeta->Fill();
    double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "makedf: file " << fileId << " " << path << " mc=" << isMC << " POT=" << potUsed
              << " reco " << nRecoPass << "/" << nReco << " truth " << nTruthPass << "/" << nTruth
              << " in " << dt << " s" << std::endl;
    ++fileId;
  }

  fout.cd();
  treco->Write();
  if (firstMC == 1) ttruth->Write();
  tmeta->Write();
  fout.Close();
  if (nfail > 0) {
    std::cerr << "makedf: " << nfail << " input file(s) failed" << std::endl;
    return 6;
  }
  return 0;
}
