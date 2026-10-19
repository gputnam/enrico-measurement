#include "Columns.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "TBranch.h"
#include "TLeaf.h"
#include "TTree.h"

namespace makedf {

size_t InBranch::size() const {
  if (!count) return lenStatic;
  long n = std::lround(count->get(0)) * (long)lenStatic;
  if (n < 0) n = 0;
  if ((size_t)n > capacity) n = capacity;
  return n;
}

double InBranch::get(size_t i) const {
  const char* p = buf.data() + i * elsize;
  switch (type) {
    case 'd': return *reinterpret_cast<const double*>(p);
    case 'f': return *reinterpret_cast<const float*>(p);
    case 'i': return *reinterpret_cast<const int*>(p);
    case 'b': return *reinterpret_cast<const bool*>(p);
    case 's': return *reinterpret_cast<const short*>(p);
    case 'l': return (double)*reinterpret_cast<const long long*>(p);
  }
  return NAN;
}

bool BranchSet::has(const std::string& name) const { return fTree->GetBranch(name.c_str()) != nullptr; }

InBranch* BranchSet::request(const std::string& name) {
  auto it = fBranches.find(name);
  if (it != fBranches.end()) return it->second.get();
  TBranch* br = fTree->GetBranch(name.c_str());
  if (!br) return nullptr;
  TLeaf* leaf = br->GetLeaf(name.c_str());
  if (!leaf) leaf = (TLeaf*)br->GetListOfLeaves()->At(0);
  if (!leaf) return nullptr;

  auto ib = std::make_unique<InBranch>();
  ib->name = name;
  std::string tn = leaf->GetTypeName();
  if (tn == "Double_t") { ib->type = 'd'; ib->elsize = 8; }
  else if (tn == "Float_t") { ib->type = 'f'; ib->elsize = 4; }
  else if (tn == "Int_t" || tn == "UInt_t") { ib->type = 'i'; ib->elsize = 4; }
  else if (tn == "Bool_t") { ib->type = 'b'; ib->elsize = 1; }
  else if (tn == "Short_t") { ib->type = 's'; ib->elsize = 2; }
  else if (tn == "Long64_t" || tn == "ULong64_t") { ib->type = 'l'; ib->elsize = 8; }
  else throw std::runtime_error("makedf: unsupported branch type " + tn + " for " + name);

  ib->lenStatic = std::max(1, leaf->GetLenStatic());
  size_t cap = ib->lenStatic;
  if (TLeaf* lc = leaf->GetLeafCount()) {
    ib->count = request(lc->GetBranch()->GetName());
    long mx = std::max(1, lc->GetMaximum());
    // generous margin: ROOT writes the whole entry into our buffer
    cap = (size_t)(2 * mx + 16) * ib->lenStatic;
  }
  ib->capacity = cap;
  ib->buf.assign(cap * ib->elsize, 0);
  InBranch* raw = ib.get();
  fBranches[name] = std::move(ib);
  return raw;
}

InBranch* BranchSet::require(const std::string& name) {
  InBranch* b = request(name);
  if (!b) throw std::runtime_error("makedf: required branch missing: " + name);
  return b;
}

void BranchSet::bind() {
  fTree->SetBranchStatus("*", 0);
  for (auto& kv : fBranches) {
    fTree->SetBranchStatus(kv.first.c_str(), 1);
    fTree->SetBranchAddress(kv.first.c_str(), kv.second->buf.data());
  }
  const char* env = std::getenv("MAKEDF_CACHE_MB");
  const long mb = env ? std::atol(env) : 100;
  fTree->SetCacheSize(mb * 1024 * 1024);
  if (mb <= 0) return;
  for (auto& kv : fBranches) fTree->AddBranchToCache(kv.first.c_str(), false);
  fTree->StopCacheLearningPhase();
}

void OutColumn::fill() {
  bool ok = src && index < src->size();
  if (isInt) ival = ok ? (int)std::lround(src->get(index)) : -9999;
  else fval = ok ? (float)src->get(index) : NAN;
}

static std::string trim(const std::string& s) {
  size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

std::vector<ColumnSpec> ReadColumnFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("makedf: cannot open column file " + path);
  std::vector<ColumnSpec> out;
  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    line = trim(line);
    if (line.empty()) continue;
    std::istringstream ss(line);
    std::string name, src;
    ss >> name >> src;
    if (src.empty()) throw std::runtime_error("makedf: bad line " + std::to_string(lineno) + " in " + path);
    ColumnSpec c;
    c.out = name;
    auto lb = src.find('[');
    if (lb == std::string::npos) {
      c.branch = src;
    } else {
      c.branch = src.substr(0, lb);
      std::string idx = src.substr(lb + 1, src.find(']') - lb - 1);
      if (!idx.empty() && idx[0] == ':') c.K = std::stoul(idx.substr(1));
      else c.index = std::stoul(idx);
    }
    out.push_back(c);
  }
  return out;
}

void BindOutput(TTree* t, const std::string& name, void* addr, char type) {
  if (TBranch* b = t->GetBranch(name.c_str())) {
    b->SetAddress(addr);
    return;
  }
  if (t->GetEntries() > 0)
    throw std::runtime_error("makedf: column " + name + " absent for earlier input files (inconsistent inputs)");
  t->Branch(name.c_str(), addr, (name + "/" + type).c_str());
}

std::vector<std::unique_ptr<OutColumn>> MakeColumns(const std::vector<ColumnSpec>& specs, BranchSet& bs,
                                                     TTree* out, bool strict) {
  std::vector<std::unique_ptr<OutColumn>> cols;
  for (const auto& s : specs) {
    InBranch* b = bs.request(s.branch);
    if (!b) {
      if (strict) throw std::runtime_error("makedf: missing branch " + s.branch);
      std::cerr << "makedf: WARNING branch not found, skipping column " << s.out << " (" << s.branch << ")\n";
      continue;
    }
    bool isInt = (b->type != 'd' && b->type != 'f');
    size_t n = s.K > 0 ? s.K : 1;
    for (size_t k = 0; k < n; ++k) {
      auto c = std::make_unique<OutColumn>();
      c->name = s.K > 0 ? s.out + "_" + std::to_string(k) : s.out;
      c->src = b;
      c->index = s.K > 0 ? k : s.index;
      c->isInt = isInt;
      if (isInt) BindOutput(out, c->name, &c->ival, 'I');
      else BindOutput(out, c->name, &c->fval, 'F');
      cols.push_back(std::move(c));
    }
  }
  return cols;
}

}  // namespace makedf
