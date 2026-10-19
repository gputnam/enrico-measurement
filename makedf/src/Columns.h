// Config-driven column extraction from flat MINERvA AnaTuple branches.
//
// A column file has one column per line:
//     <out_name>  <source>
// where <source> is one of
//     branch          scalar branch (or element 0 of an array)
//     branch[i]       element i of a fixed or variable length array
//     branch[:K]      first K elements -> out_name_0 .. out_name_{K-1}, padded
// Blank lines and '#' comments are ignored. Double/Float branches are written as
// Float_t, Int/Bool/Short branches as Int_t. Padding is NaN (float) / -9999 (int).
#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

class TTree;

namespace makedf {

// A bound input branch: owns the read buffer.
struct InBranch {
  std::string name;
  char type = 0;            // 'd' double, 'f' float, 'i' int, 'b' bool, 's' short, 'l' long64
  size_t elsize = 0;
  size_t capacity = 0;      // number of elements the buffer can hold
  size_t lenStatic = 1;     // static leaf length (e.g. 4 for vtx[4])
  InBranch* count = nullptr;// count branch for variable length arrays
  std::vector<char> buf;

  size_t size() const;      // current number of elements
  double get(size_t i) const;
};

// Registry of input branches for one TTree. Branches are bound lazily on request.
class BranchSet {
 public:
  explicit BranchSet(TTree* t) : fTree(t) {}
  // Returns nullptr if the branch does not exist in the tree.
  InBranch* request(const std::string& name);
  // Like request() but throws if missing.
  InBranch* require(const std::string& name);
  bool has(const std::string& name) const;
  // Enable/disable branch statuses and set addresses. Call once after all requests.
  void bind();
  TTree* tree() const { return fTree; }

 private:
  TTree* fTree;
  std::map<std::string, std::unique_ptr<InBranch>> fBranches;
};

struct OutColumn {
  std::string name;
  InBranch* src = nullptr;
  size_t index = 0;
  bool isInt = false;
  float fval = 0;
  int ival = 0;
  void fill();
};

struct ColumnSpec {
  std::string out;
  std::string branch;
  size_t index = 0;
  size_t K = 0;   // >0 means slice [:K]
};

std::vector<ColumnSpec> ReadColumnFile(const std::string& path);

// Creates branch `name` in `t` pointing at `addr`, or re-points it if it already exists.
// Creating a new branch in a tree that already has entries is an error (inconsistent inputs).
void BindOutput(TTree* t, const std::string& name, void* addr, char type);

// Binds specs to input branches in `bs` and binds output branches in `out`.
// Specs whose branch is missing are skipped (with a warning) unless `strict`.
std::vector<std::unique_ptr<OutColumn>> MakeColumns(const std::vector<ColumnSpec>& specs,
                                                     BranchSet& bs, TTree* out, bool strict);

}  // namespace makedf
