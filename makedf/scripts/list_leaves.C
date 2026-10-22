// Leaf listing of an AnaTuple for makedf/scripts/source_branches.py:
//   <tree> <leaf> <type> <static length> <count leaf or ->
void list_leaves(const char* fn) {
  auto f = TFile::Open(fn);
  for (auto tn : {"MasterAnaDev", "Truth"}) {
    auto t = (TTree*)f->Get(tn);
    if (!t) continue;
    for (auto b : *t->GetListOfLeaves()) {
      auto l = (TLeaf*)b;
      printf("%s %s %s %d %s\n", tn, l->GetName(), l->GetTypeName(), l->GetLenStatic(),
             l->GetLeafCount() ? l->GetLeafCount()->GetName() : "-");
    }
  }
}
