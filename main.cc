#include <iostream>
#include <vector>

#include "solver.h"

namespace {
std::vector<int32_t> coast, used;
constexpr int64_t kInf = static_cast<int32_t>(1e9);
constexpr int64_t kMaxCoast = 1001;
}  // namespace

void read_and_build(Network<int64_t>& Net) {
  size_t nTask, nVertex;
  std::cin >> nTask;
  Net = Network<int64_t>(nTask + 2, 0, 1);
  coast.assign(nTask + 2, 0);
  nVertex = nTask + 2;
  TVertex source = 0, sink = 1;

  for (TVertex vertex = 2; vertex < nVertex; ++vertex) {
    std::cin >> coast[vertex];
    Net.addEdge(source, vertex, kMaxCoast + coast[vertex]);
    Net.addEdge(vertex, sink, kMaxCoasts);
  }

  for (TVertex from = 2; from < nVertex; ++from) {
    size_t k;
    std::cin >> k;
    for (size_t j = 0; j < k; ++j) {
      TVertex to;
      std::cin >> to;
      ++to;
      Net.addEdge(from, to, kInf);
    }
  }
}

int64_t getAns(Network<int64_t>& Net, TVertex vertex) {
  int64_t ans = 0;
  ans += coast[vertex];
  used[vertex] = 1;
  for (auto it = Net.begin(vertex); it.isValid(); it.goNext()) {
    if (it.getResidualCapacity() != 0 && !used[it.getFinish()])
      ans += getAns(Net, it.getFinish());
  }
  return ans;
}

void solver(Network<int64_t>& Net) {
  MalhotraKumarMaheshwari<int64_t> MKM(Net);
  int64_t flow1 = MKM.findMaxFlow();

  Net.clearFlow();

  Preflow<int64_t> preF(Net);
  int64_t flow2 = preF.findMaxFlow();

  if (flow2 != flow1) {
    std::cout << "ERROR";
    return;
  }

  used.assign(Net.getNumVertices(), 0);
  std::cout << getAns(Net, 0) << "\n";
}

int main() {
  std::vector<int32_t> coast;
  Network<int64_t> Net;
  read_and_build(Net);
  solver(Net);
  return 0;
}
