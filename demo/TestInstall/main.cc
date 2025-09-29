#include <iostream>
#include <vector>
#include <lcs_solver/util/ParamGenerator.h>

int main() {
  using lcs_solver::util::AlgoParam;
  using lcs_solver::util::ParamGenerator;
  const unsigned int seed = 123456789;
  std::vector<AlgoParam> param = ParamGenerator::genWithoutSol(
      lcs_solver::constraints::ConstraintType::MC,
      seed,
      5,
      {// Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
       {5, 10},
       {5, 10}},
      {0, 3},// bounds for gap lengths
      {'a', 'b', 'c', 'd', 'e'});

  if (!param.empty()) {
    std::cout << "Library works correctly!" << std::endl;
    return 0;
  } else {
    std::cerr << "Something went wrong!" << std::endl;
    return 1;
  }
}
