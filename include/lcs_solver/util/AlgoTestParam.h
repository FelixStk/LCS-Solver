#ifndef LCS_SOLVER_UTIL_ALGO_TEST_PARAM_H_
#define LCS_SOLVER_UTIL_ALGO_TEST_PARAM_H_

#include <string>
#include "algorithms/BaseSolution.h"
#include "algorithms/BaseAlgorithm.h"
#include "util/CommonTypes.h"

namespace lcs_solver::util {
struct AlgoParam {
  // Abbreviations
  using BaseSolution = lcs_solver::algorithms::BaseSolution;
  using Map = lcs_solver::algorithms::BaseAlgorithm::ConstraintMap;
  using StrPtrVector = lcs_solver::algorithms::BaseAlgorithm::StrPtrVector;
  using StringVec = std::vector<String>;

  // Members of the struct
  std::string name;
  StrPtrVector pointers;
  Map map;
  BaseSolution *sol = nullptr;

  // Methods
  AlgoParam(std::string s, const StringVec &vec, Map m, BaseSolution *p);
  AlgoParam(std::string s, StrPtrVector vec, Map m, BaseSolution *p);
  AlgoParam(const AlgoParam &rhs);
  AlgoParam &operator=(const AlgoParam &rhs);
  AlgoParam &operator=(AlgoParam &&rhs) noexcept;
  ~AlgoParam();

  static StrPtrVector convertToPtr(const StringVec& strings);

  friend std::ostream& operator<<(std::ostream& os, const AlgoParam& param) {
    //return os << param.DebugString();
    return os << param.name;
  }

  [[nodiscard]] std::string DebugString() const;
};

} // end of namespace
#endif //LCS_SOLVER_UTIL_ALGO_TEST_PARAM_H_
