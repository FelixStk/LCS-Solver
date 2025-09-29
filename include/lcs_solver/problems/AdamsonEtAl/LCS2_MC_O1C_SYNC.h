#ifndef LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_O1C_SYNC_H_
#define LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_O1C_SYNC_H_

#include <array>
#include <string_view>
#include "problems/BaseProblem.h"

namespace lcs_solver::problems::adamson {

/**
 * @brief LCS with a Constant Number of Synchronized Gap Constraints
 * @details The Synchronization property of the gap length bounds stored in gc is
 * defined in (see `Constraint_MC_O1C_SYNC`):
 * 1.) If gap[i] = (li,ui) than for ever string w and every gap_e(w,j) = w[e[j]+1 :e[j] -1] : li <= |gap_e(w,j)| <= ui.
 * 2.) The t-tuple for the gap constraints is decomposable into O(1) equivalence classes via gap[i]==gap[j] (elementwise)
 * 3.) the t-tuple is synchronized: For i,j in [0..t], if gap[i]==gap[j] and i <= j then gap[i+e] is a subset of gap[j+e]
 *     for all e >= 0 such that: i+e <= j+e <= t
 */
class LCS2_MC_O1C_SYNC : public BaseProblem {
 public:
  static constexpr const char *kName = "LCS2_MC_O1C_SYNC";
  static constexpr const char *kWhat =
      "LCS Problems in which the number of gap constraint bounds is fixed constant and the gaps are assumed to be synchronized";

  static constexpr const std::array<ConstraintType, 3> kRequiredConstraints = {
      ConstraintType::MC_O1C_SYNC,
      ConstraintType::STRINGS_2,
      ConstraintType::CONST_SIG
  };

  static constexpr const std::array<AlgorithmType, 3> kAvailableAlgorithms = {
      AlgorithmType::LLCS2_MC_O1_SYNC,
      AlgorithmType::LLCS2_MC,
      AlgorithmType::LCS2_RT
  };

  LCS2_MC_O1C_SYNC(
      std::string_view name,
      std::string_view description,
      const StrPtrVec &&spv,
      const ConstraintMap &&map
  );
  static constexpr ProblemType GetType();
  [[nodiscard]] ConstraintSpan GetReqConstraints() const override;
  [[nodiscard]] AlgorithmSpan GetAlgoSpan() const override;
};

}
#endif //LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_O1C_SYNC_H_
