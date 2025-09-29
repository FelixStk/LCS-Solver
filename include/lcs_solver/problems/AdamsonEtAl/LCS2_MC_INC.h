#ifndef LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_MC_INC_H_
#define LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_MC_INC_H_

#include <array>
#include <string_view>
#include "problems/BaseProblem.h"

namespace lcs_solver::problems::adamson {

/**
 * @brief LCS with Multiple Gap Constraints given by m-1 increasing tuples
 */
class LCS2_MC_INC : public BaseProblem {
 public:
  static constexpr const char *kName = "LCS2_MC_INC";
  static constexpr const char *kWhat =
      "LCS in which successive gaps bounds are increasingly larger supersets of each other";

  static constexpr const std::array<ConstraintType, 3> kRequiredConstraints = {
      ConstraintType::MC_INC,
      ConstraintType::STRINGS_2,
      ConstraintType::CONST_SIG
  };

  static constexpr const std::array<AlgorithmType, 3> kAvailableAlgorithms = {
      AlgorithmType::LLCS2_MC_INC,
      AlgorithmType::LLCS2_MC,
      AlgorithmType::LCS2_RT
  };

  LCS2_MC_INC(
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
#endif //LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_MC_INC_H_
