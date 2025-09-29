#ifndef LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_MC_H_
#define LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_MC_H_

#include <array>
#include <string_view>
#include "problems/BaseProblem.h"

namespace lcs_solver::problems::adamson {

/**
 * @brief LCS with 2 Strings and Multiple Gap Constraints given by m-1 tuples
 */
class LCS2_MC : public BaseProblem {
 public:
  static constexpr const char *kName = "LCS2_MC";
  static constexpr const char *kWhat =
      "LCS with alphabet independent Gap-Constraint";

  static constexpr const std::array<ConstraintType, 3> kRequiredConstraints = {
      ConstraintType::MC,
      ConstraintType::STRINGS_2,
      ConstraintType::CONST_SIG
  };

  static constexpr const std::array<AlgorithmType, 2> kAvailableAlgorithms = {
      AlgorithmType::LLCS2_MC,
      AlgorithmType::LCS2_RT
  };

  LCS2_MC(
      std::string_view name,
      std::string_view description,
      const StrPtrVec &&spv,
      const ConstraintMap &&map
  );

  static constexpr ProblemType getType();

  [[nodiscard]] ConstraintSpan GetReqConstraints() const override;
  [[nodiscard]] AlgorithmSpan GetAlgoSpan() const override;
};

}

#endif //LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_MC_H_
