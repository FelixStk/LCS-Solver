#ifndef LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_BR_H_
#define LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_BR_H_

#include <array>
#include <string_view>
#include "problems/BaseProblem.h"

namespace lcs_solver::problems::adamson {

/**
 * @brief LCS with the Global Constraint that the Length of the Factors containing the LCS are not greater than B.
 */
class LCS2_BR : public BaseProblem {
 public:
  static constexpr const char *kName = "LCS2_BR";
  static constexpr const char *kWhat =
      "LCS Problem in which the LCS is assumed to be in a factor of length B";

  static constexpr const std::array<ConstraintType, 3> kRequiredConstraints = {
      ConstraintType::BR,
      ConstraintType::STRINGS_2,
      ConstraintType::CONST_SIG
  };

  static constexpr const std::array<AlgorithmType, 0> kAvailableAlgorithms = {

  };

  LCS2_BR(
      std::string_view name,
      std::string_view description,
      const StrPtrVec &&spv,
      const ConstraintMap &&map
  );
  static constexpr ProblemType getType();
  [[nodiscard]] const ConstraintSpan getReqConstraints() const override;
  [[nodiscard]] const AlgorithmSpan getAlgoSpan() const override;
};

}
#endif //LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_BR_H_