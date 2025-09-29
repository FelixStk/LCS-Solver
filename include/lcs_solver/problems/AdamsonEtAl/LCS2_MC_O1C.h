#ifndef LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_O1C_H_
#define LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_O1C_H_

#include <array>
#include <string_view>
#include "problems/BaseProblem.h"

namespace lcs_solver::problems::adamson {

/**
 * @brief LCS with 2 Strings and a Constant Number of Gap Constraints
 */
class LCS2_MC_O1C : public BaseProblem {
 public:
  static constexpr const char *kName = "LCS2_MC_O1C";
  static constexpr const char *kWhat =
      "LCS Problems in which the number of gap constraint bounds is fixed constant";

  static constexpr const std::array<ConstraintType, 3> kRequiredConstraints = {
      ConstraintType::MC_O1C,
      ConstraintType::STRINGS_2,
      ConstraintType::CONST_SIG
  };

  static constexpr const std::array<AlgorithmType, 2> kAvailableAlgorithms = {
    AlgorithmType::LLCS2_MC,
    AlgorithmType::LCS2_RT
  };

  LCS2_MC_O1C(
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
#endif //LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_O1C_H_
