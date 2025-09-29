#ifndef LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_MC_1C_H_
#define LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_MC_1C_H_

#include "problems/BaseProblem.h"

#include <array>
#include <string_view>

namespace lcs_solver::problems::adamson {

/**
 * @brief LCS with 2 strings and identical gap length bounds
 */
class LCS2_MC_1C final : public BaseProblem {
 public:
  static constexpr auto *kName = "LCS2_MC_1C";
  static constexpr auto *kWhat =
      "LCS Problem in which all gaps share the same gap length bounds (l,u)";

  static constexpr const std::array<ConstraintType, 3> kRequiredConstraints = {
      ConstraintType::MC_1C,
      ConstraintType::STRINGS_2,
      ConstraintType::CONST_SIG
  };

  static constexpr const std::array<AlgorithmType, 3> kAvailableAlgorithms = {
      AlgorithmType::LLCS2_MC_1C,
      AlgorithmType::LLCS2_MC,
      AlgorithmType::LCS2_RT
  };

  LCS2_MC_1C(
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

#endif //LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_MC_1C_H_
