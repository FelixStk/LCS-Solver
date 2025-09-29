#ifndef LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_SIGMA_H_
#define LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_SIGMA_H_

#include <array>
#include <string_view>
#include "problems/BaseProblem.h"

namespace lcs_solver::problems::adamson {

/**
 * @brief LCS with symbol dependent gap boundaries
 */
class LCS2_Sigma : public BaseProblem {
 public:
  static constexpr const char *kName = "LCS2_SIGMA";
  static constexpr const char *kWhat =
      "LCS Problems with alphabet dependent left and right gap constraints";

  static constexpr const std::array<ConstraintType, 3> kRequiredConstraints = {
      ConstraintType::SIGMA,
      ConstraintType::STRINGS_2,
      ConstraintType::CONST_SIG
  };

  static constexpr const std::array<AlgorithmType, 3> kAvailableAlgorithms = {
      AlgorithmType::LLCS2_SA_MQ,
      AlgorithmType::LLCS2_SA_RMQ,
      AlgorithmType::LCS2_RT
  };

  LCS2_Sigma(
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
#endif //LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_SIGMA_H_
