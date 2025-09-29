#ifndef LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_SIGMA_L_H_
#define LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_SIGMA_L_H_

#include <array>
#include <string_view>
#include "problems/BaseProblem.h"

namespace lcs_solver::problems::adamson {

/**
 * @brief LCS with 2 strings and symbols defining constraints for gap to the left of them
 */
class LCS2_Sigma_L : public BaseProblem {
 public:
  static constexpr const char *kName = "LCS2_Sigma_L";
  static constexpr const char *kWhat =
      "LCS Problems with an alphabet dependent left gap constraints";

  static constexpr const std::array<ConstraintType, 3> kRequiredConstraints = {
      ConstraintType::SIGMA_L,
      ConstraintType::STRINGS_2,
      ConstraintType::CONST_SIG
  };

  static constexpr const std::array<AlgorithmType, 3> kAvailableAlgorithms = {
      AlgorithmType::LLCS2_SL_MQ,
      AlgorithmType::LLCS2_SL_RMQ,
      AlgorithmType::LCS2_RT
  };

  LCS2_Sigma_L(
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
#endif //LCS_SOLVER_PROBLEMS_ADAMSONETAL_LCS2_SIGMA_L_H_
