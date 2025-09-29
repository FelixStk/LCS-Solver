#ifndef LCS_SOLVER_PROBLEMS_LCS_Classic_H_
#define LCS_SOLVER_PROBLEMS_LCS_Classic_H_

#include <array>
#include <string_view>

#include "algorithms/AlgoType.h"
#include "problems/BaseProblem.h"

namespace lcs_solver::problems {


/*******************************************************************************
 * @brief Standard LCS (Longest Common Subsequence) Problem
 * This class inherits from the BaseProblem class and represents the classic
 * Longest Common Subsequence problem with two string. It provides member
 * variables for the input sequence.
 ******************************************************************************/
class LCS_Classic final : public BaseProblem {
 public:
  static constexpr auto kName = "LCS_Classic";
  static constexpr auto kWhat =
      "The well-known Longest Common Subsequence Problem in its standard form";

  static constexpr std::array<ConstraintType, 2> kRequiredConstraints = {
      ConstraintType::STRINGS_2, ConstraintType::CONST_SIG
  };

  static constexpr std::array<AlgorithmType, 4> kAvailableAlgorithms = {
      AlgorithmType::LLCS_STD_FL,
      AlgorithmType::LLCS2_STD_FL,
      AlgorithmType::LCS2_STD_S,
      AlgorithmType::LCS2_RT
  };

  static constexpr ProblemType GetType();

  LCS_Classic(
    std::span<const util::String> s,
    algorithms::AlgoType at,
    bool calc_llcs
  );

  LCS_Classic(
      std::string_view name,
      std::string_view description,
      StrPtrVec &&spv,
      ConstraintMap &&map
  );

  [[nodiscard]] ConstraintSpan GetReqConstraints() const override;
  [[nodiscard]] AlgorithmSpan GetAlgoSpan() const override;
};

}  // namespace lcs_solver::problems
#endif /*LCS_SOLVER_PROBLEMS_LCS_Classic_H_*/