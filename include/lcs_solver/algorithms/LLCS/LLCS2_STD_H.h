#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_STD_H_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_STD_H_H_

#include "algorithms/LLCS/LLCS2_Algorithm.h"

namespace lcs_solver::algorithms::llcs {

struct LLCS2_STD_H final : LLCS2_Algorithm {
  static constexpr auto name = "LLCS2_STD_H";
  static const std::string_view description;

  explicit LLCS2_STD_H(const StrPtrVector &vec, const ConstraintMap &map);

  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  void doPreprocessing() override;
  void reset(ResetLevel lvl) override;

  [[nodiscard]] const Matrix & getMatrix() const override;
  [[nodiscard]] bool isExtensible(Pair a, Pair b, uint llcs_of_a) const override;
  [[nodiscard]] Window getPrevRange(const Pair &pair, uint llcs) const override;

private:
  // Row top_row, bottom_rwo;
  uint llcs;
};

} // namespace lcs_solver::algorithms::llcs
#endif /*LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_STD_FL_H_*/