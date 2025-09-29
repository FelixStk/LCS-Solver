#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_STD_FL_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_STD_FL_H_

#include "algorithms/LLCS/LLCS2_Algorithm.h"

namespace lcs_solver::algorithms::llcs {

struct LLCS2_STD_FL final : public LLCS2_Algorithm {
  static constexpr auto name = "LLCS2_STD_FL";
  static const std::string_view description;

  explicit LLCS2_STD_FL(const StrPtrVector &vec, const ConstraintMap &map);

  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  void doPreprocessing() override;
  void reset(ResetLevel l) override;

  [[nodiscard]] const Matrix & getMatrix() const override;
  [[nodiscard]] bool isExtensible(Pair a, Pair b, uint llcsOfA) const override;
  [[nodiscard]] Window getPrevRange(const Pair &pair, uint llcs) const override;

private:
  Matrix dp;
};

} // namespace lcs_solver::algorithms::llcs
#endif /*LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_STD_FL_H_*/