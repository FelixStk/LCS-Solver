#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS_STD_FL_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS_STD_FL_H_

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/LLCS/LLCS_Algorithm.h"
#include "algorithms/BaseSolution.h"
#include "structures/Matrix.h"

namespace lcs_solver::algorithms::llcs {

struct LLCS_STD_FL : public LLCS_Algorithm {
  using uint = util::uint;
  using Matrix = structures::Matrix<util::uint>;
  static constexpr const char *name = "LLCS_STD_FL";

  explicit LLCS_STD_FL(const StrPtrVector &vec, const ConstraintMap &map);

  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  [[nodiscard]] const Matrix &getMatrix();
  void doPreprocessing() override;
  void reset(ResetLevel l) override;

 private:
  static bool matchingSymbol(const std::vector<size_t> &, const StringViewVector &);
  Matrix dp;
};

} // namespace lcs_solver::algorithms::llcs
#endif /*LCS_SOLVER_ALGORITHMS_LLCS_LLCS_STD_FL_H_*/