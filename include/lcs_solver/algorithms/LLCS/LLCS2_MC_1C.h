#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_1C_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_1C_H_

#include "algorithms/LLCS/LLCS2_MC_Algorithm.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/BaseSolution.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::llcs {

struct LLCS2_MC_1C final : public LLCS2_MC_Algorithm {
  static constexpr const char *name = "LLCS2_MC_1C";

  LLCS2_MC_1C(const StrPtrVector &vec, const ConstraintMap &map);

  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  [[nodiscard]] const Matrix &getMatrix() const override;
  void doPreprocessing() override;
  void reset(ResetLevel lvl) override;

 private:
  const lcs_solver::util::StringView v, w; ///< Input strings
  const uint m, n; ///< Length of input strings m = |v| <= |w| = n
  const uint l, u; ///< lower and upper bound for the length of a gap in v (or w)
  std::vector<std::vector<uint>> M; ///< M[i][j]==p means there is a gc[1:p-1] subsequence of v and w with length p
};

} // lcs_solver::algorithms::llcs
#endif /* LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_1C_H_ */