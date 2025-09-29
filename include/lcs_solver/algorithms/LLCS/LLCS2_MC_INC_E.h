#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_INC_E_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_INC_E_H_

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/LLCS/LLCS2_MC_Algorithm.h"
#include "algorithms/BaseSolution.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::llcs {

struct LLCS2_MC_INC_E final: public LLCS2_MC_Algorithm {
  static constexpr const char *name = "LLCS2_MC_INC_E";

  LLCS2_MC_INC_E(const StrPtrVector &vec, const ConstraintMap &map);

  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  [[nodiscard]] const Matrix &getMatrix() const override;
  void doPreprocessing() override;
  void reset(ResetLevel lvl) override;

 private:
  using Matrix = std::vector<std::vector<uint>>;

  const lcs_solver::util::StringView v, w; ///< Input strings
  const uint m, n; ///< Length of input strings m = |v| <= |w| = n
  Matrix M; ///< M[i][j]==p iff there is a gc[1:p-1] subsequence of v and w with length p
};

}  // namespace lcs_solver::algorithms::llcs

#endif //LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_INC_E_H_
