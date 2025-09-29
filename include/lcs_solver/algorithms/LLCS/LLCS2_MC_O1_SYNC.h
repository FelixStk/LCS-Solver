#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_O1_SYNC_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_O1_SYNC_H_

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_MC_Algorithm.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::llcs {
struct LLCS2_MC_O1_SYNC final : public LLCS2_MC_Algorithm {
  static constexpr const char *name = "LLCS2_MC_O1_SYNC";

  LLCS2_MC_O1_SYNC(const StrPtrVector &vec, const ConstraintMap &map);

  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  [[nodiscard]] const Matrix &getMatrix() const override;
  void doPreprocessing() override;
  void reset(ResetLevel lvl) override;

 private:
  using Row = std::vector<uint>;
  using Matrix2D = std::vector<Row>;
  using Matrix3D = std::vector<Matrix2D>;

  [[nodiscard]] static uint findMaxIn3DMatrix(const Matrix3D &matrix);
  [[nodiscard]] bool Phi(uint i, uint j) const;

  const util::StringView v, w;///< Input strings
  const uint m, n;            ///< Length of input strings m = |v| <= |w| = n
  Matrix3D M;                 ///< M[r][i][j] == x means: there exists a mc subsequence in v[1:m] and w[1:n] with length x that is extendable with Cp[r]. Cp is the set of unique gaps.
  Matrix2D R;                 ///< R[i][j] == x means: there exists a mc subsequence of length x in v[1:m] and w[1:n]
};
}// namespace lcs_solver::algorithms::llcs
#endif /* LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_O1_SYNC_H_ */