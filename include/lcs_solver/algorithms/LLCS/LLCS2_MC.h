#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_H_

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_MC_Algorithm.h"

namespace lcs_solver::algorithms::llcs {
struct LLCS2_MC final : public LLCS2_MC_Algorithm {
  using Matrices = std::vector<Matrix>;
  static constexpr auto name = "LLCS2_MC";

  LLCS2_MC(const StrPtrVector &vec, const ConstraintMap &map);

  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  [[nodiscard]] const Matrix &getMatrix() const override;
  void doPreprocessing() override;
  void reset(ResetLevel lvl) override;

  [[nodiscard]] const Matrices &getMatrices() const;

private:
  uint k, m, n;// k is the placeholder for the llcs, m = |s[0]|, n = |s[1]|
  Matrices M;  ///< dp[i][j]==1 iff there is a gc[1:k-1] subsequence s of v and w
  Matrix A;    ///< Stores the sum of d consecutive entries M_{p-1}[i][j-d+1]...M_{p-1}[i][j]
  Matrix B;    ///< Stores the sum of all entries M_{p-1}[i'][j'] with 0 <= (i - i') <d and 0 <= (j - j') <d
  void doSetupDynProgramming();
  void updateA(uint p);
  void updateB(uint p);
  bool updateM(uint p);
};
}// namespace lcs_solver::algorithms::llcs
#endif /* LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_H_ */