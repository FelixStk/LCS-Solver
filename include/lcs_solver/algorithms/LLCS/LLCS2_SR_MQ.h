#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SR_MQ_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SR_MQ_H_

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"
#include "algorithms/BaseSolution.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::llcs {
struct LLCS2_SR_MQ final : public LLCS2_SIG_Algorithm {
  static constexpr const char name[] = "LLCS2_SR_MQ";

  LLCS2_SR_MQ(const StrPtrVector &vec,const ConstraintMap &map);

  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  void doPreprocessing() override;
  void reset(ResetLevel lvl) override;

  [[nodiscard]] const Matrix &getMatrix() const override;
  [[nodiscard]] bool isExtensible(Pair a, Pair b, uint llcsOfA) const override;
  [[nodiscard]] Window getPrevRange(const Pair &pair, uint llcs) const override;

 private:
  using Row = std::vector<uint>;
  using Matrix = std::vector<Row>;

  [[nodiscard]] bool Phi(uint i, uint j) const;

  const lcs_solver::util::StringView v, w; ///< Input strings
  const size_t m, n;  ///< Length of input strings m = |v| <= |w| = n
  Matrix M;
};
}
#endif /* LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SR_MQ_H_ */