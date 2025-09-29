#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SR_RMQ_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SR_RMQ_H_

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"
#include "algorithms/BaseSolution.h"
#include "util/CommonTypes.h"
#include "util/MetaTableCalc.h"

namespace lcs_solver::algorithms::llcs {
struct LLCS2_SR_RMQ : public LLCS2_SIG_Algorithm {
  static constexpr const char name[] = "LLCS2_SR_RMQ";

  LLCS2_SR_RMQ(const StrPtrVector &vec,const ConstraintMap &map);

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
  using Matrix2D = std::vector<Row>;
  using Matrix3D = std::vector<Matrix2D>;
  using Log2Tabel = ::lcs_solver::util::Log2_table<uint>;
  using Pow2Tabel = ::lcs_solver::util::Pow2_table<uint>;

  [[nodiscard]] uint rmqQuery(uint s0, uint e0, uint s1, uint e1) const;
  void rmqUpdate(uint i, uint j);
  void setFirstRowColum();

  const util::StringView v, w; ///< Input strings
  const uint m, n; ///< Length of input strings m = |v| <= |w| = n
  const uint maxq; //<< Determines the number of calculated rmq answers: max M[i:i+2^q][j:j+2^q] with q in [0:maxq-1], i in [0:m] and j in [0:n]

  Matrix2D M;     ///< M[i][j]==p, iff there's a LCS that satisfies the sigma right constraint and ends at v[i-1]==w[j-1], otherwise 0
  Matrix3D rmq;   ///< m x n x ( 1 + ceil( log m) ) matrix, rmq[i,j,q] stores the maximum of M[I][J] with I = [i-2^q+1:i], J=[j-2^q+1:j], otherwise 0
  Log2Tabel LOG2;
  Pow2Tabel POW2;
};
}
#endif /* LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SR_RMQ_H_ */