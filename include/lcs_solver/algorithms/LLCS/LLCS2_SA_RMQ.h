#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SA_RMQ_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SA_RMQ_H_

#include <memory>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"
#include "util/CommonTypes.h"
#include "util/MetaTableCalc.h"

namespace lcs_solver::algorithms::llcs {
struct LLCS2_SA_RMQ : public LLCS2_SIG_Algorithm {
  static constexpr auto name = "LLCS2_SA_RMQ";

  LLCS2_SA_RMQ(const StrPtrVector &vec, const ConstraintMap &map);

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
  using SigmaMatrix2DMap = std::unordered_map<Symbol, Matrix2D>;
  using SigmaMatrix3DMap = std::unordered_map<Symbol, Matrix3D>;
  using SigmaSigmaGapMap = std::unordered_map<Symbol, std::unordered_map<Symbol, Pair>>;
  using WindowBounds = std::tuple<uint, uint, uint, uint>;

  static WindowBounds getWindowOffSet(uint i, uint j, Pair gap);
  void rmqUpdate(uint i, uint j);
  void fillFirstRow();
  void fillFirstColumn();
  SigmaSigmaGapMap setup();
  [[nodiscard]] uint rmqQuery(const Matrix3D &mat, uint s0, uint e0, uint s1, uint e1) const;

  const util::StringView v, w;///< Input strings
  const uint m, n;            ///< Length of input strings m = |v| <= |w| = n
  const uint maxq;            ///< Number of answers to be calculated answers: max M[b][i-2^q+1:i][j-2^q+1:j] with q in [0:maxq-1], i in [0:m] and j in [0:n]
  SigmaMatrix2DMap Mp;        ///< if v[i-1]==w[j-1]==b then M[b][i][j] := length of longest (left,right)-subsequence of v[0:i-1] and w[0:j-1], otherwise M[b][i][j] := 0
  SigmaMatrix3DMap RMQp;      ///< rmq[b][i][j][q] := max_{x in [i-2^q+1:i], y in [j-2^q+1:j]} M[b][x][y]
  Matrix2D M;                 ///< if v[i-1]==w[j-1] then M[i][j] := length of longest (left,right)-subsequence of v[0:i-1] and w[0:j-1], otherwise M[i][j] := 0

  util::Log2_table<uint> LOG2;
  util::Pow2_table<uint> POW2;
};

}// namespace lcs_solver::algorithms::llcs
#endif /* LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SA_RMQ_H_ */