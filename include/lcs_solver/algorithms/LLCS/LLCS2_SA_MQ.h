#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SA_MQ_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SA_MQ_H_

#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"
#include "structures/MaxQueueN.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::llcs {
struct LLCS2_SA_MQ final : public LLCS2_SIG_Algorithm {
  static constexpr auto name = "LLCS2_SA_MQ";

  LLCS2_SA_MQ(const StrPtrVector &vec, const ConstraintMap &map);

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
  using SigmaMatrixMap = std::unordered_map<Symbol, Matrix>;
  using SymbSymbTupleMap = std::unordered_map<Symbol, std::unordered_map<Symbol, std::pair<uint, uint>>>;
  using SymbSymbQueueMap = std::unordered_map<Symbol, std::unordered_map<Symbol, lcs_solver::structures::MaxQueue2D<uint>>>;

  static uint getMaxMatrixElement(const Matrix &matrix);
  static void doLineSetups(uint currentLine, SymbSymbQueueMap &D, SymbSymbTupleMap &gapLR);
  [[nodiscard]] SymbSymbTupleMap calcGapIntersection() const;
  SymbSymbQueueMap getMaxQueue(SymbSymbTupleMap &gapLR);
  [[nodiscard]] bool Phi(uint i, uint j) const;

  const util::StringView v, w;///< Input strings
  const uint m, n;            ///< Length of input strings m = |v| <= |w| = n
  const std::vector<Symbol> alphabet;

  SigmaMatrixMap Mp;///< if v[i-1]==w[j-1]==b then M[b][i][j] := length of
                    ///< longest (left,right)-subsequence of v[0:i-1] and
                    ///< w[0:j-1], otherwise M[b][i][j] := 0
  Matrix M;         ///< if v[i-1]==w[j-1] then M[i][j] := length of longest
                    ///< (left,right)-subsequence of v[0:i-1] and w[0:j-1], otherwise
                    ///< M[i][j] := 0
};

}// namespace lcs_solver::algorithms::llcs
#endif// LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SA_MQ_H_