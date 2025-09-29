#ifndef LCS_SOLVER_INCLUDE_ALGORITHMS_LLCS_LLCS2_ALGORITHM_H_
#define LCS_SOLVER_INCLUDE_ALGORITHMS_LLCS_LLCS2_ALGORITHM_H_

#include <string>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"

namespace lcs_solver::algorithms::llcs {

class LLCS2_Algorithm : public BaseAlgorithm {
 public:
  using uint = util::uint;
  using Row = std::vector<uint>;
  using Matrix = std::vector<Row>;
  using Pair = std::pair<uint, uint>;
  using KeyPointMatrix = std::vector<std::vector<Pair>>;
  using Window = std::pair<Pair, Pair>;

  ~LLCS2_Algorithm() override = default;
  LLCS2_Algorithm(AlgoType algo, const StrPtrVector &vec, const ConstraintMap &map, bool doTracking = false, bool oneBased = true);

  [[nodiscard]] virtual const Matrix &getMatrix() const = 0;
  [[nodiscard]] virtual bool isExtensible(Pair a, Pair b, uint llcsOfA) const = 0;
  [[nodiscard]] virtual Window getPrevRange(const Pair &pair, uint llcs) const = 0;

  void reset(ResetLevel lvl) override;

  [[nodiscard]] bool isMatched(const Pair &pair) const;
  [[nodiscard]] bool isMatched(const Pair &pair, bool oneBased) const;
  [[nodiscard]] std::vector<Pair> genAllKeyPairs(bool oneBased) const;
  [[nodiscard]] std::vector<Pair> genAllKeyPairsNaive(bool oneBased) const;

  [[nodiscard]] static uint lenKPtM(const KeyPointMatrix &kpm);
  [[nodiscard]] static uint max(const Matrix &matrix);

  void track(const Pair &pair, uint llcs);
  void track(const Pair &&pair, uint llcs);
  void setTracking(bool tracking);///< Vector of zero-based matching Indices
  [[nodiscard]] bool getTrackingFlag() const;
  [[nodiscard]] const KeyPointMatrix &getKeyPairs() const;

 protected:
  [[nodiscard]] static std::string toString(
      const Matrix &m,
      StringViewVector s,
      bool trim = true,
      const std::string &name = "",
      bool printStrings = true);

  [[nodiscard]] static std::string toString(const KeyPointMatrix &kpm, bool trackKeyPairs);

  const bool areIndicesOneBased;
  bool trackKeyPairs;
  KeyPointMatrix keyPairs;
};

}// namespace lcs_solver::algorithms::llcs

#endif//LCS_SOLVER_INCLUDE_ALGORITHMS_LLCS_LLCS2_ALGORITHM_H_