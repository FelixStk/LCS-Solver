#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SIG_ALGORITHM_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SIG_ALGORITHM_H_

#include <unordered_map>
#include <utility>
#include <vector>

#include "algorithms/LLCS/LLCS2_Algorithm.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::llcs {
class LLCS2_SIG_Algorithm : public LLCS2_Algorithm {
 public:
  using Pair = std::pair<uint, uint>;
  using Gap = std::pair<size_t, size_t>;
  using SymbolEqual = ::lcs_solver::util::SymbolEqual;
  using SymbolPerfectHash = ::lcs_solver::util::SymbolPerfectHash;
  using SigmaTupleMap = std::unordered_map<Symbol, Gap, SymbolPerfectHash, SymbolEqual>;
  using SymbSymbTupleMap = std::unordered_map<Symbol, SigmaTupleMap, SymbolPerfectHash, SymbolEqual>;

  ~LLCS2_SIG_Algorithm() override = default;

 protected:
  LLCS2_SIG_Algorithm(
      AlgoType algo,
      const StrPtrVector &vec,
      const ConstraintMap &map,
      const SigmaTupleMap &left,
      const SigmaTupleMap &right
  );

  static const SigmaTupleMap &getSigLMap(
      const ConstraintMap &map,
      ConstraintType type
  );
  static const SigmaTupleMap &getSigRMap(
      const ConstraintMap &map,
      ConstraintType type
  );
  static std::vector<Symbol> getAlphabet(const SigmaTupleMap &left);
  static SymbSymbTupleMap calcGapIntersection(
      const SigmaTupleMap &left,
      const SigmaTupleMap &right
  );

  [[nodiscard]] bool isExtensibleHelper(Pair a, Pair b, bool ignoreRight, bool ignoreLeft) const;

  static Window getPrefRange(const Pair &pair, Symbol c, const SigmaTupleMap &right);

 public:
  const SigmaTupleMap &left; ///< left[b] is gc tuple such that b is left to the gap and part of the lcs
  const SigmaTupleMap &right; ///< right[a] is gc tuple such that a is right to the gap and part of the lcs
  const uint sig; ///< number of different symbols found in the strings

};

} // end of namespace
#endif //LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SIG_ALGORITHM_H_