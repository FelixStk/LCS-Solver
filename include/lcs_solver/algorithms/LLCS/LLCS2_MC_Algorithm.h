#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_ALGORITHM_HPP_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_ALGORITHM_HPP_

#include <vector>
#include <utility>
#include "algorithms/AlgoType.h"
#include "algorithms/LLCS/LLCS2_Algorithm.h"

namespace lcs_solver::algorithms::llcs {
class LLCS2_MC_Algorithm : public LLCS2_Algorithm {
 public:
  using Pair = std::pair<uint, uint>;
  using GapVector = std::vector<std::pair<size_t, size_t>>;

  static std::string StringFrom(const GapVector& gc);

  ~LLCS2_MC_Algorithm() override = default;
  [[nodiscard]] const GapVector &getGaps() const;

 protected:
  const uint nMaxGaps;     ///< Maximal number of gaps (|s[0]| -1)
  const GapVector &C; ///< C[k] are the bounds for the gap after the kth lcs symbol

  LLCS2_MC_Algorithm(
      AlgoType algo,
      const StrPtrVector &vec,
      const ConstraintMap &map,
      const GapVector &gap,
      bool doTracking = false,
      bool oneBased = true
  );

  static const GapVector &getGaps(
      const ConstraintMap &map,
      ConstraintType type
  );

  [[nodiscard]] bool isExtensible(Pair a, Pair b, uint llcsOfA) const override;
  [[nodiscard]] Window getPrevRange(const Pair &pair, uint llcs) const override;

};

} //end of namespace
#endif //LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_MC_ALGORITHM_HPP_
