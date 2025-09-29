#ifndef LCS_SOLVER_ALGORITHMS_SOLUTIONS_EMPTYSOLUTION_H_
#define LCS_SOLVER_ALGORITHMS_SOLUTIONS_EMPTYSOLUTION_H_

#include "algorithms/BaseSolution.h"

namespace lcs_solver::algorithms::solutions {
class EmptySolution : public BaseSolution {
 public:
  EmptySolution();

  [[nodiscard]] SolutionType getType() const override;
  [[nodiscard]] bool isEqual(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessThan(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessEqualThan(const BaseSolution &rhs) const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] BaseSolution *clone() const override;
  [[nodiscard]] bool empty() const override;
};

}
#endif //LCS_SOLVER_ALGORITHMS_SOLUTIONS_EMPTYSOLUTION_H_
