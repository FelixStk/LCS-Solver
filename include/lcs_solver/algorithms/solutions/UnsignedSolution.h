#ifndef LCS_SOLVER_ALGORITHMS_SOLUTIONS_UNSIGNEDSOLUTION_H_
#define LCS_SOLVER_ALGORITHMS_SOLUTIONS_UNSIGNEDSOLUTION_H_
#include "algorithms/BaseSolution.h"

namespace lcs_solver::algorithms::solutions {

class UnsignedSolution final : public BaseSolution {
 public:
  explicit UnsignedSolution(uint number);
  [[nodiscard]] SolutionType getType() const override;
  [[nodiscard]] bool isEqual(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessThan(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessEqualThan(const BaseSolution &rhs) const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] BaseSolution *clone() const override;
  [[nodiscard]] bool empty() const override;
  [[nodiscard]] uint GetNumber() const;

 private:
  uint number_;
};

}  // namespace lcs_solver::algorithms::solutions
#endif  // LCS_SOLVER_ALGORITHMS_SOLUTIONS_UNSIGNEDSOLUTION_H_
