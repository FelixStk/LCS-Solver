#ifndef LCS_SOLVER_ALGORITHMS_SOLUTIONS_STRINGSOLUTION_H_
#define LCS_SOLVER_ALGORITHMS_SOLUTIONS_STRINGSOLUTION_H_

#include "algorithms/BaseSolution.h"

namespace lcs_solver::algorithms::solutions {
class StringSolution : public BaseSolution {
public:
  StringSolution();
  explicit StringSolution(std::string && s);
  explicit StringSolution(std::string_view v);
  explicit StringSolution(const util::String & s);

  [[nodiscard]] SolutionType getType() const override;
  [[nodiscard]] bool isEqual(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessThan(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessEqualThan(const BaseSolution &rhs) const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] BaseSolution *clone() const override;
  [[nodiscard]] bool empty() const override;
private:
  std::string s_;
};

}
#endif //LCS_SOLVER_ALGORITHMS_SOLUTIONS_STRINGSOLUTION_H_
