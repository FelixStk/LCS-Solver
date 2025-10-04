#ifndef LCS_SOLVER_ALGORITHMS_SOLUTIONS_BASECOLLECTOR_H_
#define LCS_SOLVER_ALGORITHMS_SOLUTIONS_BASECOLLECTOR_H_

#include <set>
#include "algorithms/BaseSolution.h"
#include "algorithms/solutions/BaseIterator.h"
#include "algorithms/solutions/AnyIterator.h"

namespace lcs_solver::algorithms::solutions {

class BaseCollector : public BaseSolution {

 public:
  using AnyIterator = ::lcs_solver::algorithms::solutions::AnyIterator;

  [[nodiscard]] static std::set<AnyIterator::value_type> genSet(const BaseCollector *p);

  [[nodiscard]] virtual AnyIterator begin() const = 0 ;
  [[nodiscard]] virtual AnyIterator end() const = 0;
  //[[nodiscard]] virtual BaseSolution * clone() const override = 0; // in BaseSolution

  [[nodiscard]] std::string DebugString() const final;
  [[nodiscard]] bool empty() const final;
  [[nodiscard]] SolutionType getType() const final;
  [[nodiscard]] bool isSubset(const BaseSolution &rhs) const;

  private:
  // static std::set<AnyIterator::value_type> genSet(BaseCollector* p);
  [[nodiscard]] bool isEqual(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessThan(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessEqualThan(const BaseSolution &rhs) const override;

  std::set<AnyIterator::value_type> values_;
};

}
#endif //LCS_SOLVER_ALGORITHMS_SOLUTIONS_BASECOLLECTOR_H_
