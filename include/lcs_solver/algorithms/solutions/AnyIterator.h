#ifndef LCS_SOLVER_ALGORITHMS_SOLUTIONS_ITERATOR_H_
#define LCS_SOLVER_ALGORITHMS_SOLUTIONS_ITERATOR_H_

#include "algorithms/solutions/BaseIterator.h"

// https://stackoverflow.com/questions/13670671/abstract-iterator-for-underlying-collections

namespace lcs_solver::algorithms::solutions {

template <typename T>
concept DerivedIterator = std::is_base_of_v<BaseIterator, T>;

class AnyIterator {

  std::unique_ptr<BaseIterator> impl;

 public:
  using iterator_category = BaseIterator::iterator_category;
  using value_type = BaseIterator::value_type;
  using difference_type = BaseIterator::difference_type;
  using pointer = BaseIterator::pointer;
  using reference = BaseIterator::reference;

  template <DerivedIterator T>
  explicit AnyIterator(T && x)
  : impl(new std::decay_t<T>(std::forward<T>(x)))
  {
  }

  AnyIterator(AnyIterator && rhs) = default;

  AnyIterator(AnyIterator const & rhs) : impl(rhs.impl->clone()) { }

  AnyIterator& operator++() {
    impl->advance();
    return *this;
  }

  BaseIterator::reference operator*() { return impl->operator*(); }
  BaseIterator::pointer operator->() { return impl->operator->(); }
  bool operator==(const AnyIterator & rhs) const {
    return *impl == *rhs.impl;
  }
  bool operator!=(const AnyIterator & rhs) const {
    bool res = *impl != *rhs.impl;
    return res;
  }
};

} // lcs_solver::algorithms::solutions
#endif //LCS_SOLVER_ALGORITHMS_SOLUTIONS_ITERATOR_H_
