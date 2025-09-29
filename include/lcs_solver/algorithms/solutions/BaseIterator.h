#ifndef LCS_SOLVER_ALGORITHMS_SOLUTIONS_BASEITERATOR_H_
#define LCS_SOLVER_ALGORITHMS_SOLUTIONS_BASEITERATOR_H_

#include <memory>

#include "algorithms/solutions/Points.h"

namespace lcs_solver::algorithms::solutions {

class BaseIterator {
public:
  using iterator_category = std::input_iterator_tag;
  using value_type = ::lcs_solver::algorithms::solutions::Points;
  using difference_type = std::ptrdiff_t;
  using pointer = const value_type*;
  using reference = const value_type &;

  virtual ~BaseIterator() = default;
  virtual void advance() = 0;
  [[nodiscard]] virtual std::unique_ptr<BaseIterator> clone() const = 0;
  [[nodiscard]] virtual std::string DebugString() const;

  reference operator*() const;
  pointer operator->() const ;
  BaseIterator& operator++();
  bool operator==(const BaseIterator& other) const;
  bool operator!=(const BaseIterator& other) const;

protected:
  pointer current = nullptr;
  bool is_end = false;
};

}// namespace lcs_solver::algorithms::solutions
#endif//LCS_SOLVER_ALGORITHMS_SOLUTIONS_BASEITERATOR_H_
