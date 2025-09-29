#ifndef LCS_SOLVER_STRUCTURES_MATRIXCOSTITERATOR_H_
#define LCS_SOLVER_STRUCTURES_MATRIXCOSTITERATOR_H_

#include <iterator>
#include <vector>
#include <string>
#include "util/CommonTypes.h"

namespace lcs_solver::structures {

template<typename T> class Matrix;

template<typename T>
class MatrixConstIter {
public:
  using uint = std::size_t;
  using iterator_category = std::forward_iterator_tag;
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using pointer = const T*;
  using reference = const T&;

  MatrixConstIter();
  MatrixConstIter(const Matrix<T>& matrix, const std::vector<uint>& axis_order, bool reverse=false);
  MatrixConstIter(const MatrixConstIter&) = default;
  MatrixConstIter& operator=(const MatrixConstIter&) = default;

  reference operator*() const;
  pointer operator->() const;
  MatrixConstIter& operator++(); // Pre-increment operator
  MatrixConstIter operator++(int); // Post-increment operator
  MatrixConstIter& operator--(); // Pre-increment operator
  MatrixConstIter operator--(int); // Post-increment operator

  bool operator==(const MatrixConstIter& other) const;
  bool operator!=(const MatrixConstIter& other) const;

  std::string DebugString() const;

private:
  void incrementIndices();
  void decrementIndices();

  const Matrix<T>* matrix_;
  const bool reverse_;
  const std::vector<uint> axis_order_;
  uint linear_index_;
  bool is_end_;
};

extern template class MatrixConstIter<int>;
extern template class MatrixConstIter<util::uint>;

} // namespace lcs_solver::structures
#endif //LCS_SOLVER_STRUCTURES_MATRIXCOSTITERATOR_H_
