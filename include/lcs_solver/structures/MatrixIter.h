#ifndef LCS_SOLVER_STRUCTURES_MATRIXITERATOR_H_
#define LCS_SOLVER_STRUCTURES_MATRIXITERATOR_H_

#include <iterator>
#include <vector>
#include <string>
#include "util/CommonTypes.h"

namespace lcs_solver::structures {

template<typename T> class Matrix;

template<typename T>
class MatrixIter {
public:
  using uint = std::size_t;
  using iterator_category = std::forward_iterator_tag;
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using pointer = T*;
  using reference = T&;

  MatrixIter();
  MatrixIter(Matrix<T>& matrix, const std::vector<uint>& axis_order, bool reverse=false);
  MatrixIter(const MatrixIter&) = default;
  MatrixIter& operator=(const MatrixIter&) = default;

  reference operator*() const;
  pointer operator->() const;
  MatrixIter& operator++(); // Pre-increment operator
  MatrixIter operator++(int); // Post-increment operator
  MatrixIter& operator--(); // Pre-increment operator
  MatrixIter operator--(int); // Post-increment operator

  bool operator==(const MatrixIter& other) const;
  bool operator!=(const MatrixIter& other) const;

  std::string DebugString() const;

private:
  void incrementIndices();
  void decrementIndices();

  Matrix<T>* matrix_;
  const bool reverse_;
  const std::vector<uint> axis_order_;
  uint linear_index_;
  bool is_end_;
};

extern template class MatrixIter<int>;
extern template class MatrixIter<util::uint>;

} // namespace lcs_solver::structures
#endif //LCS_SOLVER_STRUCTURES_MATRIXITERATOR_H_
