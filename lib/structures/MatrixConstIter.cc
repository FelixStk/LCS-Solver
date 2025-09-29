#include "structures/MatrixConstIter.h"
#include "structures/Matrix.h"
#include "util/CommonTypes.h"
#include <sstream>
namespace lcs_solver::structures {
// Default constructor
template<typename T>
MatrixConstIter<T>::MatrixConstIter()
    : matrix_(nullptr),
      reverse_(false),
      axis_order_(),
      linear_index_(0),
      is_end_(true) {}

// Parameterized constructor
template<typename T>
MatrixConstIter<T>::MatrixConstIter(
    const Matrix<T> &matrix,
    const std::vector<uint> &axis_order,
    bool reverse)
    : matrix_(&matrix),
      reverse_(reverse),
      axis_order_(axis_order),
      linear_index_(!reverse_ ? 0 : matrix_->size() - 1),
      is_end_(false) {
  if (axis_order_.size() != matrix_->dim.size()) {
    throw std::invalid_argument("Axis order must match the number of dimensions in the matrix.");
  }
}

// Dereference operators
template<typename T>
typename MatrixConstIter<T>::reference MatrixConstIter<T>::operator*() const {
  if (!matrix_) {
    throw std::runtime_error("Dereferencing a default-constructed iterator.");
  }
  return (*matrix_)[linear_index_];
}

template<typename T>
typename MatrixConstIter<T>::pointer MatrixConstIter<T>::operator->() const {
  if (!matrix_) {
    throw std::runtime_error("Dereferencing a default-constructed iterator.");
  }
  return &(*matrix_)[linear_index_];
}

// Pre-increment operator
template<typename T>
MatrixConstIter<T> &MatrixConstIter<T>::operator++() {
  if (!matrix_) {
    throw std::runtime_error("Incrementing a default-constructed iterator.");
  }
  if (!reverse_)
    incrementIndices();
  else
    decrementIndices();
  return *this;
}

// Post-increment operator
template<typename T>
MatrixConstIter<T> MatrixConstIter<T>::operator++(int) {
  MatrixConstIter temp = *this;
  ++(*this);
  return temp;
}

// Pre-decrement operator
template<typename T>
MatrixConstIter<T> &MatrixConstIter<T>::operator--() {
  if (!matrix_) {
    throw std::runtime_error("Incrementing a default-constructed iterator.");
  }
  if (!reverse_)
    decrementIndices();
  else
    incrementIndices();
  return *this;
}

// Post-decrement operator
template<typename T>
MatrixConstIter<T> MatrixConstIter<T>::operator--(int) {
  MatrixConstIter temp = *this;
  ++(*this);
  return temp;
}

// Equality operators
template<typename T>
bool MatrixConstIter<T>::operator==(const MatrixConstIter &other) const {
  const bool sameMatrix = matrix_ == other.matrix_;
  const bool sameLinearIndex = linear_index_ == other.linear_index_;
  // const bool sameAxis_order = axis_order_ == other.axis_order_;
  const bool sameEnd = is_end_ == other.is_end_;
  if (sameEnd) {
    if (is_end_) {
      return true;
    }
    return sameMatrix && sameLinearIndex;
  }
  return false;
}

template<typename T>
bool MatrixConstIter<T>::operator!=(const MatrixConstIter &other) const {
  return !(*this == other);
}

template<typename T>
void MatrixConstIter<T>::incrementIndices() {
  if (!matrix_) {
    throw std::runtime_error("Incrementing a default-constructed iterator.");
  }

  for (int i = axis_order_.size() - 1; i >= 0; --i) {
    uint axis = axis_order_[i];
    if (matrix_->vectorizeIndex(linear_index_, axis) < matrix_->dim[axis] - 1) {
      linear_index_ = matrix_->getIdxRel(linear_index_, axis, 1);
      return;
    }
    else {
      linear_index_ = matrix_->getIdxAbs(linear_index_, axis, 0);
    }
  }
  is_end_ = true;
}

template<typename T>
void MatrixConstIter<T>::decrementIndices() {
  if (!matrix_) {
    throw std::runtime_error("Decrementing a default-constructed iterator.");
  }

  for (int i = axis_order_.size() - 1; i >= 0; --i) {
    uint axis = axis_order_[i];
    if (matrix_->vectorizeIndex(linear_index_, axis) > 0) {
      linear_index_ = matrix_->getIdxRel(linear_index_, axis, -1);
      return;
    }
    else {
      linear_index_ = matrix_->getIdxAbs(linear_index_, axis, matrix_->dim[axis] - 1);
    }
  }
  is_end_ = true;
}


template<typename T>
std::string MatrixConstIter<T>::DebugString() const {
  std::ostringstream oss;
  oss << "is_end_: " << is_end_ << " - ";
  oss << "linear_index_: " << linear_index_;
  if (matrix_) {
    oss << " (";
    for (const auto &el : matrix_->vectorizeIndex(linear_index_)) {
      oss << el << " ";
    }
    oss << ")";
  }
  return oss.str();
}

//=== Explicit instantiation ===================================================
template class MatrixConstIter<int>;
template class MatrixConstIter<util::uint>;
}// namespace lcs_solver::structures