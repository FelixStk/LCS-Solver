/*********************************************************************
 * @file Matrix.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of Matrix: N-Dim, dynamic Matrix based on std::vector
 *
 * @details Example the M(dim = {2,3,4}) is a container<T> of size 2x3x4:
 *
 *{M[{0,0,0}], M[{0,0,1}], M[{0,0,2}], M[{0,0,3}],
 * M[{0,1,0}], M[{0,1,2}], M[{0,1,2}], M[{0,1,3}],
 * M[{0,2,0}], M[{0,2,2}], M[{0,2,2}], M[{0,2,3}],
 *
 * M[{1,0,0}], M[{1,0,1}], M[{1,0,2}], M[{1,0,3}],
 * M[{1,1,0}], M[{1,1,2}], M[{1,1,2}], M[{1,1,3}],
 * M[{1,2,0}], M[{1,2,2}], M[{1,2,2}], M[{1,2,3}]}
 * =
 *{M[ 0], M[ 1], M[ 2], M[ 3],
 * M[ 4], M[ 5], M[ 6], M[ 7],
 * M[ 8], M[ 0], M[10], M[11],
 *
 * M[12], M[13], M[14], M[15],
 * M[16], M[17], M[18], M[19],
 * M[20], M[21], M[22], M[23]}
 *
 * lastIndex = {1,2,3}
 * stride = {1, 4, 4*3, 4*3*2}
 * size() = 4*3*2
 ********************************************************************/

#include "structures/Matrix.h"
#include "structures/MatrixIter.h"
#include "structures/MatrixConstIter.h"
#include "util/CommonTypes.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iostream>
#include <cassert>
#include <type_traits>

namespace lcs_solver::structures {

//=== Constructors =============================================================
template<typename T>
Matrix<T>::Matrix()
    : dim{}, lastIndex{}, stride{}, strideDiag{}, data{} {

}

template<typename T>
Matrix<T>::Matrix(const std::initializer_list<uint> dimensions)
    : Matrix<T>::Matrix(std::vector<uint>(dimensions)) {
  (void) dimensions; // d is used by the other constructor and copied into `dim`
}

template<typename T>
Matrix<T>::Matrix(const std::vector<uint> &dimensions, T init)
    : dim(dimensions),
      lastIndex(createLastIndex(dim)),
      stride(calcProduct(dim)),
      strideDiag(dimensions.empty()? 0: std::accumulate(stride.begin(),
                                 stride.end() - 1,
                                 0,
                                 [](size_t a, size_t b) { return a + b; })),
      data(dimensions.empty()? 0: std::accumulate(dimensions.begin(),
                           dimensions.end(),
                           1,
                           [](size_t a, size_t b) { return a * b; }), init),
      axis_order(getDefaultAxisOrder()){
  for (const auto &x : dim)
    if (x == 0)
      throw std::invalid_argument("Matrix dimensions must be non-zero.");
}

template<typename T>
std::vector<size_t> Matrix<T>::calcProduct(const std::vector<uint> &dim) {
  if (dim.empty()) return {};
  std::vector<size_t> result(dim.size() + 1);

  const size_t d = dim.size();
  result[0] = 1;
  for (size_t i = 0; i < dim.size(); ++i) {
    result[i + 1] = result[i] * dim[d - (i + 1)];
  }
  return result;
}

template<typename T>
std::vector<size_t> Matrix<T>::createLastIndex(const std::vector<uint> &dim) {
  if (dim.empty())
    return {};
  std::vector<uint> result(dim.size());
  std::transform(dim.begin(),
                 dim.end(),
                 result.begin(),
                 [](uint n) { return (n == 0) ? 0 : n - 1; });
  return result;
}

template<typename T>
std::vector<typename Matrix<T>::uint> Matrix<T>::getDefaultAxisOrder() const {
  std::vector<uint> defaultAxisOrder(dim.size());
  iota(defaultAxisOrder.begin(), defaultAxisOrder.end(), 0);
  return defaultAxisOrder;
}

//=== linearIndex ==============================================================
template<typename T>
size_t Matrix<T>::linearIndex(const std::span<const size_t> idx) const {
  if (idx.size() != dim.size()) {
    throw std::invalid_argument("Incorrect number of indices.");
  }

  uint index = 0;
  const uint K = dim.size();
  for (size_t i = 0; i < K; ++i) {
    index += idx[i] * stride[K - i - 1];
  }
  return index;
}

template<typename T>
size_t Matrix<T>::linearIndex(const std::initializer_list<size_t> indices) const {
  return linearIndex(std::span<const size_t>(indices.begin(), indices.size()));
}

//=== vectorizeIndex ===========================================================
template<typename T>
std::vector<size_t> Matrix<T>::vectorizeIndex(size_t pos) const {
  const size_t K = dim.size();
  auto v = std::vector<size_t>(K, 0);
  for (size_t k = K - 1, i = 0; i < K; ++i, --k) {
    const size_t r = pos % dim[k];
    v[k] = r;
    pos /= dim[k];  // removes the data of the current dimension
  }
  return v;
}

//=== Bracket Operator[] =======================================================
template<typename T>
T &Matrix<T>::operator[](size_t index) {
  return data[index];
}

template<typename T>
const T &Matrix<T>::operator[](size_t index) const {
  return data[index];
}

template<typename T>
T &Matrix<T>::operator[](std::span<const size_t> indices) {
  return data[linearIndex(indices)];
}

template<typename T>
const T &Matrix<T>::operator[](const std::span<const size_t> indices) const {
  const size_t idx = linearIndex(indices);
  return data[idx];
}

//=== empty() ==================================================================
template<typename T>
bool Matrix<T>::empty() const {
  return data.empty();
}

//=== size() ===================================================================
template<typename T>
size_t Matrix<T>::size() const {
  return data.size();
}

//=== clear() ==================================================================
template<typename T>
void Matrix<T>::clear() {
  data = std::vector<T>(size());
}

//=== getZero() ================================================================
template<typename T>
T &Matrix<T>::getZero() {
  zero = 0;
  return zero;
}

template<typename T>
const T &Matrix<T>::getZero() const {
  return zero;
}

//== stepDown ==================================================================
template<typename T>
T &Matrix<T>::stepDown(const size_t pos, const size_t d) {
  return data[getStepDownIndex(pos, d)];
}
template<typename T>
const T &Matrix<T>::stepDown(const size_t pos, const size_t d) const {
  return data[getStepDownIndex(pos, d)];
}
template<typename T>
size_t Matrix<T>::getStepDownIndex(const size_t pos, const size_t d) const {
  if (pos % stride[d + 1] < stride[d])
    return pos;
  return pos - stride[d];
}

//== stepDownDiag ==============================================================
template<typename T>
T &Matrix<T>::stepDownDiag(const size_t pos) {
  return data[getStepDownDiagIndex(pos)];
}
template<typename T>
const T &Matrix<T>::stepDownDiag(const size_t pos) const {
  return data[getStepDownDiagIndex(pos)];
}
template<typename T>
size_t Matrix<T>::getStepDownDiagIndex(const size_t pos) const {
  for (size_t d = 0; d < dim.size(); ++d)
    if (pos % stride[d + 1] < stride[d])
      return pos;
  return pos - strideDiag;
}

//=== Incr in a Mixed Radix System =============================================
template<typename T>
bool Matrix<T>::incr(std::vector<uint> &vec) const {
  if (vec.empty() || dim.empty() || vec.size() != dim.size())
    return false;

  uint k = dim.size() - 1;
  while (k < dim.size()) {
    if (vec[k] < dim[k] - 1) {
      vec[k]++;
      return true;
    } else {
      vec[k] = 0;
      if (k == 0) {
        return false;
      }
      k--;
    }
  }
  return false;
}
template<typename T>
void Matrix<T>::setAxisOrder(const std::vector<uint> &vec) {
  axis_order = vec;
}
template<typename T>
const std::vector<typename Matrix<T>::uint> &Matrix<T>::getAxisOrder() const {
  return axis_order;
}

/*******************************************************************************
 * Getter: getIdxAbs. Gets an index by moving to an absolut position along a dim
 * @tparam T Typename of data in the matrix
 * @param startIdx Index as the starting point for a move along a dim
 * @param d Identifier of the dim to move along. Zero Based. Row-major order.
 * @param i Indexes in [0:dim[d]-1] to move to
 * @return uint linearized index specified by the absolute movement
 ******************************************************************************/
template<typename T>
typename Matrix<T>::uint Matrix<T>::getIdxAbs(const uint startIdx, const uint d, const uint i) const {
  const uint K = dim.size();
  uint m = (startIdx % stride[K - d]) / stride[K - d - 1];
  assert(i < dim[d] + m && "getIdxAbs: Error i >= dim[d] - i out of bound");
  return startIdx + (i - m) * stride[K - d - 1];
}

/*******************************************************************************
 * Getter: getIdxRel. Gets an index by moving relative along a dimension
 * @tparam T Typename of data in the matrix
 * @param startIdx Index as the starting point for a move along a dim
 * @param d Identifier of the dim to move along. Zero Based. Row-major order.
 * @param i Number of Indexes to move along the dimension d
 * @return uint linearized index specified by the relative movement
 ******************************************************************************/
template<typename T>
typename Matrix<T>::uint Matrix<T>::getIdxRel(const uint startIdx, const uint d, const int i) const {
  const uint K = dim.size();
  assert((startIdx % stride[K - d]) / stride[K - d - 1] + i < dim[d]
             && "getIdxRel: i >= dim[d]");
  return startIdx + i * stride[K - d - 1];
}

/*******************************************************************************
 * Getter: First index of the d dimensional sub-matrix that contains idx
 * @tparam T Typename of data in the matrix
 * @param idx linearized index used for accessing the matrices' data
 * @param d number of the dimension (major-row order)
 * @return i = linearIndex(i1,i2, ..., id, 0, ..., 0)
 * @note Calculates return value with a preprocessed stride value
 ******************************************************************************/
template<typename T>
size_t Matrix<T>::getIdxFst(const size_t idx, const size_t d) const{
  const size_t n = dim.size()-1;
  return idx - idx % stride[n-d];
}

/*******************************************************************************
 * Getter: Given a linearized index it calculates the position along dimension d
 * @tparam T  Typename of data in the matrix
 * @param idx linearized index used for accessing the matrices' data
 * @param d number of the dimension (major-row order)
 * @return size_t position = vectorizeIndex(idx).at(d)
 * @note Calculates return value with a preprocessed stride value
 ******************************************************************************/
template<typename T>
size_t Matrix<T>::vectorizeIndex(const size_t idx, const size_t d) const {
  const size_t n = dim.size() - 1;
  assert(d <= n && "dimensionality d > matrix dimensionality");
  if (d == n)// For the last dim, return the index modulo the dim size
    return idx % dim[d];
  return (idx / stride[n - d]) % dim[d];
}

template<typename T>
typename Matrix<T>::uint Matrix<T>::getMax() const {
  size_t max = 0;
  for(auto elem : data){
    max = std::max<uint>(max, elem);
  }
  return max;
}



//=== Printing =================================================================
template<typename T>
std::string Matrix<T>::DebugString() const {
  std::stringstream ss;
  if (empty()) {
    return {"Empty Matrix"};
  }
  if (dim.size() == 1) {
    ss << "{";
    for (size_t i = 0; i < data.size(); i++) {
      ss << data[i];
      if (i != data.size() - 1)
        ss << ", ";
    }
    ss << "}";
    return ss.str();
  }
  const size_t M = dim[dim.size() - 2];
  const size_t N = dim[dim.size() - 1];
  size_t pos = 0;
  while (pos < data.size()) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        if (i == 0 && j == 0) {
          auto idx = vectorizeIndex(pos);
          ss << "M[";
          for (auto it = idx.begin(); it != idx.end() - 2; ++it)
            ss << *it << ",";
          ss << "0:" << M - 1 << ",0:" << N - 1 << "]: \n";
        }
        ss << data[pos] << " ";
        ++pos;
      }
      ss << "\n";
    }
    ss << "\n";
  }
  return ss.str();
}

template<typename T>
bool Matrix<T>::incrementIndex(uint &idx) const {
  for (int i = axis_order.size() - 1; i >= 0; --i) {
    uint axis = axis_order[i];
    if (this->vectorizeIndex(idx, axis) < this->dim[axis] - 1) {
      idx = this->getIdxRel(idx, axis, 1);
      return true;
    }
    else {
      idx = this->getIdxAbs(idx, axis, 0);
    }
  }
  return false;
}

template<typename T>
bool Matrix<T>::decrementIndex(uint &idx) const {
  for (int i = axis_order.size() - 1; i >= 0; --i) {
    uint axis = axis_order[i];
    if (this->vectorizeIndex(idx, axis) > 0) {
      idx = this->getIdxRel(idx, axis, -1);
      return true;
    }
    else {
      idx = this->getIdxAbs(idx, axis, this->dim[axis] - 1);
    }
  }
  return false;
}

//=== ostream operator =========================================================
template<typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &matrix) {
  out << matrix.DebugString() << std::endl;
  return out;
}

//=== cbegin() =================================================================
template<typename T>
MatrixConstIter<T> Matrix<T>::cbegin() const {
  return MatrixConstIter<T>(*this, axis_order, false);
}

//=== cend() ====================================================================
template<typename T>
MatrixConstIter<T> Matrix<T>::cend() const {
  return MatrixConstIter<T>();
}
//=== crbegin() ================================================================
template<typename T>
MatrixConstIter<T> Matrix<T>::crbegin() const {
  return MatrixConstIter<T>(*this, axis_order, true);
}
//=== crend() ===================================================================
template<typename T>
MatrixConstIter<T> Matrix<T>::crend() const {
  return MatrixConstIter<T>();
}

//=== begin() ==================================================================
template<typename T>
MatrixIter<T> Matrix<T>::begin() {
  return MatrixIter<T>(*this, axis_order, false);
}
//=== end() ====================================================================
template<typename T>
MatrixIter<T> Matrix<T>::end() {
  return MatrixIter<T>();
}
//=== rbegin() =================================================================
template<typename T>
MatrixIter<T> Matrix<T>::rbegin() {
  return MatrixIter<T>(*this, axis_order, true);
}
//=== rrend() ==================================================================
template<typename T>
MatrixIter<T> Matrix<T>::rend() {
  return MatrixIter<T>();
}

//=== Explicit instantiation ===================================================

template class Matrix<util::uint>;
template std::ostream &operator<<<util::uint>(std::ostream &, const Matrix<util::uint> &);
template class Matrix<int>;
template std::ostream &operator<<<int>(std::ostream &, const Matrix<int> &);


} // namespace lcs_solver::structures