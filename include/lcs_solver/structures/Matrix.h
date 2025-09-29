#ifndef LCS_SOLVER_STRUCTURES_MATRIX_H_
#define LCS_SOLVER_STRUCTURES_MATRIX_H_

#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <span>
#include <string>
#include <util/Logger.hpp>
#include <vector>

namespace lcs_solver::structures {

template<typename T>
class MatrixConstIter;

template<typename T>
class MatrixIter;

// === Matrix Class Definition =================================================
template<typename T>
class Matrix {
 public:
  using uint = std::size_t;

  Matrix();
  Matrix(std::initializer_list<size_t> dimensions);
  explicit Matrix(const std::vector<size_t> &dimensions, T init = T());

  void clear();
  [[nodiscard]] uint linearIndex(std::span<const uint> idx) const;
  [[nodiscard]] uint linearIndex(std::initializer_list<uint> indices) const;
  [[nodiscard]] std::vector<uint> vectorizeIndex(uint pos) const;
  [[nodiscard]] bool empty() const;
  [[nodiscard]] uint size() const;

  [[nodiscard]] MatrixConstIter<T> cbegin() const;
  [[nodiscard]] MatrixConstIter<T> cend() const;
  [[nodiscard]] MatrixConstIter<T> crbegin() const;
  [[nodiscard]] MatrixConstIter<T> crend() const;

  [[nodiscard]] MatrixIter<T> begin();
  [[nodiscard]] MatrixIter<T> end();
  [[nodiscard]] MatrixIter<T> rbegin();
  [[nodiscard]] MatrixIter<T> rend();

  T &operator[](uint index);
  const T &operator[](uint index) const;
  T &operator[](std::span<const uint> indices);
  const T &operator[](std::span<const uint> indices) const;

  T &stepDown(uint pos, uint d);
  [[nodiscard]] const T &stepDown(uint pos, uint d) const;
  T &stepDownDiag(uint pos);
  [[nodiscard]] const T &stepDownDiag(uint pos) const;

  [[nodiscard]] uint getStepDownIndex(uint pos, uint d) const;
  [[nodiscard]] uint getStepDownDiagIndex(uint pos) const;
  [[nodiscard]] uint getMax() const;
  [[nodiscard]] uint getIdxAbs(uint startIdx, uint d, uint i) const;
  [[nodiscard]] uint getIdxRel(uint startIdx, uint d, int i) const;
  [[nodiscard]] uint getIdxFst(uint idx, uint d) const;
  [[nodiscard]] uint vectorizeIndex(uint idx, uint d) const;
  [[nodiscard]] std::string DebugString() const;
  [[maybe_unused]] bool incrementIndex(uint &idx) const;
  [[maybe_unused]] bool decrementIndex(uint &idx) const;
  bool incr(std::vector<uint> &vec) const;

  void setAxisOrder(const std::vector<uint> &vec);
  [[nodiscard]] const std::vector<uint> &getAxisOrder() const;

  std::vector<size_t> dim;
  std::vector<uint> lastIndex;
  std::vector<uint> stride;///< Def: stride[i] = (i==0) ? 1 : 1*dim[K]*dim[K-1]*...*dim[K-(i+1)] with K=dim.size()-1
  uint strideDiag;
  std::vector<T> data;
  std::vector<uint> axis_order;

 private:
  static std::vector<uint> calcProduct(const std::vector<uint> &dim);///< calculates product used in constructors
  static std::vector<uint> createLastIndex(const std::vector<uint> &dim);
  [[nodiscard]] std::vector<uint> getDefaultAxisOrder() const;
  T zero = T();
  T &getZero();
  [[nodiscard]] const T &getZero() const;
};

//=== ostream operator =========================================================
template<typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &matrix);

// Declare explicit instantiations
extern template class Matrix<int>;
extern template class Matrix<util::uint>;

}// namespace lcs_solver::structures
#endif// LCS_SOLVER_STRUCTURES_MATRIX_H_