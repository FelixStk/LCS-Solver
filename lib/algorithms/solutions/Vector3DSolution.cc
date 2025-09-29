/*******************************************************************************
 * @file Vector3DSolution.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Impl. of a  BaseCollector. Stores a std::vector<solutions::Points>
 ******************************************************************************/

#include "algorithms/solutions/Vector3DSolution.h"

namespace lcs_solver::algorithms::solutions {

//== ListIterator ==============================================================
Vector3DSolution::VecIterator::VecIterator(const VecIter iter, const VecIter end)
    : mCurrent(iter), mEnd(end) {
  if (iter == end) {
    is_end = true;
    current = nullptr;
  } else {
    current = iter.operator->();
  }
}

void Vector3DSolution::VecIterator::advance() {
  if (is_end) {
    return;
  }
  ++mCurrent;
  is_end = mCurrent == mEnd;
  if (is_end) {
    current = nullptr;
    return;
  }
  const Points &temp = *mCurrent;
  current = &temp;
}

std::unique_ptr<BaseIterator> Vector3DSolution::VecIterator::clone() const {
  return std::make_unique<VecIterator>(mCurrent, mEnd);
}

std::string Vector3DSolution::VecIterator::DebugString() const {
  return BaseIterator::DebugString();
}

//== EmbeddingsList ============================================================
Vector3DSolution::Vector3DSolution() = default;

Vector3DSolution::Vector3DSolution(const StrPtrVec &spv, const Matrix3D &matrix) {
  for (const auto &mat2d : matrix) {
    mSolutions.emplace_back(spv, mat2d);
  }
}

Vector3DSolution::Vector3DSolution(const std::vector<Points> &vec){
  // If an element in vec is empty, it will be ignored. The reason for this is
  // that during BaseCollector::isEqual we compare sizes of containers.
  // Currently, an empty Vector3DSolution is always an empty `mSolutions` and
  // `empty()` is defined by integrators (on `mSolution` via an `AnyIterator` in
  // `BaseCollector`) that call the overloaded `begin()` and `end()` functions
  for (const auto &points : vec) {
    if (!points.empty()) {
      mSolutions.emplace_back(points);
    }
  }
}

AnyIterator Vector3DSolution::begin() const {
  // const auto iter = mSolutions.cbegin();
  // const auto end = mSolutions.cend();
  // auto res = ListIterator(iter, end);
  // return res;
  auto temp = AnyIterator(VecIterator(mSolutions.cbegin(), mSolutions.cend()));
  return temp;
}

AnyIterator Vector3DSolution::end() const {
  return AnyIterator(VecIterator(mSolutions.cend(), mSolutions.cend()));
}

BaseSolution *Vector3DSolution::clone() const {
  return new Vector3DSolution(mSolutions);
}

template<typename... Args>
void Vector3DSolution::emplace_back(Args &&...args) {
  mSolutions.emplace_back(Points{std::forward<Args>(args)...});
}

void Vector3DSolution::push_back(const Points &element) {
  mSolutions.push_back(element);
}

void Vector3DSolution::push_back(const Points &&element) {
  mSolutions.push_back(element);
}

void Vector3DSolution::modify(const uint solIdx, const uint lcsIdx, const uint strIdx, const uint pos, const bool reverse) {
  mSolutions[solIdx].modify(lcsIdx, strIdx, pos, reverse);
}

}// namespace lcs_solver::algorithms::solutions