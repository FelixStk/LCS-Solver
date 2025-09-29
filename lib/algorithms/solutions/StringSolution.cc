/*******************************************************************************
 * @file StringSolution.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Implementation of a StringSolution type
 ******************************************************************************/

#include "algorithms/solutions/StringSolution.h"

namespace lcs_solver::algorithms::solutions {

StringSolution::StringSolution() = default;
StringSolution::StringSolution(std::string &&s) : s_(std::move(s)){}
StringSolution::StringSolution(std::string_view v) : s_(v.begin(), v.end()) {}
StringSolution::StringSolution(const util::String &s) : s_(util::to_string(s)) {}

BaseSolution::SolutionType StringSolution::getType() const {
  return BaseSolution::SolutionType::UFT8String;
}

bool StringSolution::isEqual(const BaseSolution &rhs) const {
  if (getType() == rhs.getType()) {
    if (const auto p = dynamic_cast<const StringSolution *>(&rhs)) {
      return s_ == p->s_;
    }
  }
  return false;
}

bool StringSolution::isLessThan(const BaseSolution &rhs) const {
  if (getType() == rhs.getType()) {
    if (const auto p = dynamic_cast<const StringSolution *>(&rhs)) {
      return s_ < p->s_;
    }
  }
  return getType() < rhs.getType();
}

bool StringSolution::isLessEqualThan(const BaseSolution &rhs) const {
  if (getType() == rhs.getType()) {
    if (const auto p = dynamic_cast<const StringSolution *>(&rhs)) {
      return s_ <= p->s_;
    }
  }
  return getType() < rhs.getType();
}

std::string StringSolution::DebugString() const {
  return util::to_string(s_);
}

BaseSolution *StringSolution::clone() const {
  return new StringSolution(s_);
}

bool StringSolution::empty() const {
  return s_.empty();
}
}// namespace lcs_solver::algorithms::solutions