/*******************************************************************************
 * @file UnsignedSolution.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Solution Class for returning the length of the LCS (stores an uint)
 ******************************************************************************/
#include "algorithms/solutions/UnsignedSolution.h"

namespace lcs_solver::algorithms::solutions {
UnsignedSolution::UnsignedSolution(const uint number) : number_(number) {}

BaseSolution::SolutionType UnsignedSolution::getType() const {
  return BaseSolution::SolutionType::Unsigned;
}

bool UnsignedSolution::isEqual(const BaseSolution &rhs) const {
  if (getType() == rhs.getType()) {
    if (const auto ptr = dynamic_cast<const UnsignedSolution *>(&rhs)) {
      return number_ == ptr->number_;
    }
  }
  return false;
}

bool UnsignedSolution::isLessThan(const BaseSolution &rhs) const {
  if (getType() == rhs.getType()) {
    if (const auto ptr = dynamic_cast<const UnsignedSolution *>(&rhs)) {
      return number_ < ptr->number_;
    }
  }
  return getType() < rhs.getType();
}

bool UnsignedSolution::isLessEqualThan(const BaseSolution &rhs) const {
  if (getType() == rhs.getType()) {
    if (const auto ptr = dynamic_cast<const UnsignedSolution *>(&rhs)) {
      return number_ <= ptr->number_;
    }
  }
  return getType() < rhs.getType();
}

std::string UnsignedSolution::DebugString() const {
  return "Unsigned: " + std::to_string(number_);
}

BaseSolution *UnsignedSolution::clone() const {
  return new UnsignedSolution(number_);
}

bool UnsignedSolution::empty() const { return number_ == 0; }

uint UnsignedSolution::GetNumber() const {
  return number_;
}

}// end of namespace