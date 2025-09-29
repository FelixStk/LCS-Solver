/*******************************************************************************
 * @file BaseSolution.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of a BaseSolution for solutions
 ******************************************************************************/

#include "algorithms/BaseSolution.h"

#include <iostream>

namespace lcs_solver::algorithms {

//bool BaseSolution::empty() const {
//  return getType() == BaseSolution::SolutionType::Empty;
//}

bool BaseSolution::operator==(const BaseSolution &rhs) const {
  if(empty() && rhs.empty())
    return true;
  return isEqual(rhs);
}

bool BaseSolution::operator!=(const BaseSolution &rhs) const {
  return !isEqual(rhs);
}

bool BaseSolution::operator<(const BaseSolution &rhs) const {
  return isLessThan(rhs);
}

bool BaseSolution::operator>(const BaseSolution &rhs) const {
  return !isLessEqualThan(rhs);
}

bool BaseSolution::operator<=(const BaseSolution &rhs) const {
  return isLessEqualThan(rhs);
}

bool BaseSolution::operator>=(const BaseSolution &rhs) const {
  return !isLessThan(rhs);
}
void BaseSolution::print() const {
  print(std::cout);
}
void BaseSolution::print(std::ostream &os) const {
  os << DebugString();
}

}// end of namespace
