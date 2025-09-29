/*******************************************************************************
 * @file Constraint_LCS_2.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Constraint to specify that a problem has only two Strings
 ******************************************************************************/

#include "constraints/input/Constraint_LCS_2.h"
#include <sstream>

namespace lcs_solver::constraints::input {

Constraint_LCS_2::Constraint_LCS_2()
    : Constraint_LCS_K(ConstraintType::STRINGS_2, 2) {

}

std::string_view Constraint_LCS_2::GetName() const {
  return kName;
}

std::string_view Constraint_LCS_2::GetDescription() const {
  return kDescription;
}

std::string Constraint_LCS_2::DebugString() const {
  std::ostringstream oss;
  oss << kName;
  return oss.str();
}

}// end of namespace