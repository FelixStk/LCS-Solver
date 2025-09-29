/*******************************************************************************
 * @file Constraint_LCS_K.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of constraint for the problem having K strings
 ******************************************************************************/

#include "constraints/input/Constraint_LCS_K.h"
#include <sstream>

namespace lcs_solver::constraints::input {

Constraint_LCS_K::Constraint_LCS_K(const Unsigned k)
    : BaseConstraint(ConstraintCategory::INPUT, ConstraintType::STRINGS_K),
      k_(k) {}

std::string_view Constraint_LCS_K::GetName() const {
  return {kName};
}

std::string_view Constraint_LCS_K::GetDescription() const {
  return {kDescription};
}

bool Constraint_LCS_K::IsEmbeddingValid(const Embedding &e) const {
  return true;
}

bool Constraint_LCS_K::IsConstraintValid(const StrPtrVector &strPtrVec) const {
  return strPtrVec.size() == k_;
}

std::string Constraint_LCS_K::DebugString() const {
  std::ostringstream oss;
  oss << kName << std::endl;
  oss << k_;
  return oss.str();
}

void Constraint_LCS_K::SetK(const Unsigned k) {
  k_ = k;
}

Constraint_LCS_K::Constraint_LCS_K(ConstraintType t, Unsigned k)
    : BaseConstraint(ConstraintCategory::INPUT, t),
      k_(k) {

}
BaseConstraint::Unsigned Constraint_LCS_K::GetK() const {
  return k_;
}

}// end of namespace