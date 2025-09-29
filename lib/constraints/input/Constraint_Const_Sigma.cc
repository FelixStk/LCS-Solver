/*******************************************************************************
 * @file Constraint_Const_Sigma.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Constraint that determines the problems' strings alphabet size
 ******************************************************************************/

#include <memory>
#include <sstream>
#include <unordered_set>
#include "util/CommonTypes.h"
#include "constraints/input/Constraint_Const_Sigma.h"

namespace lcs_solver::constraints::input {

Constraint_Const_Sigma::Constraint_Const_Sigma(size_t alphabetSize)
    : BaseConstraint(ConstraintCategory::INPUT, ConstraintType::CONST_SIG),
      sigma_(alphabetSize) {

}

Constraint_Const_Sigma::Constraint_Const_Sigma(const BaseConstraint::StrPtrVector &strPtrVec)
    : Constraint_Const_Sigma(CountDifferentSymbolsInVector(strPtrVec)) {

}

std::string_view Constraint_Const_Sigma::GetName() const {
  return {kName};
}

std::string_view Constraint_Const_Sigma::GetDescription() const {
  return {kDescription};
}

bool Constraint_Const_Sigma::IsConstraintValid(const BaseConstraint::StrPtrVector &strPtrVec) const {
  return CountDifferentSymbolsInVector(strPtrVec) <= sigma_;
}

bool Constraint_Const_Sigma::IsEmbeddingValid(const Embedding &e) const {
  return false;
}

//*** Constraint_Const_Sigma: Private Methods for Counting Symbols *************

size_t Constraint_Const_Sigma::CountDifferentSymbolsInVector(const BaseConstraint::StrPtrVector &strPtrVec) {
  using Alphabet = std::unordered_set<
      ::lcs_solver::util::Symbol,
      ::lcs_solver::util::SymbolPerfectHash,
      ::lcs_solver::util::SymbolEqual,
      std::allocator<::lcs_solver::util::Symbol>
  >;
  Alphabet alphabet;
  for(const auto& sp : strPtrVec){
    for(const auto symbol : *sp){
      alphabet.insert(symbol);
    }
  }
  return alphabet.size();
}

std::string Constraint_Const_Sigma::DebugString() const {
  std::ostringstream oss;
  oss << kName << std::endl;
  oss << sigma_;
  return oss.str();
}

size_t Constraint_Const_Sigma::GetSigmaSize() const {
  return sigma_;
}
void Constraint_Const_Sigma::SetSigmaSize(const size_t sigma) {
  sigma_ = sigma;
}

}// end of namespace