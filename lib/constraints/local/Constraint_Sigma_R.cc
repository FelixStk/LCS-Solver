/******************************************************************************
 * @file Constraint_Sigma_R.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of local constraint: right symbol gap constraint map
 * @details See Adamson et al. "Longest Common Subsequence with Gap Constraints"
 *  url: http://arxiv.org/abs/2304.05270
 *****************************************************************************/

#include "constraints/local/Constraint_Sigma_R.h"
#include <set>
#include <sstream>
#include "util/Logger.hpp"

namespace lcs_solver::constraints::local {

Constraint_Sigma_R::Constraint_Sigma_R(const SigmaTupleMap &r)
    : Constraint_Sigma(SigmaTupleMap(), r, ConstraintType::SIGMA_R) {}

Constraint_Sigma_R::Constraint_Sigma_R(
    const SymbolVector &alphabet,
    const GapVector &right)
    : Constraint_Sigma(alphabet, {}, right, ConstraintType::SIGMA_R) {}

std::string_view Constraint_Sigma_R::GetName() const {
  return kName;
}

std::string_view Constraint_Sigma_R::GetDescription() const {
  return kDescription;
}

bool Constraint_Sigma_R::IsConstraintValid(const StrPtrVector &spv) const {
  return IsSigmaTupleMapValid(GetRight(), spv);
}

bool Constraint_Sigma_R::IsEmbeddingValid(const Embedding &e) const {
  return Constraint_Sigma::IsEmbeddingValid(e);
}

std::string Constraint_Sigma_R::DebugString() const {
  std::ostringstream oss;
  oss << kName << std::endl;
  oss << GetRight().size();
  for (const auto &[key, value] : GetRight()) {
    oss << "\n" << util::to_string(key) << " " << value.first << " "<< value.second;
  }
  return oss.str();
}

Constraint_Sigma_R Constraint_Sigma_R::CreateRelaxed(const StrPtrVector &strPtrVec) {
  auto r = SigmaTupleMap();
  std::set<Symbol> alphabet;
  Unsigned maxStrLength = 0;  // Saves longest string length seen

  for (const auto &p : strPtrVec) {
    for (const auto &symbol : *p) {
      alphabet.insert(symbol);
    }
    if (maxStrLength < p->size()) {
      maxStrLength = p->size();
    }
  }

  for (const auto &symbol : alphabet) {
    r.insert({symbol, std::make_pair(0, maxStrLength)});
  }
  return Constraint_Sigma_R(r);
}

} // end of namespace