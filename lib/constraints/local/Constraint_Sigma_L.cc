/******************************************************************************
 * @file Constraint_Sigma_L.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of local constraint: left symbol gap constraint map
 * @details See Adamson et al. "Longest Common Subsequence with Gap Constraints"
 *  url: http://arxiv.org/abs/2304.05270
 *****************************************************************************/

#include "constraints/local/Constraint_Sigma_L.h"
#include <algorithm>
#include <set>
#include <sstream>
#include <utility>
#include <vector>
#include "util/Logger.hpp" // already in Constraint_Sigma header file

namespace lcs_solver::constraints::local {

Constraint_Sigma_L::Constraint_Sigma_L(const SigmaTupleMap &left)
    : Constraint_Sigma(left, SigmaTupleMap(), ConstraintType::SIGMA_L) {}

Constraint_Sigma_L::Constraint_Sigma_L(
    const SymbolVector &symbolVec,
    const GapVector &left)
    : Constraint_Sigma(symbolVec, left, {}, ConstraintType::SIGMA_L) {}

std::string_view Constraint_Sigma_L::GetName() const {
  return kName;
}

std::string Constraint_Sigma_L::DebugString() const {
  std::stringstream oss;
  oss << kName << std::endl;
  oss << GetLeft().size();
  for (const auto &[key, value] : GetLeft()) {
    oss << "\n" << util::to_string(key) << " " << value.first << " " << value.second;
  }
  return oss.str();
}

std::string_view Constraint_Sigma_L::GetDescription() const {
  return kDescription;
}

bool Constraint_Sigma_L::IsConstraintValid(const StrPtrVector &spv) const {
  return IsSigmaTupleMapValid(GetLeft(), spv);
}

bool Constraint_Sigma_L::IsEmbeddingValid(const Embedding &e) const {
  return Constraint_Sigma::IsEmbeddingValid(e);
}

Constraint_Sigma_L Constraint_Sigma_L::CreateRelaxed(const StrPtrVector &spv) {
  auto l = SigmaTupleMap();
  std::set<Symbol> alphabet;
  uint maxStrLength = 0;  // Saves longest string length seen

  for (const auto &p : spv) {
    for (const auto &symbol : *p) {
      alphabet.insert(symbol);
    }
    maxStrLength = std::max<util::uint>(maxStrLength, p->size());
  }

  for (const auto &symbol : alphabet) {
    l.insert({symbol, std::make_pair(0, maxStrLength)});
  }
  return Constraint_Sigma_L(l);
}

} // end of namespace