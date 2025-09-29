/******************************************************************************
 * @file Constraint_MC_O1C.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of local constraint: Const number for the gap bound tuples
 * @details See Adamson et al. "Longest Common Subsequence with Gap Constraints"
 *  url: http://arxiv.org/abs/2304.05270
 *****************************************************************************/

#include "constraints/local/Constraint_MC_O1C.h"
#include <sstream>
#include "util/Logger.hpp"

namespace lcs_solver::constraints::local {

Constraint_MC_O1C::Constraint_MC_O1C(const GapVector &gc)
    : Constraint_MC(gc, ConstraintType::MC_O1C) {
}

std::string_view Constraint_MC_O1C::GetName() const {
  return kName;
}

std::string_view Constraint_MC_O1C::GetDescription() const {
  return kDescription;
}

std::string Constraint_MC_O1C::DebugString() const {
  std::stringstream oss;
  oss << kName << std::endl;
  const auto &gc = GetGapVector();
  oss << gc.size();
  for (const auto &it : gc) {
    oss << "\n" << it.first << " " << it.second;
  }
  return oss.str();
}

bool Constraint_MC_O1C::IsConstraintValid(const StrPtrVector &spv) const {
  return IsGapVectorValid(spv);
}

bool Constraint_MC_O1C::IsEmbeddingValid(const Embedding &e) const {
  return Constraint_MC::IsEmbeddingValid(e);
}

//=== Protected Methods ========================================================
Constraint_MC_O1C::Constraint_MC_O1C(const GapVector &gc, ConstraintType type)
    : Constraint_MC(gc, type) {

}

} // end of namespace