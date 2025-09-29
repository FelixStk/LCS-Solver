/******************************************************************************
* @file Constraint_LLCS_Known.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Constraint for when the Length of the LCS is known
 *****************************************************************************/

#include "constraints/global/Constraint_LLCS_Known.h"
#include <algorithm>
#include <sstream>
#include "util/Logger.hpp"

namespace lcs_solver::constraints::global {

//=== Constraint LLCS Known ============================================================================================
Constraint_LLCS_Known::Constraint_LLCS_Known(Unsigned lengthOfLongestSubsequence)
    : BaseConstraint(ConstraintCategory::GLOBAL, ConstraintType::LLCS_KNOWN),
      llcs_(lengthOfLongestSubsequence) {

}

std::string_view Constraint_LLCS_Known::GetName() const {
  return {kName};
}

std::string_view Constraint_LLCS_Known::GetDescription() const {
  return kDescription;
}

bool Constraint_LLCS_Known::IsConstraintValid(const BaseConstraint::StrPtrVector &strPtrVec) const {
  using lcs_solver::util::Logger;
  auto isSizeGreaterThanLLCS =
      [this](const auto &strPtr) { return strPtr->size() > llcs_; };
  if (std::ranges::any_of(strPtrVec, isSizeGreaterThanLLCS)) {
    Logger::Warning() << kName
                      << " constraint: LLCS is larger than the some string!"
                      << std::endl;
    return false;
  }
  return true;
}

bool Constraint_LLCS_Known::IsEmbeddingValid(const Embedding &e) const {
  return e.size() == llcs_;
}

std::string Constraint_LLCS_Known::DebugString() const {
  std::ostringstream oss;
  oss << kName << std::endl;
  oss << "llcs == " << GetLlcs() << std::endl;
  return oss.str();
}

size_t Constraint_LLCS_Known::GetLlcs() const {
  return llcs_;
}
void Constraint_LLCS_Known::SetLlcs(const Unsigned llcs) {
  llcs_ = llcs;
}

}