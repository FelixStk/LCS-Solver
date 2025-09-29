/******************************************************************************
 * @file Constraint_BR.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of global constraints (constructors, isEmbeddingValid)
 * @details See "Longest Common Subsequence with Gap Constraints" for LCS-BR
 *  http://arxiv.org/abs/2304.05270
 *****************************************************************************/

#include "constraints/global/Constraint_BR.h"
#include <algorithm>
#include <sstream>
#include "util/Logger.hpp"

namespace lcs_solver::constraints::global {

Constraint_BR::Constraint_BR(Unsigned bWindowSize)
    : BaseConstraint(ConstraintCategory::GLOBAL, ConstraintType::BR),
      b_(bWindowSize) {

}

std::string_view Constraint_BR::GetName() const {
  return {kName};
}

std::string_view Constraint_BR::GetDescription() const {
  return {kDescription};
}

bool Constraint_BR::IsConstraintValid(const BaseConstraint::StrPtrVector &strPtrVec) const {
  auto isStrSizeGreaterThanB = [this](const auto &strPtr) { return strPtr->size() > b_; };
  if (std::ranges::any_of(strPtrVec, isStrSizeGreaterThanB)) {
    using ::lcs_solver::util::Logger;
    Logger::Warning() << kName
                      << " constraint: B is larger than the some input string!"
                      << std::endl;
    return false;
  }
  return true;
}

bool Constraint_BR::IsEmbeddingValid(const Embedding &e) const {
  using lcs_solver::util::Logger;
  const auto minSizeOfFactor = e[e.size() - 1] - e[0] + 1;
  if (minSizeOfFactor > b_) {
    Logger::Warning() << "Constraint_BR: Invalid, first symbol is at position"
                      << e[0]
                      << "last symbol of the subsequence is at position: "
                      << e[e.size() - 1] << "B is " << b_ << std::endl;
    return false;
  }
  return true;
}

std::string Constraint_BR::DebugString() const {
  std::ostringstream oss;
  oss << kName << std::endl;
  oss << "BR == " << b_ << std::endl;
  return oss.str();
}

util::uint Constraint_BR::GetB() const {
  return b_;
}

void Constraint_BR::SetB(Unsigned b) {
  b_ = b;
}

//=== Dev Section ======================================================================================================
// Old Code for checking if a sub is a valid BR-subsequence of str
// bool LCS_BR::isSubsequence(std::string sub, std::string str){
//     // If str has no factor of length B return false
//     if(b > (int) str.size()){
//         return false;
//     }
//     // Use a sliding window technique to examine every possible substring
//     bool res = false;
//     for(int i = 0; i <= (int)str.size() - b; ++i) {
//         std::string substring = str.substr(i, b);
//         // Check Subsequence
//         int subIndex = 0;
//         int strIndex = 0;
//         while(subIndex < (int)sub.size() && strIndex < (int)str.size()) {
//             if(sub[subIndex] == str[strIndex]) {
//                 ++subIndex;
//             }
//             ++strIndex;
//         }
//         res = res || subIndex == (int)sub.size();
//     }
//     return res;
// }

}
