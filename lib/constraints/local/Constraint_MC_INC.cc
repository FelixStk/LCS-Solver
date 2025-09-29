/******************************************************************************
 * @file Constraint_MC_INC.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of local constraint: Multiple and Increasing Gap Constraint
 * @details See Adamson et al. "Longest Common Subsequence with Gap Constraints"
 *  url: http://arxiv.org/abs/2304.05270
 *****************************************************************************/
#include "constraints/local/Constraint_MC_INC.h"
#include <sstream>
#include "util/Logger.hpp"
namespace lcs_solver::constraints::local {


//== Constraint_MC_INC =================================================================================================

Constraint_MC_INC::Constraint_MC_INC(const GapVector &gc)
    : Constraint_MC(gc, ConstraintType::MC_INC) {}

std::string_view Constraint_MC_INC::GetName() const {
  return {kName};
}

std::string_view Constraint_MC_INC::GetDescription() const {
  return kDescription;
}

/*******************************************************************************
 * isConstraintValid
 * @details The increasing property between two gaps g1 = (l1,u1) & g2 =(u2,l2)
 *  means that g1 is a sub-interval of g2. Here is a visual representation:
 *  l2--l1+++++u1----u2
 * @param spv StrPtrVector to test
 * @return true if `gap` is a valid gap vector and it satisfies the property
 *  INC(`gap[i]`,`gap[i+1]`) for all i in [0:`gap.size()`-2]
 ******************************************************************************/
bool Constraint_MC_INC::IsConstraintValid(const StrPtrVector &spv) const {
  if (!IsGapVectorValid(spv))
    return false;

  // Check the increasing property of the gaps
  const auto gc = GetGapVector();
  for (uint i = 0; i < gc.size() - 1; i++) {
    const auto &[l1, u1] = gc[i];
    const auto &[l2, u2] = gc[i + 1];
    const bool isIncreasing = l2 <= l1 && u1 <= u2;    // visual representation:
    if (!isIncreasing) {
      util::Logger::Warning() << "Constraint_MC_INC: Invalid Constraint "
                              << "[" << gc[i].first << "," << gc[i].second
                              << "] is not subset of"
                              << "[" << gc[i + 1].first << "," << gc[i + 1].second
                              << "]"
                              << "(gap[" << i << "] & gap[" << i + 1
                              << "] violates increasing gap constraint property)!"
                              << "\n";
      return false;
    }
  }/*i*/
  return true;
}

bool Constraint_MC_INC::IsEmbeddingValid(const Embedding &e) const {
  return Constraint_MC::IsEmbeddingValid(e);
}

std::string Constraint_MC_INC::DebugString() const {
  std::stringstream oss;
  oss << kName << std::endl;
  const auto& gc = GetGapVector();
  oss << gc.size();
  for (const auto &it : gc) {
    oss << "\n" << it.first << " " << it.second;
  }
  return oss.str();
}

} // end of namespace