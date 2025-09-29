/******************************************************************************
 * @file Constraint_MC_O1C_SYNC.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of local constraint: Const Number of synchronized gc-tuples
 * @details See Adamson et al. "Longest Common Subsequence with Gap Constraints"
 *  url: http://arxiv.org/abs/2304.05270
 *****************************************************************************/

#include "constraints/local/Constraint_MC_O1C_SYNC.h"
#include <sstream>
#include "util/Logger.hpp"

namespace lcs_solver::constraints::local {


//== Constraint_O1C_SYNC =======================================================

Constraint_MC_O1C_SYNC::Constraint_MC_O1C_SYNC(const GapVector &gc)
    : Constraint_MC_O1C(gc, ConstraintType::MC_O1C_SYNC) {
}

std::string_view Constraint_MC_O1C_SYNC::GetName() const {
  return {kName};
}

std::string_view Constraint_MC_O1C_SYNC::GetDescription() const {
  return kDescription;
}

std::string Constraint_MC_O1C_SYNC::DebugString() const {
  std::stringstream oss;
  oss << kName << std::endl;
  const auto &gc = GetGapVector();
  oss << gc.size();
  for (const auto &it : gc) {
    oss << "\n" << it.first << " " << it.second;
  }
  return oss.str();
}

bool Constraint_MC_O1C_SYNC::IsConstraintValid(const StrPtrVector &spv) const {
  if (!Constraint_MC_O1C::IsConstraintValid(spv))
    return false;

  // Checking synchronization by following the property definition O(n^3)
  const auto &gc = GetGapVector();
  for (Unsigned i = 0; i < gc.size(); i++) {
    for (Unsigned j = i + 1; j < gc.size(); j++) {
      if (gc[i].first == gc[j].first && gc[i].second == gc[j].second) {
        for (Unsigned e = 0; e + i < gc.size() && e + j < gc.size(); e++) {
          const Unsigned l1 = gc[i + e].first;
          const Unsigned u1 = gc[i + e].second;
          const Unsigned l2 = gc[j + e].first;
          const Unsigned u2 = gc[j + e].second;
          if (!(l2 <= l1 && u1 <= u2)) { // [l1:u1] \not\subseteq [l2:u2]
            util::Logger::Warning() << "Invalid O(1)C-SYNC Constraint: "
                                    << "Synchronization property does not hold!"
                                    << " (i,j,e) = " << i << j << e << "\n";
            return false;
          }
        } /*e*/
      }
    }/*j*/
  }/*i*/

  return true;
}

bool Constraint_MC_O1C_SYNC::IsEmbeddingValid(const Embedding &e) const {
  return Constraint_MC_O1C::IsEmbeddingValid(e);
}

} // end of namespace