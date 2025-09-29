/******************************************************************************
 * @file Constraint_MC_1C.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of Constraint_MC_1C: All gap-lengths share the same (l,u) bound
 * @details See Adamson et al. "Longest Common Subsequence with Gap Constraints"
 *  url: http://arxiv.org/abs/2304.05270
 *****************************************************************************/
#include "constraints/local/Constraint_MC_1C.h"

#include <algorithm>
#include <sstream>

#include "util/Logger.hpp"

namespace {
using ::lcs_solver::util::Logger;
}

namespace lcs_solver::constraints::local {

//=== Public Methods ===========================================================
Constraint_MC_1C::Constraint_MC_1C(const GapVector &gapLengthVector)
    : Constraint_MC(gapLengthVector, ConstraintType::MC_1C) {

}

Constraint_MC_1C::Constraint_MC_1C(const StrPtrVector &spv, uint lower, uint upper)
    : Constraint_MC(CreateGapVector1C(GetMinStrLen(spv) > 0 ? GetMinStrLen(spv) - 1 : 0,
                                      lower,
                                      upper),
                    ConstraintType::MC_1C) {

}

Constraint_MC_1C::Constraint_MC_1C(
    uint gap_vector_length,
    uint lower,
    uint upper
) : Constraint_MC(CreateGapVector1C(gap_vector_length, lower, upper), ConstraintType::MC_1C) {

}

Constraint_MC::GapVector Constraint_MC_1C::CreateGapVector1C(
    uint size,
    uint l,
    uint u
) {
  return std::vector<Pair>(size, {l, u});
}

std::string_view Constraint_MC_1C::GetName() const {
  return {kName};
}

std::string_view Constraint_MC_1C::GetDescription() const {
  return kDescription;
}

std::string Constraint_MC_1C::DebugString() const {
  std::stringstream oss;
  oss << kName << "\n";
  oss << GetLower() << " " << GetUpper();
  return oss.str();
}

bool Constraint_MC_1C::IsConstraintValid(const StrPtrVector &spv) const {
  if (!IsGapVectorValid(spv))
    return false;

  auto isConstant = [this](const Pair &tuple) {
    return tuple.first == GetLower() || tuple.second == GetUpper();
  };
  if (!std::ranges::all_of(GetGapVector(), isConstant)) {
    util::Logger::Warning() << "Constraint_MC_1C: t-tuple for the gap constraint is not constant!\n";
    return false;
  }
  return true;
}

bool Constraint_MC_1C::IsEmbeddingValid(const Embedding &e) const {
  return Constraint_MC::IsEmbeddingValid(e);
}

Constraint_MC::uint Constraint_MC_1C::GetMinStrLen(const StrPtrVector &spv) {
  auto comp = [](const auto &a, const auto &b) { return a->size() < b->size(); };
  auto shortestString = std::min_element(spv.begin(), spv.end(), comp);
  return (*shortestString)->size();
}

Constraint_MC::uint Constraint_MC_1C::GetLower() const {
  return GetGapsLowerBound(); // from Constraint_MC is equal to min gap[i].first for all i
}

Constraint_MC::uint Constraint_MC_1C::GetUpper() const {
  return GetGapsUpperBound();// from Constraint_MC is equal to max gap[i].second for all i
}

} // end of namespace