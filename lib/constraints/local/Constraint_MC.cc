/******************************************************************************
 * @file Constraint_MC.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of local constraint: Multiple Gap Constraint
 * @details See Adamson et al. "Longest Common Subsequence with Gap Constraints"
 *  url: http://arxiv.org/abs/2304.05270
 *****************************************************************************/

#include "constraints/local/Constraint_MC.h"

#include <algorithm>
#include <set>
#include <sstream>

#include "util/Logger.hpp"

namespace {
using ::lcs_solver::util::Logger;
}

namespace lcs_solver::constraints::local {

//=== Public Functions =========================================================
/*******************************************************************************
 * Constraint_MC - Constructor
 * @details [li,ui] = gc[i] implies li <= |gap_e(w,j)| <= ui for all strings
 * @param gc GapVector
 ******************************************************************************/
Constraint_MC::Constraint_MC(const GapVector &gc)
    : BaseConstraint(ConstraintCategory::LOCAL, ConstraintType::MC),
      gap_(gc),
      t_(gc.size()),
      gap_lower_bound_(CalcLowerBound(gc)),
      gap_upper_bound_(CalcUpperBound(gc)),
      gap_max_length_(CalcMaxLength(gc)),
      gap_num_uniques_(CalcNumUniques(gc)) {}

/*******************************************************************************
 * Constraint_MC - Constructor
 * @details [li,ui] = gc[i] implies li <= |gap_e(w,j)| <= ui for all strings
 * @param gc GapVector
 * @param type ConstraintType
 ******************************************************************************/
Constraint_MC::Constraint_MC(const Constraint_MC::GapVector &gc,
                             ConstraintType type)
    : BaseConstraint(ConstraintCategory::LOCAL, type),
      gap_(gc),
      t_(gc.size()),
      gap_lower_bound_(CalcLowerBound(gc)),
      gap_upper_bound_(CalcUpperBound(gc)),
      gap_max_length_(CalcMaxLength(gc)),
      gap_num_uniques_(CalcNumUniques(gc)) {}

/*******************************************************************************
 * getName
 * @return name
 ******************************************************************************/
std::string_view Constraint_MC::GetName() const { return {kName}; }

/*******************************************************************************
 * getDescription
 * @return description
 ******************************************************************************/
std::string_view Constraint_MC::GetDescription() const {
  return {kDescription};
}

/*******************************************************************************
 * DebugString
 * @return std::string allowing a full reconstruction of the constraint
 ******************************************************************************/
std::string Constraint_MC::DebugString() const {
  std::ostringstream oss;
  oss << kName << std::endl;
  oss << gap_.size();
  for (const auto &it : gap_) {
    oss << "\n" << it.first << " " << it.second;
  }
  return oss.str();
}

/*******************************************************************************
 * isEmbeddingValid
 * @param e Embedding to test
 * @return true iff e complies with the constraint
 ******************************************************************************/
bool Constraint_MC::IsEmbeddingValid(const Embedding &e) const {
  for (uint i = 0; i < e.size() - 1; i++) {
    // Calculate the number of Symbols between e.strView[i] and e.strView[i+1]
    const uint gapLength = (e[i + 1] - 1) - (e[i] + 1) + 1;
    if (gapLength < gap_[i].first || gapLength > gap_[i].second) {
      return false;
    }
  }
  return e.isValidEmbedding();
}

/*******************************************************************************
 * isConstraintValid
 * @param spv Vector with share pointers of the problem
 * @return
 ******************************************************************************/
bool Constraint_MC::IsConstraintValid(const StrPtrVector &spv) const {
  return IsGapVectorValid(spv);
}

/*******************************************************************************
 * getGapsLowerBound
 * @return gapLowerBound
 ******************************************************************************/
Constraint_MC::uint Constraint_MC::GetGapsLowerBound() const {
  return gap_lower_bound_;
}

/*******************************************************************************
 * getGapsUpperBound
 * @return gapUpperBound
 ******************************************************************************/
Constraint_MC::uint Constraint_MC::GetGapsUpperBound() const {
  return gap_upper_bound_;
}

/*******************************************************************************
 * getGapsMaxLength
 * @return gapMaxLength
 ******************************************************************************/
Constraint_MC::uint Constraint_MC::GetGapsMaxLength() const {
  return gap_max_length_;
}

/*******************************************************************************
 * getGapsNumUniques
 * @return gapNumUniques
 ******************************************************************************/
Constraint_MC::uint Constraint_MC::GetGapsNumUniques() const {
  return gap_num_uniques_;
}

/*******************************************************************************
 * getGaps
 * @return const reference to the gc-tuple vector `gap`
 ******************************************************************************/
const std::vector<Constraint_MC::Pair> &Constraint_MC::GetGapVector() const {
  return gap_;
}

void Constraint_MC::SetGapVector(const GapVector &gc) {
  gap_ = gc;
  t_ = gap_.size();
  gap_lower_bound_ = CalcLowerBound(gap_);
  gap_upper_bound_ = CalcUpperBound(gap_);
  gap_max_length_ = CalcMaxLength(gap_);
  gap_num_uniques_ = CalcNumUniques(gap_);
}

Constraint_MC::GapVector Constraint_MC::GenRelaxedGapVector(
    const std::span<const util::String> &s) {
  uint max = 0;
  uint min = s.empty()? 0 : s.front().size();
  for (const auto& str : s) {
    max = std::max<uint>(max, str.size());
    min = std::min<uint>(min, str.size());
  }
  return GapVector{min, std::make_pair(0, max)};
}

Constraint_MC::GapVector Constraint_MC::GenRelaxedGapVector(
    const StrPtrVector &spv) {
  GapVector vec;
  const uint min = util::CalcMinStrLen(spv);
  const uint max = util::CalcMaxStrLen(spv);
  return GapVector{min, std::make_pair(0, max)};
}
// bool Constraint_MC::operator==(const BaseConstraint &other) const {
//   if (other.GetType() != kType) return false;
//   if (other.GetCategory() != kCategory) return false;
//   auto p = dynamic_cast<const Constraint_MC*>(&other);
//
//   if (gap_ != p->gap_) return false;
//   if (t_ != p->t_) return false;
//   if (gap_lower_bound_ != p->gap_lower_bound_) return false;
//   if (gap_max_length_ != p->gap_max_length_) return false;
//   if (gap_num_uniques_ != p->gap_num_uniques_) return false;
//   if (gap_upper_bound_ != p->gap_upper_bound_) return false;
//   if (gap_num_uniques_ != p->gap_num_uniques_) return false;
//
//   return true;
// }

//=== Protected Methods ========================================================
/*******************************************************************************
 * isGapVectorValid
 * @param spv StrPtrVector containing the problems strings
 * @return true iff there is no problem between spv and the `gap` vector
 ******************************************************************************/
bool Constraint_MC::IsGapVectorValid(const StrPtrVector &spv) const {
  // find the longest and smallest strings in spv
  const auto *const maxStringPtr =
      std::ranges::max_element(spv, [](const auto &a, const auto &b) {
        return a->size() < b->size();
      })->get();
  const auto *const minStringPtr =
      std::ranges::min_element(spv, [](const auto &a, const auto &b) {
        return a->size() < b->size();
      })->get();

  // Warn when shortest string prevents any gap
  if (minStringPtr->size() < 2) {
    Logger::Detail()
        << "Constraint_MC: Shortest String is very short (length <= 1)!.\n";
  }

  // Check size of gap-length constraint
  if (gap_.size() != minStringPtr->size() - 1) {
    util::Logger::Warning()
        << "Constraint_MC: The size of the tuple of gap-length constraints gc "
           "does not match the length of the shortest string. (gap.size() != "
           "|shortest string| - 1)\n";
    return false;
  }

  // Check if the tuples in the gap bounds are in an appropriate range and order
  auto isInvalid = [maxStringPtr](const Pair &p) {
    return p.first > maxStringPtr->size() || p.second > maxStringPtr->size() ||
           p.first > p.second;
  };
  if (std::ranges::any_of(gap_, isInvalid)) {
    util::Logger::Warning()
        << "Constraint_MC: Some gap tuple is invalid! (a,b) => (a <= b and a,b "
           "<= |longest string| \n";
    return false;
  }
  return true;
}

//=== Private Functions ========================================================
/*******************************************************************************
 * calcLowerBound
 * @param gc GapVector
 * @return the smallest lower bound for the length of the gaps
 ******************************************************************************/
Constraint_MC::uint Constraint_MC::CalcLowerBound(const GapVector &gc) {
  if (!gc.empty()) {
    const auto min = std::ranges::min_element(
        gc, [](const auto &a, const auto &b) { return a.first < b.first; });
    return min->first;
  }
  return 0;
}

/*******************************************************************************
 * calcLowerBound
 * @param gc GapVector
 * @return the largest upper bound for the length of the gaps
 ******************************************************************************/
Constraint_MC::uint Constraint_MC::CalcUpperBound(const GapVector &gc) {
  if (!gc.empty()) {
    const auto max = std::ranges::max_element(
        gc, [](const auto &a, const auto &b) { return a.second < b.second; });
    return max->second;
  }
  return 0;
}

/*******************************************************************************
 * calcLowerBound
 * @param gc StrPtrVector
 * @return max { ui-li | for all (li,ui) in g) }
 ******************************************************************************/
Constraint_MC::uint Constraint_MC::CalcMaxLength(const GapVector &gc) {
  if (!gc.empty()) {
    auto len = *std::ranges::max_element(gc, [](const auto &a, const auto &b) {
      return (a.second - a.first) < (b.second - b.first);
    });
    return len.second - len.first;
  }
  return 0;
}

/*******************************************************************************
 * calcNumUniques
 * @param gc StrPtrVector
 * @return number of unique elements in gc
 ******************************************************************************/
Constraint_MC::uint Constraint_MC::CalcNumUniques(const GapVector &gc) {
  std::set<Pair> const uniqueElements(gc.begin(), gc.end());
  return uniqueElements.size();
}

}  // namespace lcs_solver::constraints::local