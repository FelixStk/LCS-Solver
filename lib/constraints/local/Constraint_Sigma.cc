/******************************************************************************
 * @file Constraint_Sigma.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. of local constraint: left and right symbol gap constraint map
 * @details See Adamson et al. "Longest Common Subsequence with Gap Constraints"
 *  url: http://arxiv.org/abs/2304.05270
 *****************************************************************************/

#include "constraints/local/Constraint_Sigma.h"

#include <algorithm>
#include <set>
#include <sstream>
#include <utility>

#include "util/Logger.hpp"

namespace lcs_solver::constraints::local {

//=== Public Methods ===========================================================
Constraint_Sigma::Constraint_Sigma(
    const SymbolVector &symbol_vec,
    GapVector left_gaps,
    GapVector right_gaps) : BaseConstraint(ConstraintCategory::LOCAL, ConstraintType::SIGMA) {
  if (!left_gaps.empty()) {
    uint const shortest_length = std::min({symbol_vec.size(), left_gaps.size()});
    for (uint i = 0; i < shortest_length; i++)
      left_.emplace(symbol_vec[i], left_gaps[i]);
  }
  if (!right_gaps.empty()) {
    uint const shortest_length = std::min({symbol_vec.size(), right_gaps.size()});
    for (uint i = 0; i < shortest_length; i++)
      right_.emplace(symbol_vec[i], right_gaps[i]);
  }
}

Constraint_Sigma::Constraint_Sigma(SigmaTupleMap l, SigmaTupleMap r)
    : BaseConstraint(ConstraintCategory::LOCAL, ConstraintType::SIGMA),
      left_(std::move(l)),
      right_(std::move(r)) {}

std::string_view Constraint_Sigma::GetName() const {
  return {kName};
}

std::string_view Constraint_Sigma::GetDescription() const {
  return {kDescription};
}

std::string Constraint_Sigma::DebugString() const {
  std::stringstream oss;
  oss << kName << std::endl;
  oss << left_.size() << std::endl;
  for (const auto &[key, value] : left_) {
    oss << util::to_string(key) << " " << value.first << " " << value.second << std::endl;
  }
  oss << right_.size();
  for (const auto &[key, value] : right_) {
    oss << std::endl
        << util::to_string(key) << " " << value.first << " " << value.second;
  }
  return oss.str();
}

bool Constraint_Sigma::IsConstraintValid(const StrPtrVector &spv) const {
  if (left_.size() != right_.size()) {
    util::Logger::Warning() << "Constraint_Sigma: SigmaTupleMaps have different size!\n";
    return false;
  }
  if (!IsSigmaTupleMapValid(left_, spv)) {
    util::Logger::Warning() << "Constraint_Sigma: There's a problem with the left SigmaTupleMap\n";
    return false;
  }
  if (!IsSigmaTupleMapValid(right_, spv)) {
    util::Logger::Warning() << "Constraint_Sigma: There's a problem with the right SigmaTupleMap\n";
    return false;
  }
  return true;
}

bool Constraint_Sigma::IsEmbeddingValid(const Embedding& e) const {
  bool foo = DoLeftCheck(left_, e);
  foo = foo && DoRightCheck(right_, e);
  foo = foo && e.isValidEmbedding();
  return foo;
}

Constraint_Sigma::SigmaTupleMap Constraint_Sigma::GenSigmaTupleMap(
    SymbolSpan alph, GapSpan gc) {
  SigmaTupleMap gc_map;
  for (const auto&[symb, gap] : std::views::zip(alph, gc)) {
    gc_map[symb] = gap;
  }
  return gc_map;
}

Constraint_Sigma::SigmaTupleMap Constraint_Sigma::GenSigmaTupleMap(
    const SymbolVector &alph, const GapVector &gc) {
  SigmaTupleMap gc_map;
  for (const auto&[symb, gap] : std::views::zip(alph, gc)) {
    gc_map[symb] = gap;
  }
  return gc_map;
}

Constraint_Sigma::GapVector Constraint_Sigma::GenRelaxedGapVector(
    const SymbolVector& alph,
    const StrPtrVector &spv) {
  uint max = util::CalcMaxStrLen(spv);
  return {alph.size(), std::make_pair(0, max)};
}

Constraint_Sigma::GapVector Constraint_Sigma::GenRelaxedGapVector(
    const SymbolVector& alph,
    const std::span<const util::String> &s) {
  uint max = 0;
  for (const auto &p : s) {
    max = std::max<uint>(max, p.length());
  }
  return {alph.size(), std::make_pair(0, max)};
}

Constraint_Sigma Constraint_Sigma::CreateRelaxed(const StrPtrVector &spv) {
  auto l = SigmaTupleMap();
  auto r = SigmaTupleMap();

  std::set<Symbol> alphabet;
  for (const auto &p : spv) {
    for (const auto &symbol : *p) {
      alphabet.insert(symbol);
    }
  }

  uint max = util::CalcMaxStrLen(spv);
  for (const auto &symbol : alphabet) {
    l[symbol] = std::make_pair(0, max);
    r[symbol] = std::make_pair(0, max);
  }
  return {Constraint_Sigma(l, r)};
}

const Constraint_Sigma::SigmaTupleMap &Constraint_Sigma::GetLeft() const {
  return left_;
}
const Constraint_Sigma::SigmaTupleMap &Constraint_Sigma::GetRight() const {
  return right_;
}
void Constraint_Sigma::SetLeft(const SigmaTupleMap &l) {
  left_ = l;
}
void Constraint_Sigma::SetRight(const SigmaTupleMap &r) {
  right_ = r;
}

//=== Protected Methods ========================================================

Constraint_Sigma::Constraint_Sigma(
    SigmaTupleMap l,
    SigmaTupleMap r,
    const ConstraintType t)
    : BaseConstraint(ConstraintCategory::LOCAL, t)
    , left_(std::move(l)), right_(std::move(r)) {

}

Constraint_Sigma::Constraint_Sigma(
    const SymbolVector &symbol_vec,
    GapVector left_gaps,
    GapVector right_gaps,
    const ConstraintType t)
    : BaseConstraint(ConstraintCategory::LOCAL, t) {
  if (!left_gaps.empty()) {
    uint const shortest_length = std::min({symbol_vec.size(), left_gaps.size()});
    for (uint i = 0; i < shortest_length; i++)
      left_.emplace(symbol_vec[i], left_gaps[i]);
  }
  if (!right_gaps.empty()) {
    uint const shortest_length = std::min({symbol_vec.size(), right_gaps.size()});
    for (uint i = 0; i < shortest_length; i++)
      right_.emplace(symbol_vec[i], right_gaps[i]);
  }
}

/*******************************************************************************
 * isSigmaTupleMapValid
 * @param map SigmaTupleMap left or right map
 * @param spv StrPtrVector std::vector to pointers to the strings of the problem
 * @return true iff all strings in `spv` are consistent with the defined symbols
 *  in `map` and the values of map are within a valid range in relation to the
 *  strings
 ******************************************************************************/
bool Constraint_Sigma::IsSigmaTupleMapValid(
    const SigmaTupleMap &map,
    const StrPtrVector &spv) {
  if (spv.empty()) {
    return true;
  }

  if (map.empty()) {
    util::Logger::Warning() << kName << " Constraint: SigmaTupleMap is empty!\n";
    return false;
  }

  //Check: For every symbol there is a tuple in the map *p
  for (const auto &str_ptr : spv) {
    for (const auto &c : *str_ptr) {
      if (!map.contains(c)) {
        util::Logger::Warning() << kName
                                << " Constraint: For the symbol '"
                                << util::to_string(c)
                                << "' there is a gap constraint missing";
        return false;
      }
    }
  }

  // Calc the longest String length
  uint max = spv.front()->size();
  for (const auto &p : spv) {
    max = std::max<util::uint>(max, p->size());
  }

  //Check: tuples are in range
  for (const auto &[key, value] : map) {
    if (const auto &[l, u] = value; l > u || l > max || u > max) {
      util::Logger::Warning() << "isSigmaTupleMapValid: " << util::to_string(key)
                              << ", Pair: (" << l << ", " << u << ") is bad";
    }
  }
  return true;
}

/*******************************************************************************
 * doRightCheck
 * @param l SigmaTupleMap left map
 * @param e Embedding to check
 * @return true iff e satisfies all gap length constraints set by l
 * @note does not check whether e is a valid embedding
 ******************************************************************************/
bool Constraint_Sigma::DoLeftCheck(const SigmaTupleMap &l, const Embedding &e) {
  using util::Logger;
  if (l.empty()) {
    return true;
  }
  for (uint i = 0; i < e.size() - 1; i++) {
    // Calculate the number of Symbols between e.strView[i] and e.strView[i+1]
    const uint gap_length = (e[i + 1] - 1) - (e[i] + 1) + 1;
    const Symbol left_symbol = e.getSymbolAt(i);
    if (const auto &[lower, upper] = l.at(left_symbol);
        !(lower <= gap_length && gap_length <= upper)) {
      Logger::Warning() << " Embedding doLeftCheck: Failed! gap length (" << gap_length << ") "
                        << "does not match left[" << util::to_string(left_symbol)
                        << "] = (" << lower << ", " << upper << "), "
                        << "gap[" << i << "] = " << util::to_string(e.getGapStrView(i))
                        << "( [" << e[i] << ":" << e[i + 1] << "] in String: "
                        << util::to_string(e.getStr()) << ")";
      return false;
    }
  } /*i*/
  return true;
}

/*******************************************************************************
 * doRightCheck
 * @param r SigmaTupleMap right map
 * @param e Embedding to check
 * @return true iff e satisfies all gap length constraints set by r
 * @note does not check whether e is a valid embedding
 ******************************************************************************/
bool Constraint_Sigma::DoRightCheck(const SigmaTupleMap &r, const Embedding &e) {
  using util::Logger;
  if (r.empty()) return true;
  for (uint i = 0; i < e.size() - 1; i++) {
    // Calculate the number of Symbols between e.strView[i] and e.strView[i+1]
    const uint gap_length = (e[i + 1] - 1) - (e[i] + 1) + 1;
    const Symbol right_symbol = e.getSymbolAt(i + 1);
    if (const auto &[lower, upper] = r.at(right_symbol);
        !(lower <= gap_length && gap_length <= upper)) {
      Logger::Warning() << " Embedding doRightCheck: Failed! gap length (" << gap_length << ") "
                        << "does not match right[" << util::to_string(right_symbol) << "] = ("
                        << lower << ", " << upper << "), "
                        << "gap[" << i << "] = " << util::to_string(e.getGapStrView(i))
                        << "( [" << e[i] << ":" << e[i + 1] << "] in String: "
                        << util::to_string(e.getStr()) << ")";
      return false;
    }
  }
  return true;
}


}// namespace lcs_solver::constraints::local
//== Old Code for subsequence checking =========================================
// bool LCS_Sigma::isSubsequence(std::string sub, std::string str) {
//     int subLength = sub.length();
//     int strLength = str.length();
//     int i = 0; // Index for the subsequence
//     int j = 0; // Index for the string
//     int g = 0; // Count of characters between aligned characters
//     bool isFirstMatch = true; // Indicator for the first match
//     while (i < subLength && j < strLength) {
//         if (sub[i] == str[j]) {
//             if(!isFirstMatch){
//                 auto rGap = right.find(sub[i]);
//                 auto lGap = left.find(sub[i]);
//                 if(rGap != right.end() && (g < rGap->second.first || g > rGap->second.second) )
//                     return false;
//                 if(lGap != right.end() && (g < lGap->second.first || g > lGap->second.second) )
//                     return false;
//             } else {
//                 isFirstMatch = false;
//             }
//             ++i; // Move to the next character in the subsequence
//             g = 0; // Reset the gap
//         }
//         else {
//             if (!isFirstMatch) { // Count gaps only after the first match
//                 ++g; // Increase the gap
//             }
//         }
//         ++j; // Always move to the next character in the string
//     }
//     return i == subLength;  // Does not consider the gap after the last matching
// }