/******************************************************************************
 * @file LLCS2_SIG_Algorithm.cc
 * @author Steinkopp:Felix
 * @version 1.2
 * @brief Specification of BaseAlgorithm for LLCS Problems
 *****************************************************************************/
#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"

#include <algorithm>
#include <ranges>
#include <set>
#include <stdexcept>

#include "algorithms/AlgoType.h"
#include "algorithms/LLCS/LLCS_Algorithm.h"
#include "constraints/local/Constraint_Sigma.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::llcs {

/*******************************************************************************
 * @brief Constructor for LLCS2_SIG_Algorithm objects
 * @param algo AlgoType identifies the algorithm uniquely
 * @param vec StrPtrVector of shared points to the constant strings of a problem
 * @param map std::map<std::string, share_ptr<BaseConstraint> for constraints
 * @param l std::unordered_map<Symbol,std::Pair<uint,unit> left mapping with the
 * meaning of l[c] = (l,u) => [l ≤ length(gap before Symbol c) ≤ u]
 * @param r std::unordered_map<Symbol,std::Pair<uint,unit> right mapping with
 * the meaning of r[c] = (l,u) => [l ≤ length(gap after Symbol c) ≤ u]
 ******************************************************************************/
LLCS2_SIG_Algorithm::LLCS2_SIG_Algorithm(
    AlgoType algo,
    const StrPtrVector &vec,
    const ConstraintMap &map,
    const SigmaTupleMap &l,
    const SigmaTupleMap &r)
    : LLCS2_Algorithm(algo, vec, map),
      left(l),
      right(r),
      sig(std::max(l.size(), r.size())) {
}

/*******************************************************************************
 * @brief Gets a reference to the left SigmaTupleMap of a local, sigma
 * dependent, gap-constraint by looking at the ConstraintMap `map`.
 * @param map std::map from ConstraintType to std::shared_ptr of the constraint
 * @param type KnownConstraint to look for in map. Should be Sigma_R or Sigma.
 * @return left[b] mapping for bounds of the gap length after a symbol b
 ******************************************************************************/
const LLCS2_SIG_Algorithm::SigmaTupleMap &LLCS2_SIG_Algorithm::getSigLMap(
    const ConstraintMap &map,
    ConstraintType type) {
  auto key = type;
  const auto &entry = map.at(key);
  using ::lcs_solver::constraints::local::Constraint_Sigma;
  auto *p = dynamic_cast<Constraint_Sigma *>(entry.get());
  if (p == nullptr) {
    throw std::runtime_error("LLCS2_SIG_Algorithm::getSigLMap - Type cast failed");
  }
  return p->GetLeft();
}

/*******************************************************************************
 * @brief Gets a reference to the right SigmaTupleMap of a local, sigma
 * dependent, gap-constraint by looking at the ConstraintMap `map`.
 * @param map std::map from ConstraintType to std::shared_ptr of the constraint
 * @param type KnownConstraint to look for in map. Should be Sigma_R or Sigma.
 * @return right[a] mapping for bounds of the gap length before a symbol a
 ******************************************************************************/
const LLCS2_SIG_Algorithm::SigmaTupleMap &LLCS2_SIG_Algorithm::getSigRMap(
    const ConstraintMap &map,
    ConstraintType type) {
  auto key = type;
  const auto &entry = map.at(key);
  using ::lcs_solver::constraints::local::Constraint_Sigma;
  auto *p = dynamic_cast<Constraint_Sigma *>(entry.get());
  if (p == nullptr) {
    throw std::runtime_error("Type cast failed");
  }
  return p->GetRight();
}

/*******************************************************************************
 * @brief Generates a std::vector of all symbols that are found in `map`
 * @param map SigmaTupleMap - the right or left function from the paper
 * @return std::vector with all keys / symbols in map
 ******************************************************************************/
std::vector<lcs_solver::util::Symbol> LLCS2_SIG_Algorithm::getAlphabet(
    const SigmaTupleMap &map) {
  std::set<lcs_solver::util::Symbol> s;
  for (const auto &key : map | std::views::keys) {
    s.insert(key);
  }
  return {s.begin(), s.end()};
}

/*******************************************************************************
 * @brief Generates a lookUP Object that combines the left and right maps
 * @details A gap may have a right and a left symbol at its borders. This
 * function calculates the resulting bounds on the gap length when both the left
 * and  right functions pose a constraint on a gap
 * @param left with left[b] = bounds of the gap length after a symbol b
 * @param right with right[a] = bounds of the gap length before a symbol a
 * @return SymbSymbTupleMap, where the first key is for the left map symbol and
 * the second key is for the right key symbol
 ******************************************************************************/
LLCS2_SIG_Algorithm::SymbSymbTupleMap LLCS2_SIG_Algorithm::calcGapIntersection(
    const SigmaTupleMap &left,
    const SigmaTupleMap &right) {
  SymbSymbTupleMap result;
  for (const auto &[lKey, lGap] : left) {
    for (const auto &[rKey, rGap] : right) {
      result[lKey][rKey] =
          std::make_pair(
              std::max(lGap.first, rGap.first),
              std::min(lGap.second, rGap.second));
    }
  }
  return result;
}

/*******************************************************************************
 * Helper to avoid code duplication in isExtensible functions
 * @param a Pair (a1,a2) a one based position: s1[a1],s2[a2]
 * @param b Pair (b1,b2) a one based position: s1[b1],s2[b2]
 * @param ignoreRight Whether to ignore a constraint posed by right
 * @param ignoreLeft Whether to ignore a constraint posed by left
 * @note Effect of left: If 'x' + gap ++ 'y', then 'x' influences len(gap)
 * @note Effect of right: If 'x' + gap ++ 'y', then 'y' influences len(gap)
 * @return Whether an LCS from position a can be extended to position b
 ******************************************************************************/
bool LLCS2_SIG_Algorithm::isExtensibleHelper(
    Pair a, Pair b, bool ignoreRight, bool ignoreLeft) const {
  const bool matched = isMatched(a) && isMatched(b);
  if (!matched || s.size() != 2)
    return false;

  const Symbol c1 = s[0].at(a.first-1);
  const Symbol c2 = s[0].at(b.first-1);
  const auto &[lr, ur] = right.at(c2);
  const auto &[ll, ul] = left.at(c1);
  const auto &[ax, ay] = a;
  const auto &[bx, by] = b;
  const uint lenX = bx > ax ? bx - ax - 1 : 0;
  const uint lenY = by > bx ? by - ay - 1 : 0;
  bool gapX = true, gapY = true;
  if (!ignoreRight && !ignoreLeft) {
    gapX = lenX >= lr && lenX <= ur && lenX >= ll && lenX <= ul;
    gapY = lenY >= lr && lenY <= ur && lenY >= ll && lenY <= ul;
  }
  else if (!ignoreRight) {
    gapX = lenX >= lr && lenX <= ur;
    gapY = lenY >= lr && lenY <= ur;
  }
  else if (!ignoreLeft) {
    gapX = lenX >= ll && lenX <= ul;
    gapY = lenY >= ll && lenY <= ul;
  }

  return gapX && gapY;
}

/**
 * Helper to avoid code duplication in getPrevRange
 * @param pair Pair (x,y) a one based position: s1[x],s2[y]
 * @param c Symbol at s[0].at(x) and s[1].at(y)
 * @param right SigmaTupleMap for the sigma right constraint
 * @return (s1Factor, s2Factor) where a factor is a one based uint pair. It
 *  represents a superset of possible positions where a position (i,j) in it
 *  can be extended to (x,y) under the right sigma constraint
 */
LLCS2_Algorithm::Window LLCS2_SIG_Algorithm::getPrefRange(const Pair &pair, Symbol c, const SigmaTupleMap &right) {
  const auto &[x,y] = pair;
  const auto &[l, u] = right.at(c);
  auto sub = [](uint a, uint b) { return a > b ? a - b : 1; };
  Pair r1 = {sub(x, u + 1), sub(x, l + 1)};
  Pair r2 = {sub(y, u + 1), sub(y, l + 1)};
  return {r1,r2};
}

}// namespace lcs_solver::algorithms::llcs