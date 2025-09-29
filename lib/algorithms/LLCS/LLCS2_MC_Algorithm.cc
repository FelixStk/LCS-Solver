/*******************************************************************************
 * @file LLCS2_MC_Algorithm.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Implementation of a superclass for LLCS problems with alphabet
 *  independent gap tuples.
 ******************************************************************************/
#include "algorithms/LLCS/LLCS2_MC_Algorithm.h"

#include <cassert>
#include <stdexcept>

#include "algorithms/AlgoType.h"
#include "algorithms/LLCS/LLCS_Algorithm.h"
#include "constraints/local/Constraint_MC.h"

#include <map>

namespace lcs_solver::algorithms::llcs {

/*******************************************************************************
 * @brief Constructor of LLCS2_MC_Algorithm
 * @details provides access to the gc tuples via `C` and the maximum number of
 *  gaps possible as determined by the shortest string length
 * @param algo AlgoType identifies the algorithm uniquely
 * @param vec StrPtrVector of shared points to the constant strings of a problem
 * @param map std::map<std::string, share_ptr<BaseConstraint> for constraints
 * @param gap std::vector<std::Pair<unit,unit>> contains the gap tuples
 * @param doTracking true, iff the LLCS2_Algorithm algorithm is supposed to do
 *  bookkeeping of points that are extendable: If Pair x is in keyPoints[l],
 *  then there exists a shared subsequence of length l that ends at x (meaning
 *  it is a subsequence of s1[0:x.first-1] and s2[0:x.second-1]
 * @param oneBased Whether the indices used in the matrices of the algorithm are
 *  oneBased (in [1:s1.size()]x[1:s2.size()] or zeroBased
 ******************************************************************************/
LLCS2_MC_Algorithm::LLCS2_MC_Algorithm(
    const AlgoType algo,
    const StrPtrVector &vec,
    const ConstraintMap &map,
    const GapVector &gap,
    const bool doTracking,
    const bool oneBased)
    : LLCS2_Algorithm(algo, vec, map, doTracking, oneBased),
      nMaxGaps(s.empty() ? 0 : s[0].size() - 1),//assumes s is sorted
      C(gap) {}

/*******************************************************************************
 * Extracts a K-Tuple for local, sigma-independent, gap-constraint
 * @param map map<std::string, std::shared_ptr<BaseConstraint> for constraints
 * @param type gap-constraint to look for
 * @return GapVector& ref to the std::vector<std::pair<uint, uint>> of
 ******************************************************************************/
const LLCS2_MC_Algorithm::GapVector &LLCS2_MC_Algorithm::getGaps(
    const ConstraintMap &map,
    ConstraintType type) {
  const auto key = type;
  const auto &entry = map.at(key);
  auto *p = dynamic_cast<constraints::local::Constraint_MC *>(entry.get());
  if (p == nullptr)
    throw std::runtime_error("LLCS2_MC_Algorithm::getGaps - Type cast failed");

  return p->GetGapVector();
}

/*******************************************************************************
 * isExtensible
 * @param a Pair (a1,a2) a one based position: s1[a1],s2[a2]
 * @param b Pair (b1,b2) a one based position: s1[b1],s2[b2]
 * @param llcsOfA Uint equal to LLCS(s1[1:a1],s2[1:a2])
 * @return Whether an LCS from position a can be extended to position b
 ******************************************************************************/
bool LLCS2_MC_Algorithm::isExtensible(Pair a, Pair b, uint llcsOfA) const {
  const auto &[ax, ay] = a;
  const auto &[bx, by] = b;

  const bool match = isMatched(a) && isMatched(b);
  const bool validPairs = bx > ax && by > ay;
  const bool llcsOfAValid = llcsOfA > 0 && (llcsOfA - 1) < static_cast<BaseAlgorithm::uint>(C.size());

  // Allow Extension of the Zero Pair
  // if (llcsOfA == 0 && ax == 0 && by == 0) {
  //   return isMatched(b) && validPairs;
  // }

  if (match && validPairs && llcsOfAValid) {
    const auto &[l, u] = C[llcsOfA - 1];
    const uint lenX = bx - ax - 1;// gap has len == 0 if b.first == (a.first - 1)
    const uint lenY = by - ay - 1;
    const bool gapX = lenX >= l && lenX <= u;
    const bool gapY = lenY >= l && lenY <= u;
    return gapX && gapY;
  }
  return false;
}

/*******************************************************************************
 * getPrevRange
 * @param pair Pair (x,y) a one based position: s1[x],s2[y]
 * @param llcs The length of the LCS at position (x,y)
 * @return (s1Factor, s2Factor) where a factor is a one based uint pair. It
 *  represents all possible positions where a positions (i,j) that can be
 *  extended to (x,y) can be found.
 ******************************************************************************/
LLCS2_Algorithm::Window LLCS2_MC_Algorithm::getPrevRange(const Pair& pair,
                                                         uint llcs) const {
  const auto& [x, y] = pair;
  Pair r1, r2;
  assert(llcs > 0);
  if (llcs == 1) {
    // The gap before the first position (with LLCS(x,y)==1) is not considered
    r1 = {1, x - 1};
    r2 = {1, y - 1};
  } else {
    // The constraint is first active for LLCS(x,y)==2 with C[0]
    assert(llcs <= C.size() + 2);
    const auto& gapBounds = C[llcs - 2];
    const auto& [l, u] = gapBounds;
    auto sub = [](uint a, uint b) { return a > b ? a - b : 1; };
    r1 = {sub(x, u + 1), sub(x, l + 1)};
    r2 = {sub(y, u + 1), sub(y, l + 1)};
  }
  return {r1, r2};
}

/*******************************************************************************
 * StringFrom
 * @return std::string representation of the gap vector
 ******************************************************************************/
std::string LLCS2_MC_Algorithm::StringFrom(const GapVector& gc) {
  std::ostringstream oss;
  for (uint i = 0; i < gc.size(); ++i) {
    oss << "C[" << i << "] = " << gc[i].first << "," << gc[i].second;
    if (i != gc.size() - 1) {
      oss << "\n";
    }
  }
  return oss.str();
}

/*******************************************************************************
 * Getter for the vector of gap length bounds
 * @return const reference to protected GapVector C
 ******************************************************************************/
const LLCS2_MC_Algorithm::GapVector &LLCS2_MC_Algorithm::getGaps() const {
  return C;
}

}// namespace lcs_solver::algorithms::llcs