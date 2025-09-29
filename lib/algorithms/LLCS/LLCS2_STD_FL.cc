/*******************************************************************************
 * @file LLCS2_STD_FL.cc
 * @author Steinkopp:Felix
 * @brief Folklore LCS Algorithm for Two Strings
 * @details Time Complexity: \f$ \mathcal O \left( |s_1| \cdot |s_2|)\f$
 *          Space Complexity: \f$ \mathcal O \left( |s_1| \cdot |s_2|)\f$
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_STD_FL.h"

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"

namespace lcs_solver::algorithms::llcs {

constexpr std::string_view LLCS2_STD_FL::description = R"DESC(Pseudocode: LCS Folklore Algorithm for n strings
> LLCS(s1,s2,...,sk):
>   let li = length(si) and li<=lj for all i<j
>   let dp[l1 + 1][l2 + 1]...[lk + 1] filled with zeros
>   for i1 = 1 to l1 do
>     for i2 = 1 to l2 do
>       if s1[i1-1] == s2[i2-1] == ... == sk[ik-1]
>         dp[i1][i2]...[ik] = dp[i1-1][i2-1] + 1
>       else
>         dp[i1][i2]...[ik] = max( dp[i1-1][i2], dp[i1][i2-1] )
>   return dp[l1][l2]...[lk]
)DESC";

/*******************************************************************************
 * Constructor for LLCS2_STD_FL
 * @param vec The vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 ******************************************************************************/
LLCS2_STD_FL::LLCS2_STD_FL(
    const BaseAlgorithm::StrPtrVector &vec,
    const ConstraintMap &map)
    : LLCS2_Algorithm(AlgoType::LLCS2_STD_FL, vec, map, false, true) {}

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS2_STD_FL::isValid() const {
  bool const resIsSorted = isSorted();
  bool const resIsFilled = isFilled();
  bool const resEachValid = isEachConstraintIndividuallyValid();
  return resIsSorted && resIsFilled && resEachValid;
}

/*******************************************************************************
 * @brief Getter for the algorithms pseudocode
 * @return std::string_view of the description string
 ******************************************************************************/
std::string_view LLCS2_STD_FL::getDescription() const {
  return description;
}

/*******************************************************************************
 * @brief Executes the query operation for the algorithm and returns a solution.
 * @details This function resets the object's state and verifies its validity.
 * The function proceeds with preprocessing and returns a solution pointer.
 * @return std::unique_ptr<BaseSolution> A unique pointer to a solution object:
 * - `EmptySolution` if the object is not valid or zero strings are provided
 * - `UnsignedSolution` with s[0].size if only one problem string is provided
 * - `UnsignedSolution` containing the llcs calculated
 ******************************************************************************/
std::unique_ptr<BaseSolution> LLCS2_STD_FL::query() {
  reset(ResetLevel::Full);
  if (s.empty() || !isValid()) {
    return std::make_unique<solutions::EmptySolution>();
  }
  if (s.size() == 1) {
    return std::make_unique<solutions::UnsignedSolution>(s[0].size());
  }
  doPreprocessing();
  // uint lastElem = dp.empty() ? 0 : dp.back().back();
  return std::make_unique<solutions::UnsignedSolution>(max(dp));
}

/*******************************************************************************
 * Resets the state of the LLCS2_STD_FL object to its initial condition.
 * @note Side effect: Erases the elements from M and frees memory
 * @param lvl BaseAlgorithm::ResetLevel to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_STD_FL::reset(BaseAlgorithm::ResetLevel lvl) {
  LLCS2_Algorithm::reset(lvl);
  dp.clear();
}

/*******************************************************************************
 * @brief Executes the folklore llcs algorithm generalized to k strings
 ******************************************************************************/
void LLCS2_STD_FL::doPreprocessing() {
  const uint l1 = s.size() < 2 ? 0 : s[0].size();
  const uint l2 = s.size() < 2 ? 0 : s[1].size();
  dp = Matrix(l1 + 1, std::vector<uint>(l2 + 1));

  // Folklore Algorithm
  Pair idx = {1, 1};
  for (uint &i = idx.first; i < l1 + 1; ++i, idx.second = 1) {
    for (uint &j = idx.second; j < l2 + 1; ++j) {
      if (isMatched(idx, true)) {
        dp[i][j] = dp[i - 1][j - 1] + 1;
        if (trackKeyPairs)
          track({i, j}, dp[i][j]);
      } else
        dp[i][j] = std::max(dp[i - 1][j], dp[i][j - 1]);
    }
  }
  // If wanted, adapt the shared dp format
  if (trackKeyPairs) {
    idx = {1, 1};
    for (uint &i = idx.first; i <= l1; ++i) {
      dp[i][1] = isMatched(idx, true);
    }
    idx = {1, 1};
    for (uint &j = idx.second; j <= l2; ++j) {
      dp[1][j] = isMatched(idx, true);
    }
    if (l1 > 1 && l2 > 1) {
      Pair end = {2, 2};
      for (uint &i = end.first; i < l1 + 1; ++i, end.second = 2) {
        for (uint &j = end.second; j < l2 + 1; ++j) {
          if (!isMatched(end, true)) {
            dp[i][j] = 0;
          }
        }
      }
    }
  }
  setState(State::Preprocessed);
}

/*******************************************************************************
 * @brief Creates a std::string to describe the state of the algorithm
 * @return std::string containing:
 * - the strings of the problem
 * - the dp matrix
 ******************************************************************************/
std::string LLCS2_STD_FL::DebugString() const {
  std::stringstream oss;
  oss << BaseAlgorithm::toString(s) << "\n";
  oss << toString(dp,s, true, "dp", true) << "\n";
  oss << toString(keyPairs, trackKeyPairs);
  return oss.str();
}

/*******************************************************************************
 * Getter for the dp-Matrix
 * @return structures::Matrix<Unsigned> dp
 * @note dp is modified in doPreprocessing() to fulfill the common format:
 *  dp[i][j] = isMatched(i,j) ? llcs(s1[:i],s2[:j]) : 0
 *  If the flag is false dp is the monotone increasing dp matrix of the folklore
 *  algorithm.
 *******************************************************************************/
const LLCS2_STD_FL::Matrix &LLCS2_STD_FL::getMatrix() const {
  return dp;
}

/**
 * isExtensible
 * @param a Start Index (to be extended south-east to b)
 * @param b End Index (b1,b2)
 * @param llcsOfA Uint equal to LLCS(s1[1:a1],s2[1:a2])
 * @return Whether an LCS from position a can be extended to position b
 */
bool LLCS2_STD_FL::isExtensible(const Pair a, const Pair b, uint llcsOfA) const {
  const bool match = isMatched(a) && isMatched(b);
  const bool gapX = (b.first > a.first) && (b.first - a.first >= 1);
  const bool gapY = (b.second > a.second) && (b.second - a.second >= 1);
  return match && gapX && gapY;
}

/**
 * getNextWindow
 * @param pair Index (x,y) for specifying the position s[1:x], s[1:y]
 * @param llcs The length of longest common subsequence of s[1:x] and s[1:y]
 * @return {[1:x-1], [1:y-1]}
 */
LLCS2_Algorithm::Window LLCS2_STD_FL::getPrevRange(const Pair &pair, uint llcs) const {
  return {{1, pair.first - 1}, {1, pair.second - 1}};
}

}// end of namespace lcs_solver::algorithms::llcs
