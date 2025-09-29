/*******************************************************************************
 * @file LLCS2_MC_INC.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Algorithm for LLCS2_MC_INC Problem. Uses MaxSegTreeND.
 * @details Time Complexity: O(log(m)*log(n)m*n)
 ******************************************************************************/
#include "algorithms/LLCS/LLCS2_MC_INC.h"

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_MC_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"
#include "structures/MaxSegTreeND.h"

#include <sstream>

namespace lcs_solver::algorithms::llcs {

/*******************************************************************************
* @brief Constructor for LLCS2_MC_INC
* @param vec The vector of shared points to the constant strings of a problem
* @param map std::map<std::string, share_ptr<BaseConstraint> for constraints
*******************************************************************************/
LLCS2_MC_INC::LLCS2_MC_INC(const StrPtrVector &vec, const ConstraintMap &map)
    : LLCS2_MC_Algorithm(
          AlgoType::LLCS2_MC_INC, vec, map,
          getGaps(map,ConstraintType::MC_INC),
          false, true),
      v(s.size() == 2 ? s[0] : util::StringView()),
      w(s.size() == 2 ? s[1] : util::StringView()),
      m(s.size() == 2 ? s[0].size() : 0),
      n(s.size() == 2 ? s[1].size() : 0) {
}

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS2_MC_INC::isValid() const {
  if (s.size() != 2) return false;
  if (!isSorted()) return false;
  if (!isFilled()) return false;
  if (!isEachConstraintIndividuallyValid()) return false;
  return true;
}

/*******************************************************************************
 * @brief Getter for the algorithms pseudocode
 * @return std::string_view of the description string
 ******************************************************************************/
std::string_view LLCS2_MC_INC::getDescription() const {
  return {R"(
LLCS Algorithm for the case |s1|-1 gc-constraint-pairs with the additional constraint that the gaps are increasing.
This algorithm uses a multidimensional segment tree. Time Complexity: O(log(m)*log(n)m*n)
/** Pseudocode (Adamson et. al, LCS wit Gap constraints page 12) ***************
 *  LLCS2_MC_1C(s1,s2,gc):
 *      // Initialize variables and datastructures
 *      let m = |s1|, n = |s2| (wlog n<=m)
 *      let M = Zeros[m + 1][n + 1]
 *      let T = Empty 2D-Max-Segment Tree over [0,1, .., m] x [0,1, .., n]
 *      let Phi(s1[i], s2[j]) = ith symbol of s1 is equal to jth symbol of s2
 *      for i in [1..l1]
 *          for j in [1..l2]
 *              if(Phi(s1[i],s2[j]) T.insert([i:i][j:j],1)
 *
 *      // Fill M using dynamic programming and O(m*n*log(n)*log(m))
 *      for i in [0..l1]
 *          for j = [0..l2]
 *              let p = T.query([i:i][j:j]) and (l,u) = gap[p]
 *              let s0=i+l+1, e0=i+u+1, s1=i+l+1, e1=i+u+1,
 *              update T[x][y] = max {M[x][y], p+1} for all (x,y) in [s0:e0]x[s1:e1]
 *              if Phi(s1[i], s2[j]) == 1 then
 *                  set M[i,j] = p
 *              else
 *                  set M[i][j] = max {M[i-1][j] , M[i][j-1])}
 *      return M[l1,l2]
 */
)"};
}

/*******************************************************************************
 * @brief Executes the query operation for the algorithm and returns a solution.
 *
 * @details This function resets the object's state and verifies its validity.
 * The function proceeds with preprocessing and returns solution pointer.
 *
 * @return std::unique_ptr<BaseSolution> A unique pointer to a solution object:
 * - `EmptySolution` if the object is not valid
 * - `UnsignedSolution` containing the lcs calculated
 ******************************************************************************/
std::unique_ptr<BaseSolution> LLCS2_MC_INC::query() {
  reset(ResetLevel::Full);
  if (!isValid()) {
    return std::make_unique<solutions::EmptySolution>();
  }
  doPreprocessing();
  //return std::make_unique<UnsignedSolution>(M[std::array<uint, 2>{m, n}]);
  return std::make_unique<solutions::UnsignedSolution>(max(M));
}

/*******************************************************************************
 * @brief Uses a MaxSegTreeND to the dp matrix `M`
 ******************************************************************************/
void LLCS2_MC_INC::doPreprocessing() {
  using std::min;
  using structures::MaxSegTreeND;
  M = Matrix(m + 1, Row(n + 1, 0));
  auto T = MaxSegTreeND({m + 1u, n + 1u});

  // Mark matching symbols in M[1:m,1:n] in the segment tree
  for (std::array<size_t, 2> idx = {1, 1}; idx[0] <= m; ++idx[0], idx[1] = 1) {
    for (; idx[1] <= n; ++idx[1]) {
      if (v[idx[0] - 1] == w[idx[1] - 1]) {
        T.update(idx, 1);
        if (trackKeyPairs) track({idx[0], idx[1]}, 1);
      }
    }
  }

  // Check if there are gc[0:p-1] subsequences & save result in M
  using Index = std::array<size_t, 2>;
  for (Index i = {1, 1}; i[0] <= m; ++i[0], i[1] = 1) {
    for (; i[1] <= n; ++i[1]) {
      if (v[i[0] - 1] == w[i[1] - 1]) {
        /* Note about the use of min:
         * If M[x]=M[y]=p such that T[i] was updated in the loop for x and y,
         * there had to be a guard against incrementing T[i] with p+1 and p+2
         */
        //uint const p = std::min(T.query(i), M[Index{i[0] - 1, i[1] - 1}] + 1);
        uint const p = T.query(i);
        M[i[0]][i[1]] = p;
        if (trackKeyPairs) track({i[0], i[1]}, M[i[0]][i[1]]);
        if (p > 0) {
          const auto &[l, u] = C[p - 1];
          std::array<std::pair<size_t, size_t>, 2> range;
          const uint s0 = i[0] + l + 1;
          const uint e0 = i[0] + u + 1;
          const uint s1 = i[1] + l + 1;
          const uint e1 = i[1] + u + 1;
          range[0] = {min(s0, m), min(e0, m)};
          range[1] = {min(s1, n), min(e1, n)};
          if (s0 <= m && s1 <= n)
            T.update(range, p + 1);
        }
      } else {
        // M[i] = std::max(M[Index{i[0] - 1, i[1]}], M[Index{i[0], i[1] - 1}]);
      }
    }// for i[1] in [1:n]
  }// for i[0] in [1:m]
  setState(State::Preprocessed);
}

/*******************************************************************************
 * @brief Resets the state of the LLCS2_MC_INC object to its initial condition
 * @details This function clears the matrix `M` and sets the object's state to
 *  `Constructed`.
 * @param lvl The level of reset to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_MC_INC::reset(BaseAlgorithm::ResetLevel lvl) {
  M.clear();
  setState(State::Constructed);
}

/*******************************************************************************
 * @brief Getter for the dp-Matrix
 * @return std::vector<std::vector<uint> M
 *******************************************************************************/
const LLCS2_Algorithm::Matrix &LLCS2_MC_INC::getMatrix() const {
  return M;
}

/*******************************************************************************
 * @brief creates a std::string to describe the state of the algorithm
 * @return std::string containing the dp matrix `M`
 ******************************************************************************/
std::string LLCS2_MC_INC::DebugString() const {
  std::ostringstream oss;
  oss << BaseAlgorithm::toString(s) << "\n";
  oss << toString(M, s, false, "M", true) << "\n";
  oss << toString(keyPairs, trackKeyPairs) << "\n";
  return oss.str();
}

}// namespace lcs_solver::algorithms::llcs