/*******************************************************************************
 * @file LLCS2_MC_INC_E.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Algorithm for LLCS2_MC_INC Problem. Uses MaxDequeCR.
 * @details Time Complexity: O(m*n)
 ******************************************************************************/
#include "algorithms/LLCS/LLCS2_MC_INC_E.h"

#include <cassert>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_MC_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"
#include "structures/MaxDequeCR.h"

namespace {
using ::lcs_solver::algorithms::solutions::EmptySolution;
using ::lcs_solver::algorithms::solutions::UnsignedSolution;
using ::lcs_solver::structures::MaxDequeCR;
}// namespace

namespace lcs_solver::algorithms::llcs {

/*******************************************************************************
* @brief Constructor for LLCS2_MC_INC_E
* @param vec The vector of shared points to the constant strings of a problem
* @param map std::map<std::string, share_ptr<BaseConstraint> for constraints
*******************************************************************************/
LLCS2_MC_INC_E::LLCS2_MC_INC_E(const StrPtrVector &vec, const ConstraintMap &map)
    : LLCS2_MC_Algorithm(
          AlgoType::LLCS2_MC_INC_E, vec, map,
          getGaps(map, ConstraintType::MC_INC),
          false,// Do not do tracking by default
          false // In this Algorithm the dp Matrix M is zero-Based
          ),
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
bool LLCS2_MC_INC_E::isValid() const {
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
std::string_view LLCS2_MC_INC_E::getDescription() const {
  return {R"(
LLCS Algorithm for the case |s1|-1 gc-constraint-pairs with the additional constraint that the gaps are increasing.
This algorithm uses a multidimensional segment tree. Time Complexity: O(log(m)*log(n)m*n)
/** Pseudocode *****************************************************************
 *  LLCS2_MC_1C(s1,s2,gc):
 *      let m = |s1|, n = |s2| (wlog n<=m)
 *      let M = Zeros[m + 1][n + 1]
 *      let Phi(s1[i], s2[j]) = ith symbol of s1 is equal to jth symbol of s2
 *      let D = MaxDequeCR(M, Phi, lowerGapBoundGetter, upperGapBoundGetter)
 *      for i in [1..l1]
 *          for j in [1..l2]
 *              D.update(i,j)
 *              let Phi(i,j) = {M[ip,jp] | ip + Tl(M[ip,jp]) + 1 <= i <= ip + Tu(M[ip,jp] + 1 and
 *                                         jp + Tl(M[ip,jp]) + 1 <= j <= ip + Tu(M[ip,jp] + 1 and
 *                                         with ip in [i] and jp in [j] }
 *              set M[i,j] := max Phi(i,j) = D.extractMax(i,j)
 *              D.insert(i,j)
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
std::unique_ptr<BaseSolution> LLCS2_MC_INC_E::query() {
  reset(ResetLevel::Full);
  if (!isValid()) {
    return std::make_unique<EmptySolution>();
  }
  doPreprocessing();
  return std::make_unique<UnsignedSolution>(max(M));
}

/*******************************************************************************
 * @brief Getter for the dp-Matrix
 * @return std::vector<std::vector<uint> M
 *******************************************************************************/
const LLCS2_Algorithm::Matrix &LLCS2_MC_INC_E::getMatrix() const {
  return M;
}

/*******************************************************************************
 * @brief Uses a MaxSegTreeND to the the dp matrix `M`
 ******************************************************************************/
void LLCS2_MC_INC_E::doPreprocessing() {
  M = Matrix(m, std::vector<uint>(n, 0));
  auto Phi = [this](uint i, uint j) {
    assert(i < m && j < n);
    return v[i] == w[j] ? 1 : 0;
  };
  auto Tl = [this](uint x) {
    assert(x < m);
    return C[x].first;
  };
  auto Tu = [this](uint x) {
    assert(x < n);
    return C[x].second;
  };

  MaxDequeCR D{M, Phi, Tl, Tu};
  for (uint i = 0; i < m; ++i) {
    for (uint j = 0; j < n; ++j) {
      D.update(i, j);
      if (Phi(i, j) == 1) {
        M[i][j] = D.extractMax(i, j);
        if (trackKeyPairs)
          track({i+1, j+1}, M[i][j]); // +1 because of zero-based indexing
        D.insert(i, j);
      }
    }
  }
  setState(State::Preprocessed);
}

/*******************************************************************************
 * @brief Resets LLCS2_MC_INC_E object to the state before preprocessing
 * @details This function clears the matrix `M` and sets the object's state to
 *  `Constructed`.
 * @param lvl The level of reset to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_MC_INC_E::reset(BaseAlgorithm::ResetLevel lvl) {
  M.clear();
  setState(State::Constructed);
}

/*******************************************************************************
 * @brief creates a std::string to describe the state of the algorithm
 * @return std::string containing the dp matrix `M`
 ******************************************************************************/
std::string LLCS2_MC_INC_E::DebugString() const {
  std::ostringstream oss;
  oss << BaseAlgorithm::toString(s) << "\n\n";
  oss << StringFrom(C) << "\n\n";
  oss << toString(M, s, false, "M", false) << "\n";// Does not support zero Based String Annotation
  oss << toString(keyPairs, trackKeyPairs) << "\n";
  return oss.str();
}

}// namespace lcs_solver::algorithms::llcs