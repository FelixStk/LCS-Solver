/*******************************************************************************
 * @file LLCS2_MC_1C.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Implementation of an algorithm for LCS2_MC_1C. Uses MaxQueue2D.
 * @details Time Complexity: O(m*n)
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_MC_1C.h"

#include <algorithm>
#include <functional>// defining phi as a std::function
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_MC_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"
#include "structures/MaxQueueN.h"

template struct lcs_solver::structures::MaxQueue2D<::lcs_solver::algorithms::BaseAlgorithm::uint>;

namespace lcs_solver::algorithms::llcs {

/*** Constructor ***************************************************************
 * @brief Constructor for LLCS2_MC_1C
 * @details Might add a trivial Constraint_1C of (0,|s[0]|) Pairs if necessary
 * @param vec The vector of shared points to the constant strings of a problem
 * @param map std::map<std::string, share_ptr<BaseConstraint> for constraints
 ******************************************************************************/
LLCS2_MC_1C::LLCS2_MC_1C(const StrPtrVector &vec, const ConstraintMap &map)
    : LLCS2_MC_Algorithm(AlgoType::LLCS2_MC_1C, vec,
                         map,
                         LLCS2_MC_Algorithm::getGaps(
                             map,
                             ConstraintType::MC_1C)),
      v(s.size() == 2 ? s[0] : util::StringView()),
      w(s.size() == 2 ? s[1] : util::StringView()),
      m(s.size() == 2 ? s[0].size() : 0),
      n(s.size() == 2 ? s[1].size() : 0),
      l(C.empty() ? 0 : C[0].first),
      u(C.empty() ? std::max(m, n) : C[0].second) {
}

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS2_MC_1C::isValid() const {
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
std::string_view LLCS2_MC_1C::getDescription() const {
  return {R"(
LLCS Algorithm for the case that all gc-constraint-pairs the same.
/** Pseudocode (Adamson et. al, LCS wit Gap constraints page 13-14) ************
 *  LLCS2_MC_1C(s1,s2,gc):
 *      let li = length(si) and li<=lj for all i<j
 *
 *      // Initialize a k-D table to store lengths of LCS for subproblems
 *      let M = Zeros[lk + 1][lk + 1]
 *      let l,u = C[0] = C[1] = ... = C[l1-2]
 *      let s0 = u+1, e0=l+1, s1=u+1, e1=l+1
 *      let Q = MaxQueue2D(const MatrixPtr = &M, xRange = [0..l1], yRange[0..l2], bufferAxisSize = u-l)
 *      let Phi(s1[i], s2[j]) = ith symbol of s1 is equal to jth symbol of s2
 *
 *      // Fill M_p using dynamic programming and a window sliding approach. M_p is computed in O(l1*l2)
 *      for i in [0..l1]
 *          update Q (set up for processing line i: M[i-e0][0:n] enters Q)
 *          for j = [0..l2]
 *              update Q (set up for computing M[i][j]: M[i-e0][j-e1] enters the Q)
 *              use Q to retrieve m, the maximum of the sub-matrix M[I][J] where I = [i-s0:i-e0] and J = [j-s1:j-e1]
 *              m is set to be 0 when I or J are empty
 *              if Phi(s1[i], s2[j]) == 1 then
 *                  set M[i,j] = m + 1
 *              else
 *                  set M[i][j] = 0
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
 * - `UnsignedSolution` containing the llcs calculated
 ******************************************************************************/
std::unique_ptr<BaseSolution> LLCS2_MC_1C::query() {
  reset(ResetLevel::Full);
  if (!isValid()) {
    return std::make_unique<solutions::EmptySolution>();
  }
  doPreprocessing();
  return std::make_unique<solutions::UnsignedSolution>(max(M));
}

/*******************************************************************************
 * @brief does the preprocessing step of the algorithm
 ******************************************************************************/
void LLCS2_MC_1C::doPreprocessing() {
  M = std::vector<std::vector<uint>>(
      m + 1,
      std::vector<uint>(n + 1, 0));
  std::function<bool(uint, uint)> const Phi = [this](uint i, uint j) -> uint {
    if (i == 0 || i > v.size()) return 0;
    if (j == 0 || j > w.size()) return 0;
    return v[i - 1] == w[j - 1];
  };

  auto MQ2D = lcs_solver::structures::MaxQueue2D<uint>(M);

  // The window(i,j) is M[i - s0 : i - e0][j - s1 : j - e1 ]
  const uint s0 = u + 1;
  const uint e0 = l + 1;
  const uint s1 = u + 1;
  const uint e1 = l + 1;
  const uint n0 = s0 - e0 + 1;// Number of elements in [i-s0:i-e0]
  // const uint n1 = s1 - e1 + 1; // Number of elements in [j-s1:j-e1]

  // Set min capacity greater zero so that Q[f] can store M[i-e0][j-e1]
  MQ2D.setYBuffer(std::max<uint>(1, n0 - 1));
  for (uint i = 0; i <= m; ++i) {
    MQ2D.doLineSetup(i, s0, e0);// Updates Max and RMQ for new row M[i-e0][0:n]
    for (uint j = 0; j <= n; ++j) {
      MQ2D.doComputingSetup(i, j, e0, e1);// Insert M[i-e0][j-e1]
      uint const llcs = MQ2D.query(i, j, {s0, e0, s1, e1});

      if (Phi(i, j)) {
        M[i][j] = llcs + 1;
        if (trackKeyPairs) track({i, j}, M[i][j]);
      } else if (i > 0 && j > 0) {
        // M[i][j] = std::max(M[i - 1][j], M[i][j - 1]);
      }
    }
  }

  setState(State::Preprocessed);
}

/*******************************************************************************
 * @brief Resets the state of the LLCS2_MC_1C object to its initial condition.
 * @details This function clears the matrix `M` and sets the object's state to
 *  `Constructed`.
 * @param lvl The level of reset to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_MC_1C::reset(BaseAlgorithm::ResetLevel lvl) {
  M.clear();
  setState(State::Constructed);
}

/*******************************************************************************
 * @brief Getter for the dp-Matrix
 * @return std::vector<std::vector<uint> M
 *******************************************************************************/
const LLCS2_Algorithm::Matrix &LLCS2_MC_1C::getMatrix() const {
  return M;
}

/*******************************************************************************
 * @brief creates a std::string to describe the state of the algorithm
 * @return std::string containing:
 * - the dp matrix `M`
 ******************************************************************************/
std::string LLCS2_MC_1C::DebugString() const {
  std::ostringstream oss;
  oss << BaseAlgorithm::toString(s) << "\n";
  oss << toString(M, s, false, "M", true) << "\n";
  oss << toString(keyPairs, trackKeyPairs) << "\n";
  return oss.str();
}

}// namespace lcs_solver::algorithms::llcs