/*******************************************************************************
 * @file LLCS2_MC.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Implementation of an algorithm for LCS2_MC
 * @details Time Complexity: \f$ \mathcal O \left( |v| \cdot |w| \cdot llcs \right) \f$
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_MC.h"

#include <cassert>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_MC_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"
#include "util/Logger.hpp"

namespace lcs_solver::algorithms::llcs {

/*******************************************************************************
 * LLCS2_MC Constructor
 * @param vec std::vector of shared points to constant strings for the problem
 * @param map A map<std::string, share_ptr<BaseConstraint> with the constraints
 *  of a problem
 *******************************************************************************/
LLCS2_MC::LLCS2_MC(const StrPtrVector &vec, const ConstraintMap &map)
    : LLCS2_MC_Algorithm(
          AlgoType::LLCS2_MC,
          vec,
          map,
          LLCS2_MC_Algorithm::getGaps(
              map,
              ConstraintType::MC),
          false, true) {
  k = 0;
  m = 0;
  n = 0;
  M.clear();
  A.clear();
  B.clear();
  setState(State::Constructed);
}

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS2_MC::isValid() const {
  return s.size() == 2
      && isSorted()// strings in s are increasing in their length
      && isFilled()// no string is empty
      && isEachConstraintIndividuallyValid();
}

/*******************************************************************************
 * @brief Getter for the algorithms pseudocode
 * @return std::string_view of the description string
 ******************************************************************************/
std::string_view LLCS2_MC::getDescription() const {
  return {R"(
LLCS Algorithm for n strings with arbitrary gap constraint tuple
/** Pseudocode ************************************************************
 *  LLCS2_MC(s1,s2,gc):
 *      let li = length(si) and li<=lj for all i<j
 *
 *      // Initialize a k-D table to store lengths of LCS for subproblems
 *      let M_1 = BinaryMatrix[lk + 1][lk + 1] such that M_1[i][j] := (s1[i]==s2[j])
 *      let M_k = Zero[lk + 1][lk + 1] filled with zeros for all k in [1..l1]
 *      let k = 1 if there is a shared symbol in s1 and s2
 *      if k == 0 return 0
 *
 *      // Fill M_p using dynamic programming and a window sliding approach. M_p is computed in O(l1*l2)
 *      for p in [2..l1]
 *          for i1 = 1 to l1 do
 *              for i2 = 1 to l2 do
 *                  if there is a one in M_{p-1}[I][J] with I=[i-u-1 : i-l-1] and J=[j-u-1 : j-l-1] then M_p[i][j] == 1
 *          if M_p[i][j] was set to 1 for some x and y
 *              set k to p
 *          else
 *              return k
 *      return k
 */
)"};
}

/*******************************************************************************
 * @brief Executes the query operation for the algorithm and returns a solution.
 *
 * @details This function resets the object's state, verifies its validity, and
 * then handles special cases related to the `gc-tuple`. The function proceeds
 * with preprocessing if no special case is encountered. After that, if the
 * object's state is set to `Preprocessed`, it returns a pointer to the solution
 *
 * @return std::unique_ptr<BaseSolution> A unique pointer to a solution object:
 * - `EmptySolution` if the object is not valid or preprocessing fails.
 * - `UnsignedSolution` with 0 or 1 if special cases for empty gc-tuples are met
 * - `UnsignedSolution` with the value of `k` if preprocessing is successful
 ******************************************************************************/
std::unique_ptr<BaseSolution> LLCS2_MC::query() {
  reset(ResetLevel::Full);
  if (!isValid()) {
    return std::make_unique<solutions::EmptySolution>();
  }

  // Handle special cases in which the gc-tuple is empty. (LLCS is 0 or 1)
  if (C.empty()) {
    if (s[0].at(0) == s[1].at(0)) {
      return std::make_unique<solutions::UnsignedSolution>(1);
    }
    return std::make_unique<solutions::UnsignedSolution>(0);
  }

  doPreprocessing();
  if (getState() != State::Preprocessed) {
    return std::make_unique<solutions::EmptySolution>();
  }
  return std::make_unique<solutions::UnsignedSolution>(k);
}

/*******************************************************************************
 * @brief Resets the state of the LLCS2_MC object to its initial condition.
 *
 * @details This function clears the matrices and vectors used by the algorithm
 * (`M`, `A`, and `B`), and sets the internal counters  to zero. After resetting,
 * the object's state is set to `Constructed`.
 *
 * @param lvl The level of reset to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_MC::reset(BaseAlgorithm::ResetLevel lvl) {
  M.clear();
  A.clear();
  B.clear();
  k = 0;
  n = 0;
  m = 0;
  setState(State::Constructed);
}

/*******************************************************************************
 * @brief Getter for the last element of the dp-Matrix vector M
 * @return Matrix M.back()
 * @note Operation is currently not used in a meaningful way
 *******************************************************************************/
const LLCS2_Algorithm::Matrix &LLCS2_MC::getMatrix() const {
  return M.back();
}

/*******************************************************************************
 * @brief Getter for the dp-Matrix vector M
 * @return Matrices M
 ******************************************************************************/
const LLCS2_MC::Matrices &LLCS2_MC::getMatrices() const {
  return M;
}

/*******************************************************************************
 * @brief Sets up the Dynamic Programming Algorithm
 *
 * @details Initializes M,A,B, m, n and k. It also calculates M[1]:
 * M[1][x][y] = 1, iff v[i - 1] == w[j - 1]. The LLCS placeholder k is set to 1
 * if there is a shared symbol in s[0] and s[1] otherwise it is zero
 ******************************************************************************/
void LLCS2_MC::doSetupDynProgramming() {
  auto v = s[0];
  auto w = s[1];
  m = v.size();// |s[0]|
  n = w.size();// |s[1]|
  k = 0;

  // Treat A[i][j], B[i][j], dp[i][j] as 0 if i<1 or j<1. Thus, add a dimension to arrays.
  M = Matrices(2);// M[0] is empty.
  M.reserve(m + 1);
  A = Matrix(m + 1, Row(n + 1, 0));
  B = Matrix(m + 1, Row(n + 1, 0));

  // M_1[i, j] = 1 iff there is a match at pos (i,j) that is v[i-1]==w[j-1]
  M[1] = Matrix(m + 1, Row(n + 1, 0));
  for (uint i = 1; i <= m; ++i) {
    for (uint j = 1; j <= n; ++j) {
      if (v[i - 1] == w[j - 1]) {
        M[1][i][j] = 1;
        k = 1;// llcs >= 1
        if (trackKeyPairs)
          track({i, j}, k);
      } else {
        M[1][i][j] = 0;
      }
    }
  }
}

/*******************************************************************************
 * @brief updates A
 * @details Sets matrix `A` such that `M[p]` can be calculated from `A` and `B`
 * @param p uint - the index in the current current loop pass
 *******************************************************************************/
void LLCS2_MC::updateA(uint p) {
  const uint l = C[p - 2].first;// When calculating M[2], the gap must fulfill C[0]
  const uint u = C[p - 2].second;
  const uint d = u - l + 1;

  // A stores the sum of d consecutive entries M_{p-1}[i][j-d+1 : j]
  for (uint i = 1; i <= m; ++i) {
    A[i][1] = M[p - 1][i][1];
    for (uint j = 2; j <= n; ++j)
      A[i][j] = A[i][j - 1] - M[p - 1][i][(j > d) ? j - d : 0] + M[p - 1][i][j];
  }
}

/*******************************************************************************
 * @brief updates B
 * @details Sets matrix `A` such that `M[p]` can be calculated from `A` and `B`
 * @param p index in the current loop pass
 *******************************************************************************/
void LLCS2_MC::updateB(uint p) {
  const uint l = C[p - 2].first;
  const uint u = C[p - 2].second;
  const uint d = u - l + 1;

  // B[i][j] stores the sum of all entries M_{p-1}[i-d+1:i][j-d+1:j]
  for (uint j = 1; j <= n; ++j) {
    B[1][j] = A[1][j];
  }
  for (uint i = 2; i <= m; ++i) {
    B[i][1] = B[i - 1][1] - M[p - 1][(i > d) ? i - d : 0][1] + M[p - 1][i][1];
  }
  for (uint i = 2; i <= m; ++i) {
    for (uint j = 2; j <= n; ++j) {
      B[i][j] = B[i - 1][j] - A[(i > d) ? i - d : 0][j] + A[i][j];
    }
  }
}

/*******************************************************************************
 * @brief updates M
 * @details Sets M[p] based on matrices A and B.
 * @param p index in the current current loop pass
 * @return true iff a subsequence has been extended
 ******************************************************************************/
bool LLCS2_MC::updateM(uint p) {
  const uint l = C[p - 2].first;
  //  const uint u = C[p - 2].second;
  //  const uint d = u - l + 1;
  /**
    * Update dp: dp[i][j] = M_p[i][j]
    * Invariant: M_p[i][j] == 1   <=> There is a one in M_{p-1}[I][J] with I=[i-u-1 : i-l-1] and J=[j-u-1 : j-l-1]
    *                              <=> B[i-l-1][j-l-1] > 0
    */
  bool continueFill = false;
  M.emplace_back(m + 1, Row(n + 1, 0));
  for (uint i = 1; i <= m; ++i) {
    for (uint j = 1; j <= n; ++j) {
      // M_p[i][j] == 1)  <==>    There is a 1 in M_{p-1}[i-u-1:i-l-1][j-u-1:j-l-1] and the lcs can be extended
      //                  <==>    B[i-l-1][j-l-1] > 0 && v[i] == w[j]
      //                  <==>    B[i-l-1][j-l-1] > 0 && M_1[i][j] >0
      uint const x = (i > l + 1) ? i - l - 1 : 0;
      uint const y = (j > l + 1) ? j - l - 1 : 0;
      if ((B[x][y] > 0) && M[1][i][j] == 1) {
        M[p][i][j] = 1;
        continueFill = true;
        if (trackKeyPairs) {
          track({i, j}, p);
        }
      } else {
        M[p][i][j] = 0;
      }
    }
  }
  return continueFill;
}

/*******************************************************************************
 * @brief does the preprocessing step of the algorithm
 ******************************************************************************/
void LLCS2_MC::doPreprocessing() {
  doSetupDynProgramming();
  if (k == 0) { return; }// no matching Symbols were found in the setup

  // Fill M_{p}[i, j] and assume M{p-1}[i,j] was already computed
  assert(C.size() + 1 == m && m <= n);
  for (uint p = 2; p <= m; ++p) {
    updateA(p);// A stores the sum of d consecutive entries M_{p-1}[i][j-d+1 : j]
    updateB(p);// B[i][j] stores the sum of all entries M_{p-1}[i-d+1:i][j-d+1:j]
    if (updateM(p)) {
      k = p;// and continue the loop
    } else {
      setState(State::Preprocessed);
      return;
    }
  }
  k = m;
  setState(State::Preprocessed);
}

/*******************************************************************************
 * @brief creates a std::string to describe the state of the algorithm
 * @return std::string containing:
 * - the individual strings in the algorithm
 * - the placeholder of the llcs `k`
 * - the dp matrices
 ******************************************************************************/
std::string LLCS2_MC::DebugString() const {
  std::ostringstream oss;
  oss << BaseAlgorithm::toString(s) << "\n";
  uint k = 0;
  for (const auto &mat : M) {
    std::string name = "M_" + std::to_string(k++);
    oss << toString(mat, s, true, name, true) << "\n";
  }
  oss << "llcs:" << k << "\n";
  oss << toString(keyPairs, trackKeyPairs) << "\n";
  return oss.str();
}

}// namespace lcs_solver::algorithms::llcs