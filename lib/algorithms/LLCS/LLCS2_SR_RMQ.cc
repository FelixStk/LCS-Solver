/******************************************************************************
 * @file LLCS2_SR_RMQ.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Implementation of an algorithm for LCS2_SR that uses RMQ
 * @details Time Complexity: O(n*m*log(n))
 *****************************************************************************/
#include "algorithms/LLCS/LLCS2_SR_RMQ.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"

namespace {
using ::lcs_solver::algorithms::solutions::EmptySolution;
using ::lcs_solver::algorithms::solutions::UnsignedSolution;
using ::lcs_solver::util::Log2_table;
using ::lcs_solver::util::Pow2_table;
using uint = ::lcs_solver::algorithms::llcs::LLCS2_SR_RMQ::uint;
}

namespace lcs_solver::algorithms::llcs {
/*******************************************************************************
 * @brief Constructor for LLCS2_SR_RMQ
 * @param vec The vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 ******************************************************************************/
LLCS2_SR_RMQ::LLCS2_SR_RMQ(const StrPtrVector &vec, const ConstraintMap &map)
    : LLCS2_SIG_Algorithm(AlgoType::LLCS2_SR_MQ, vec, map,
                          LLCS2_SIG_Algorithm::getSigLMap(
                              map,
                              ConstraintType::SIGMA_R),
                          LLCS2_SIG_Algorithm::getSigRMap(
                              map,
                              ConstraintType::SIGMA_R)),
      v(s.size() == 2 ? s[0] : util::StringView()),
      w(s.size() == 2 ? s[1] : util::StringView()),
      m(s.size() == 2 ? s[0].size() : 0),
      n(s.size() == 2 ? s[1].size() : 0),
      maxq(n == 0 ? 1 : std::max<uint>(1, static_cast<uint>(1 + ceil(std::log2((double) n))))), // n>=m
      LOG2(n + 1),
      POW2(maxq + 1) {

}

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS2_SR_RMQ::isValid() const {
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
std::string_view LLCS2_SR_RMQ::getDescription() const {
  return {R"(
LLCS Algorithm for the gc-constraint function right(a). Builds a 2D RMQ.
/** Pseudocode (Adamson et. al, LCS wit Gap constraints page 17-19) ************
 *  LLCS2_SR_RMQ(s1,s2,gc):
 *      // Initialize the relevant data structures
 *      set li = length(si) with li<=lj for all i<j
 *      set M   = Zeros[lk + 1, lk + 1]
 *      set RMQ = Zeros[lk + 1, lk + 1, ceil(log n)]
 *      if s1[1] == s2[j] set M[1,j] and RMQ[1,j,0] to 1 for all j in [1:l2]
 *      if s1[i] == s2[1] set M[j,1] and RMQ[i,1,0] to 1 for all i in [1:l1]
 *      for q = 1 to ceil(log l1)
 *          set RMQ[1,j,q] = max {RMQ[1,j,q-1], RMQ[1, j-2^(q-1), q-1]}
 *          set RMQ[i,1,q] = max {RMQ[i,1,q-1], RMQ[i-2^(q-1), 1, q-1]}
 *
 *      // DP Loop
 *      for i in [2..l1]
 *          for j in [2..l2]
 *              get the constraint right(a) = (la, ua) for a = s1[i] = s2[v];
 *              set I = [i-ua-1:i-la-1] and J_rp = [j-ua-1:j-la-1];
 *              set M[i][j] as the maximum of M[I][J], computed using RMQ;
 *              set RMQ[i,j,0] = M[i][j]
 *              for q = 1 to ceil(log l1)
 *                  set l = 2^(q-1)
 *                  set RMQ[i,j,q] = max{RMQ[i-l,j-l,q-1], RMQ[i-l,j-l,q-1], RMQ[i,j-l,q-1], RMQ[i,j,q-1] }
 *              set M_a to be 0 when Ia or Ja are empty
 *              set M[i][j] = 1 + M_a
 *      return max {M[x][y] | x in [1..l1] and y in [1..2] }
 */
)"};
}

/*******************************************************************************
 * @brief Executes the query operation for the algorithm and returns a solution.
 * @details This function resets the object's state and verifies its validity.
 * The function proceeds with preprocessing and returns a solution pointer.
 * @return std::unique_ptr<BaseSolution> A unique pointer to a solution object:
 * - `EmptySolution` if the object is not valid
 * - `UnsignedSolution` containing the llcs calculated
 ******************************************************************************/
std::unique_ptr<BaseSolution> LLCS2_SR_RMQ::query() {
  reset(ResetLevel::Full);
  if (!isValid()) { return std::make_unique<EmptySolution>(); }
  doPreprocessing();
  return std::make_unique<UnsignedSolution>(max(M));
//  return std::make_unique<UnsignedSolution>(M[m][n]);
}

/*******************************************************************************
 * @brief Executes the dp algorithm by Adamson et al. and updates status
 ******************************************************************************/
void LLCS2_SR_RMQ::doPreprocessing() {
  M = Matrix2D(m + 1, Row(n + 1, 0));
  rmq = Matrix3D(m + 1, Matrix2D(n + 1, Row(maxq, 0)));
  setFirstRowColum();

  // Dynamic programming algorithm
  for (uint i = 2; i <= m; ++i) {
    for (uint j = 2; j <= n; ++j) {
      if (v[i - 1] == w[j - 1]) {
        // Get the constraint right(a) = (l,u) for a=v[i-1]==w[j-1]
        const Symbol a = w[j - 1];
        const uint l = right.at(a).first;
        const uint u = right.at(a).second;

        // Set M[i][j] to the maximum in the square sub-matrix M[s0:e0][s1:e1]
        const uint s0 = i > u + 1 ? i - u - 1 : 0;
        const uint e0 = i > l + 1 ? i - l - 1 : 0;
        const uint s1 = j > u + 1 ? j - u - 1 : 0;
        const uint e1 = j > l + 1 ? j - l - 1 : 0;
        const uint llcsCandidate = rmqQuery(s0, e0, s1, e1);
        M[i][j] = llcsCandidate + 1;
        if (trackKeyPairs) track({i, j}, M[i][j]);
//        M[i][j] = std::max<uint>({llcsCandidate, M[i - 1][j], M[i][j - 1]});
      } else {
//        M[i][j] = std::max<uint>(M[i - 1][j], M[i][j - 1]);
      }
      rmqUpdate(i, j); // Always Update rmq
    } /*j*/
  } /*i*/
  setState(State::Preprocessed);
}

/*******************************************************************************
 * @brief Fills the first rows and columns of `M`, `rmq` with correct values
 * @details This means:
 * - M[1][j] == 1 <=> Phi(1,j)
 * - M[i][1] == 1 <=> Phi(i,1)
 * - RMQ[1][j][q] :=  max M[1][x] with x in Xq := [j-2^q+1 : j] with |Xq| = 2^q
 * - RMQ[i][1][q] :=  max M[x][1] with x in Xq := [j-2^q+1 : j] with |Xq| = 2^q
 ******************************************************************************/
void LLCS2_SR_RMQ::setFirstRowColum() {
  // Fill 1st row of RMQ and M (in O( m*N log m))
  for (uint j = 1; j <= n; ++j) {
    if (v[0] == w[j - 1]) {
      M[1][j] = 1;
      if (trackKeyPairs) track({1, j}, M[1][j]);
      rmq[1][j][0] = 1;
    }
    // Update rmq: rmq[1][j][q] = max M[1][x] with x in Xq := [j-2^q+1 : j]. Note |Xq| = 2^q
    for (uint q = 1; q < maxq; ++q) {
      if (j >= POW2[q - 1]) {
        // Recursion is well-defined
        rmq[1][j][q] = std::max<uint>(
            rmq[1][j][q - 1],               // max M[1][j - 2^{q-1} + 1 : j]
            rmq[1][j - POW2[q - 1]][q - 1]  // max M[1][j - 2^q + 1 : j - 2^{q-1}]
        );
      } else {
        // Copy M[1][2^{q-1} + 1 : j]. The value invalid ranges (e.g M[1][-1:2])
        rmq[1][j][q] = rmq[1][j][q - 1];
      }
    }
  } // end 1st row setup

  // Fill 1st column of RMQ and M (in O( m*N log m))
  for (uint i = 1; i <= m; ++i) {
    if (v[i - 1] == w[0]) {
      M[i][1] = 1;
      if (trackKeyPairs) track({i, 1}, M[i][1]);
      rmq[i][1][0] = M[i][1];
    }
    for (uint q = 1; q < maxq; ++q) {
      if (i >= POW2[q - 1]) {
        rmq[i][1][q] = std::max<uint>(
            rmq[i][1][q - 1],
            rmq[i - POW2[q - 1]][1][q - 1]
        );
      } else {
        rmq[i][1][q] = rmq[i][1][q - 1];
      }
    } // end for q
  } // end first column setup
}

/*******************************************************************************
 * @brief Updates the `rmq` object of the algorithm
 * @param i uint current row when updating RMQ[i,j,q]
 * @param j uint current colum when updating RMQ[i,j,q]
 ******************************************************************************/
void LLCS2_SR_RMQ::rmqUpdate(uint i, uint j) {
  rmq[i][j][0] = M[i][j];
  for (uint q = 1; q < maxq; ++q) {
    uint m1 = 0;
    uint m2 = 0;
    uint m3 = 0;
    uint m4 = 0;
    if (i >= POW2[q - 1] && j >= POW2[q - 1])
      m1 = rmq[i - POW2[q - 1]][j - POW2[q - 1]][q - 1];
    if (j >= POW2[q - 1])
      m2 = rmq[i][j - POW2[q - 1]][q - 1];
    if (i >= POW2[q - 1])
      m3 = rmq[i - POW2[q - 1]][j][q - 1];
    m4 = rmq[i][j][q - 1];
    rmq[i][j][q] = std::max({m1, m2, m3, m4});
  } /*q*/
}

/*******************************************************************************
 * @brief RMQ queries the area [s0:e0][s1:e1]
 * @param s0 first index in the window along axis 0
 * @param e0 last index in the window along axis 0
 * @param s1 first index in the window along axis 1
 * @param e1 last index in the window along axis 1
 * @return Maximum in RMQ[s0:e0]x[s1:e1]
 ******************************************************************************/
LLCS2_SR_RMQ::uint LLCS2_SR_RMQ::rmqQuery(const uint s0, const uint e0, const uint s1, const uint e1) const {
  assert(e0 >= s0 && e1 >= s1);
  uint result = 0;
  uint const len0 = e0 - s0 + 1;
  uint const len1 = e1 - s1 + 1;
  uint const len = std::max<uint>(len0, len1); // square length
  uint const q = LOG2[len]; // a square with length 2^q will into [s0:e0]x[s1:e1]
  result = rmq[e0][e1][q];
  if (len != POW2[q] || len0 != len1) {
    uint m1 = 0;
    uint m2 = 0;
    uint m3 = 0;
    uint const m4 = rmq[e0][e1][q];
    if (s0 + POW2[q] <= e0 && s1 + POW2[q] <= e1)
      m1 = rmq[s0 + POW2[q] - 1][s1 + POW2[q] - 1][q];
    if (s0 + POW2[q] <= e0)
      m2 = rmq[s0 + POW2[q] - 1][e1][q];
    if (s1 + POW2[q] <= e1)
      m3 = rmq[e0][s1 + POW2[q] - 1][q];
    result = std::max({m1, m2, m3, m4});
  }
  return result;
}

/*******************************************************************************
 * Resets the state of the LLCS2_SR_RMQ object to its initial condition.
 * @note Side Effect: Erases the elements from M and frees memory
 * @param lvl BaseAlgorithm::ResetLevel to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_SR_RMQ::reset(BaseAlgorithm::ResetLevel lvl) {
  M.clear();
  rmq.clear();
  setState(State::Constructed);
}

/*******************************************************************************
 * @brief Getter for the dp-Matrix
 * @return std::vector<std::vector<uint> M
 *******************************************************************************/
const LLCS2_Algorithm::Matrix &LLCS2_SR_RMQ::getMatrix() const {
  return M;
}

/*******************************************************************************
 * isExtensible
 * @param a Pair (a1,a2) a one based position: s1[a1],s2[a2]
 * @param b Pair (b1,b2) a one based position: s1[b1],s2[b2]
 * @param llcsOfA Uint equal to LLCS(s1[1:a1],s2[1:a2])
 * @return Whether an LCS from position a can be extended to position b
 ******************************************************************************/
bool LLCS2_SR_RMQ::isExtensible(Pair a, Pair b, uint llcsOfA) const {
  return LLCS2_SIG_Algorithm::isExtensibleHelper(a,b,false, true);
}

/*******************************************************************************
 * getPrevRange
 * @param pair Pair (x,y) a one based position: s1[x],s2[y]
 * @param llcs The length of the LCS at position (x,y)
 * @return (s1Factor, s2Factor) where a factor is a one based uint pair. It
 *  contains all possible positions where a position (i,j) can be
 *  extended to (x,y).
 ******************************************************************************/
LLCS2_Algorithm::Window LLCS2_SR_RMQ::getPrevRange(const Pair &pair, uint llcs) const {
  assert(s.size() == 2);
  assert(llcs > 0);
  if (llcs == 1) {
    // The gap before the first position (with LLCS(x,y)==1) is not considered
    Pair r1 = {1, pair.first - 1};
    Pair r2 = {1, pair.second - 1};
    return {r1,r2};
  }
  Symbol c = s[0].at(pair.first-1);
  return getPrefRange(pair, c, right);
}

/*******************************************************************************
 * @brief Creates a std::string to describe the state of the algorithm
 * @return std::string containing:
 * - the rmq object
 * - the dp matrix M
 ******************************************************************************/
std::string LLCS2_SR_RMQ::DebugString() const {
  std::ostringstream oss;
  oss << "===================================================================\n";
  oss << "rmq is: " << (rmq.empty() ? "empty \n" : "\n");
  for (uint q = 0; q < maxq; ++q) {
    oss << "rmq[ q == " << q << "] is:\n";
    for (uint i = 1; i <= m; ++i) {
      for (uint j = 1; j <= n; ++j) {
        if (!rmq.empty() && !rmq[i].empty() && rmq[i][j].size() >= q)
          oss << rmq[i][j][q] << " ";
        else
          oss << "n ";
      }
      oss << "\n";
    }
  }
  oss << "M is: " << (M.empty() ? "empty \n" : "\n");
  for (uint i = 1; i <= m; ++i) {
    for (uint j = 1; j <= n; ++j) {
      oss << M[i][j] << " ";
    }
    oss << "\n";
  }
  return oss.str();
}

} // end of namespace