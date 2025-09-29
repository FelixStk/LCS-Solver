/******************************************************************************
 * @file LLCS2_SR_MQ.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Implementation of an algorithm for LCS2_SR that uses MaxQueue2D
 * @details Time Complexity: O(sigma * n * m)
 *****************************************************************************/
#include "algorithms/LLCS/LLCS2_SR_MQ.h"

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"
#include "structures/MaxQueueN.h"

namespace{
using ::lcs_solver::algorithms::solutions::EmptySolution;
using ::lcs_solver::algorithms::solutions::UnsignedSolution;
using lcs_solver::structures::MaxQueue2D;
}

namespace lcs_solver::algorithms::llcs {
/*******************************************************************************
 * Constructor for LLCS2_SR_MQ objects
 * @param vec The vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 ******************************************************************************/
LLCS2_SR_MQ::LLCS2_SR_MQ(const StrPtrVector &vec, const ConstraintMap &map)
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
      n(s.size() == 2 ? s[1].size() : 0) { }

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS2_SR_MQ::isValid() const {
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
std::string_view LLCS2_SR_MQ::getDescription() const {
  return {R"(
LLCS Algorithm for the gc-constraint function right(a). Uses a MaxQueue2D.
/** Pseudocode (Adamson et. al, LCS wit Gap constraints page 16-17) ************
 *  LLCS2_SR_MQ(s1,s2,gc):
 *      // Initialize the relevant data structures
 *      let li = length(si) and li<=lj for all i<j
 *      init M, M_b with Zeros[lk + 1, lk + 1] for all b in Sigma
 *      let D_rp = MaxQueue2D(const MatrixPtr = &M, xRange = [0..l1], yRange[0..l2], bufferAxisSize = u-l) for rp in [h]
 *
 *      // Fill M_p using dynamic programming and a window sliding approach. M_p is computed in O(l1*l2)
 *      for i in [0..l1]
 *          update D_b for b in Sigma
 *          for j in [0..l2]
 *              update D_b for b in Sigma
 *              if s1[i] == s2[j], let a = s1[i], I_a = [i-ua-1:i-la-1] and J_rp = [j-ua-1:j-la-1]
 *              use D_a to retrieve M_a, the maximum of M[I_a, J_a]for b in Sigma
 *              set M_a to be 0 when Ia or Ja are empty
 *              set M[i,j] = 1 + M_a
 *      return max {M[x,y] | x in [1..l1] and y in [1..2] }
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
std::unique_ptr<BaseSolution> LLCS2_SR_MQ::query() {
  reset(ResetLevel::Full);
  if (!isValid())
    return std::make_unique<EmptySolution>();
  doPreprocessing();
  uint max = 0;
  for (const auto &row : M)
    for (const auto &val : row)
      max = std::max<uint>(max, val);
  return std::make_unique<UnsignedSolution>(max);
//  return std::make_unique<UnsignedSolution>(M[m][n]);
}

/*******************************************************************************
 * @brief Executes the dp algorithm by Adamson et al. and updates status
 ******************************************************************************/
void LLCS2_SR_MQ::doPreprocessing() {
  // Initialize Memory for Matrices
  M = Matrix(m + 1, Row(n + 1, 0));

  // Setup D['a'] with new MaxQueue2D
  std::unordered_map<Symbol, MaxQueue2D<uint >> D;
  for (const auto &[symbol, gap] : right) {
    const size_t yBuffer = std::max<uint>(1, gap.second - gap.first);
    D.emplace(symbol, MaxQueue2D<uint>(M, yBuffer, 0));
  }

  // Fill Matrix with dynamic programming approach
  for (uint i = 1; i <= m; ++i) {
    for (auto &[symbol, Qs] : D) {
      auto gap = right.at(symbol);
      const uint s0 = gap.second + 1;
      const uint e0 = gap.first + 1;
      Qs.doLineSetup(i, s0, e0);
    }
    for (uint j = 1; j <= n; ++j) {
      uint max = 0;
      for (auto &[symbol, Qs] : D) {
        auto gap = right.at(symbol);
        const uint s0 = gap.second + 1;
        const uint e0 = gap.first + 1;
        const uint s1 = gap.second + 1;
        const uint e1 = gap.first + 1;
        Qs.doComputingSetup(i, j, e0, e1); // Inserts M[i-e0][i-e1] into Qs
        if (Phi(i, j) && symbol == v[i - 1]) {
          uint res = Qs.query(i, j, {s0, e0, s1, e1});
          if (trackKeyPairs) track({i, j}, res + 1);
          max = std::max<uint>(max, res);
        }
      }
      if (Phi(i, j))
        M[i][j] = 1 + max;
      //else if (i > 0 && j > 0)
      //M[i][j] = std::max<uint>(M[i - 1][j], M[i][j - 1]);
    }
  }
  setState(State::Preprocessed);
}

/*******************************************************************************
 * Resets the state of the LLCS2_SR_MQ object to its initial condition.
 * @param lvl BaseAlgorithm::ResetLevel to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_SR_MQ::reset(BaseAlgorithm::ResetLevel lvl) {
  M.clear();
  setState(State::Constructed);
}

/*******************************************************************************
 * @brief Getter for the dp-Matrix
 * @return std::vector<std::vector<uint> M
 *******************************************************************************/
const LLCS2_Algorithm::Matrix &LLCS2_SR_MQ::getMatrix() const {
  return M;
}

/*******************************************************************************
 * isExtensible
 * @param a Pair (a1,a2) a one based position: s1[a1],s2[a2]
 * @param b Pair (b1,b2) a one based position: s1[b1],s2[b2]
 * @param llcsOfA Uint equal to LLCS(s1[1:a1],s2[1:a2])
 * @return Whether an LCS from position a can be extended to position b
 ******************************************************************************/
bool LLCS2_SR_MQ::isExtensible(Pair a, Pair b, uint llcsOfA) const {
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
LLCS2_Algorithm::Window LLCS2_SR_MQ::getPrevRange(const Pair &pair, uint llcs) const {
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
 * @return std::string containing the dp matrices `M`
 ******************************************************************************/
std::string LLCS2_SR_MQ::DebugString() const {
  std::ostringstream oss;
  oss << BaseAlgorithm::toString(s) << "\n\n";
  oss << toString(M, s, false, "M", true) << "\n";
  oss << toString(keyPairs, trackKeyPairs) << "\n";
  return oss.str();
}

/*******************************************************************************
 * @brief Checks for matching symbols in the strings
 * @note i and j are one based indices for positions in s_1 and s_2
 * @param i uint represents a one based position in the first problem string
 * @param j uint represents a one based position in the second problem string
 * @return true iff the problems strings have matching symbols at position i & j
 ******************************************************************************/
bool LLCS2_SR_MQ::Phi(uint i, uint j) const {
  if (i == 0 || i > v.size()) return false;
  if (j == 0 || j > w.size()) return false;
  return v[i - 1] == w[j - 1];
}

} // end of namespace