/******************************************************************************
 * @file LLCS2_SA_RMQ.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Implementation of an algorithm for LCS2_SR with the Log-RMQ-DS
 * @details Time Complexity: \f$\mathcal O(n*n*sigma*\log(n))\f$
 *          Space Complexity:\f$\mathcal O(n*n*sigma*\log(n))\f$
 *****************************************************************************/
#include "algorithms/LLCS/LLCS2_SA_RMQ.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <ranges>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"

namespace lcs_solver::algorithms::llcs {
/*******************************************************************************
 * Constructor for LLCS2_SA_RMQ
 * @param vec The vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 ******************************************************************************/
LLCS2_SA_RMQ::LLCS2_SA_RMQ(const StrPtrVector &vec, const ConstraintMap &map)
    : LLCS2_SIG_Algorithm(AlgoType::LLCS2_SA_MQ, vec, map,
                          LLCS2_SIG_Algorithm::getSigLMap(
                              map,
                              ConstraintType::SIGMA
                              ),
                          LLCS2_SIG_Algorithm::getSigRMap(
                              map,
                              ConstraintType::SIGMA
                              )),
      v(s.size() == 2 ? s[0] : util::StringView()),
      w(s.size() == 2 ? s[1] : util::StringView()),
      m(s.size() == 2 ? s[0].size() : 0),
      n(s.size() == 2 ? s[1].size() : 0),
      maxq(n == 0 ? 1 : std::max<uint>(1, static_cast<uint>(1 + ceil(std::log2(static_cast<double>(n)))))), // n>=m
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
bool LLCS2_SA_RMQ::isValid() const {
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
std::string_view LLCS2_SA_RMQ::getDescription() const {
  return {R"(
LLCS Algorithm for the gc-constraint function right(a). Builds a 2D RMQ.
/** Pseudocode (Adamson et. al, LCS wit Gap constraints page 17-19) ************
 *  LLCS2_SR_RMQ(s1,s2,gc):
 *      // Initialize the relevant data structures
 *      set li = length(si) with li<=lj for all i<j
 *      set M   = Zeros[lk + 1, lk + 1]
 *      set RMQ_b = Zeros[lk + 1, lk + 1, ceil(log n)]
 *      set init gapRL[a,b] = (max {left[b].l, right[a].l}, min {left[b].u, right[a].u}) for all a,b in Sigma
 *      if s1[1] == s2[j] set M[1,j] and RMQ[1,j,0] to 1 for all j in [1:l2]
 *      if s1[i] == s2[1] set M[j,1] and RMQ[i,1,0] to 1 for all i in [1:l1]
 *      for q = 1 to ceil(log l1)
 *          set RMQ[1,j,q] = max {RMQ[1,j,q-1], RMQ[1, j-2^(q-1), q-1]}
 *          set RMQ[i,1,q] = max {RMQ[i,1,q-1], RMQ[i-2^(q-1), 1, q-1]}
 *
 *      // DP Loop
 *      for i in [2..l1]
 *          for j in [2..l2]
 *              if s1[i] = s2[v]  // Calculate M[i,j]
 *                  get the constraint right(a) = (la, ua) for a = s1[i]
 *                  set m = 0,
 *                  for b in Sigma
 *                      set I = [i-u-1:i-l-1], J_rp = [j-u-1:j-l-1] where (u,l) = gaRL[a,b]
 *                      set m = max {m, RMQ_b.query(I,J)}
 *                  set M = M_a[i,j] = m for all a in Sigma
 *
 *              // Always update RMQ data structure
 *              for q = 1 to ceil(log l1)
 *                  set l = 2^(q-1)
 *                  set RMQ[i,j,q] = max{RMQ[i-l,j-l,q-1], RMQ[i-l,j-l,q-1], RMQ[i,j-l,q-1], RMQ[i,j,q-1] }
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
std::unique_ptr<BaseSolution> LLCS2_SA_RMQ::query() {
  reset(ResetLevel::Full);
  if (!isValid()) {
    return std::make_unique<solutions::EmptySolution>();
  }
  doPreprocessing();
  return std::make_unique<solutions::UnsignedSolution>(max(M));
}

/*******************************************************************************
 * @brief Executes the main dp loop of the algorithm
 * @details Builds a RMQ data structure over squares with length that are powers
 * of two.
 ******************************************************************************/
void LLCS2_SA_RMQ::doPreprocessing() {
  auto gapRL = setup();
  fillFirstRow();
  fillFirstColumn();

  // Dynamic programming algorithm
  for (uint i = 2; i <= m; ++i) {
    for (uint j = 2; j <= n; ++j) {
      if (v[i - 1] == w[j - 1]) {
        // Update length of longest (left,right)-subsequence in M[a]
        const Symbol a = w[j - 1]; // `a` is to the right of the gap
        uint llcsCandidate = 0;
        for (const auto &[b, gap] : gapRL[a]) { // `b` is to the left of the gap
          if (gap.first > gap.second)
            continue;
          const auto [s0, e0, s1, e1] = getWindowOffSet(i, j, gap);
          const uint res = rmqQuery(RMQp[b], s0, e0, s1, e1);
          llcsCandidate = std::max<uint>(llcsCandidate, res);
          if (trackKeyPairs) track({i, j}, res + 1);
        }
        M[i][j] = llcsCandidate + 1;
        Mp[a][i][j] = M[i][j];
      }

      rmqUpdate(i, j); // Always Update rmq
    } /*j*/
  } /*i*/

  setState(State::Preprocessed);
}

/*******************************************************************************
 * @brief Does the setup for the dp main loop of the algorithm
 * @note side effect: initializes R as (m+1) x (n+1) zero matrix
 * @note side effect: initializes M[c] as a (m+1) x (n+1) x maxq zero-matrix for every c
 * @note side effect: calculates the gapIntersection
 * @return SigmaSigmaGapMap gap bounds for each symbols before/after the gap
 * SigmaSigmaGapMap[a][b] will be the bounds on a gap between to the left of
 * the gap and b to the right of the gap
 ******************************************************************************/
LLCS2_SA_RMQ::SigmaSigmaGapMap LLCS2_SA_RMQ::setup() {
  M = Matrix2D(m + 1, Row(n + 1, 0));
  // Calculate Gap Intersection
  std::unordered_map<Symbol, std::unordered_map<Symbol, Pair>> gapRL;
  for (const auto &b : std::views::keys(left)) {
    Mp[b] = Matrix2D(m + 1, Row(n + 1, 0));
    RMQp[b] = Matrix3D(m + 1, Matrix2D(n + 1, Row(maxq, 0)));
    gapRL.emplace(b, std::unordered_map<Symbol, std::pair<uint, uint>>());
  }
  for (const auto &[b, gapL] : left) {
    for (const auto &[a, gapR] : right) {
      auto gap = std::make_pair(
          std::max<uint>(gapL.first, gapR.first),
          std::min<uint>(gapL.second, gapR.second)
      );
      std::unordered_map<Symbol, std::pair<uint, uint>> const x;
      gapRL.at(a).emplace(b, gap);
    }
  }
  return gapRL;
}

/*******************************************************************************
 * @brief fills the first row of `M`, `Mp[a]`, `RMQ[a][q]`(for all a and q)
 * @details Fills the 1st row of RMQ, Mp and M (in O( sigma*m*N log m))
 ******************************************************************************/
void LLCS2_SA_RMQ::fillFirstRow() {
  for (uint j = 1; j <= n; ++j) {
    Symbol const a = w[j - 1];
    if (v[0] == a) {
      M[1][j] = 1;
      if (trackKeyPairs) track({1, j}, M[1][j]);
      Mp[a][1][j] = 1;
      RMQp[a][1][j][0] = 1;
    }
    for (auto &rmq : std::views::values(RMQp)) {
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
      } // end q
    } // end b
  } // end j
}

/*******************************************************************************
 * @brief fills the first Column of `M`, `Mp[a]`, `RMQ[a][q]`(for all a and q)
 * @details Fills the 1st column of RMQ, Mp and M (in O( sigma*m*N log m))
 ******************************************************************************/
void LLCS2_SA_RMQ::fillFirstColumn() {
  for (uint i = 1; i <= m; ++i) {
    Symbol const a = v[i - 1];
    if (a == w[0]) {
      M[i][1] = 1;
      if (trackKeyPairs) track({i, 1}, M[i][1]);
      Mp[a][i][1] = 1;
      RMQp[a][i][1][0] = Mp[a][i][1];
    }
    for (auto &rmq : std::views::values(RMQp)) {
      for (uint q = 1; q < maxq; ++q) {
        if (i >= POW2[q - 1]) {
          rmq[i][1][q] = std::max<uint>(
              rmq[i][1][q - 1],
              rmq[i - POW2[q - 1]][1][q - 1]
          );
        } else {
          rmq[i][1][q] = rmq[i][1][q - 1];
        }
      } // end q
    } // end b
  } // end  i
}

/*******************************************************************************
 * @brief RMQ queries the area [s0:e0][s1:e1]
 * @param mat RMQ Matrix, where the last axis of mat stores q for log intervals
 * @param s0 first index in the window along axis 0
 * @param e0 last index in the window along axis 0
 * @param s1 first index in the window along axis 1
 * @param e1 last index in the window along axis 1
 * @return Maximum in RMQ[s0:e0]x[s1:e1]
 ******************************************************************************/
LLCS2_SA_RMQ::uint LLCS2_SA_RMQ::rmqQuery(const Matrix3D &mat, uint s0, uint e0, uint s1, uint e1) const {
  assert(e0 >= s0 && e1 >= s1);
  uint const len0 = e0 - s0 + 1;
  uint const len1 = e1 - s1 + 1;
  uint const len = std::max<uint>(len0, len1); // square length
  uint const q = LOG2[len]; // a square with length 2^q will into [s0:e0]x[s1:e1]
  uint result = mat[e0][e1][q];
  if (len != POW2[q] || len0 != len1) {
    uint m1 = 0;
    uint m2 = 0;
    uint m3 = 0;
    uint const m4 = mat[e0][e1][q];
    if (s0 + POW2[q] <= e0 && s1 + POW2[q] <= e1)
      m1 = mat[s0 + POW2[q] - 1][s1 + POW2[q] - 1][q];
    if (s0 + POW2[q] <= e0)
      m2 = mat[s0 + POW2[q] - 1][e1][q];
    if (s1 + POW2[q] <= e1)
      m3 = mat[e0][s1 + POW2[q] - 1][q];
    result = std::max({m1, m2, m3, m4});
  }
  return result;
}

/*******************************************************************************
 * Resets the state of the LLCS2_SA_RMQ object to its initial condition.
 * @note side effect: Erases the elements from M and frees memory
 * @param lvl BaseAlgorithm::ResetLevel to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_SA_RMQ::reset(BaseAlgorithm::ResetLevel lvl) {
  for (auto &mat : Mp | std::views::values) {
    mat.clear();
    mat.shrink_to_fit();
  }
  Mp.clear();

  for (auto &rmq : RMQp | std::views::values) {
    rmq.clear();
    rmq.shrink_to_fit();
  }
  RMQp.clear();

  M.clear();
  setState(State::Constructed);
}

/*******************************************************************************
 * @brief Getter for the dp-Matrix
 * @return std::vector<std::vector<uint> M
 *******************************************************************************/
const LLCS2_Algorithm::Matrix &LLCS2_SA_RMQ::getMatrix() const {
  return M;
}

/*******************************************************************************
 * isExtensible
 * @param a Pair (a1,a2) a one based position: s1[a1],s2[a2]
 * @param b Pair (b1,b2) a one based position: s1[b1],s2[b2]
 * @param llcsOfA Uint equal to LLCS(s1[1:a1],s2[1:a2])
 * @return Whether an LCS from position a can be extended to position b
 ******************************************************************************/
bool LLCS2_SA_RMQ::isExtensible(Pair a, Pair b, uint llcsOfA) const {
  return LLCS2_SIG_Algorithm::isExtensibleHelper(a,b,false, false);
}

/*******************************************************************************
 * getPrevRange
 * @param pair Pair (x,y) a one based position: s1[x],s2[y]
 * @param llcs The length of the LCS at position (x,y)
 * @return (s1Factor, s2Factor) where a factor is a one based uint pair. It
 *  contains all possible positions where a position (i,j) can be
 *  extended to (x,y).
 * @note The returned Window does not consider the left constraint, because for
 *  it one must know the Symbol to the left of a gap.
 ******************************************************************************/
LLCS2_Algorithm::Window LLCS2_SA_RMQ::getPrevRange(const Pair &pair, uint llcs) const {
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
 * - the dp matrices M_b (for all b in Sigma)
 * - the dp matrix M
 ******************************************************************************/
std::string LLCS2_SA_RMQ::DebugString() const {
  std::ostringstream oss;
//  oss << "===================================================================\n";
//  oss << "rmq is: " << (RMQp.empty() ? "empty \n" : "\n");
//  for (const auto &[b, rmq] : RMQp) {
//    for (uint q = 0; q < maxq; ++q) {
//      oss << "rmq[b ==" << b << " q == " << q << "] is:\n ";
//      for (uint i = 1; i <= m; ++i) {
//        for (uint j = 1; j <= n; ++j) {
//          if (!rmq.empty() && !rmq[i].empty() && rmq[i][j].size() >= q)
//            oss << rmq[i][j][q] << " ";
//        }
//        oss << "\n ";
//      }
//      oss << "\n";
//    }
//  }
  oss << "M is: " << (M.empty() ? "empty \n" : "\n");
  for (uint i = 1; i <= m; ++i) {
    for (uint j = 1; j <= n; ++j) {
      oss << M[i][j] << " ";
    }
    oss << "\n";
  }

  oss << "\n";
  if (Mp.empty())
    oss << "Mp is empty \n";
  for (const auto &[b, mat] : Mp) {
    oss << "M[symbol == " << util::to_string(b) << "] is: \n";
    for (const auto &row : mat) {
      for (const auto &val : row) {
        oss << val << " ";
      }
      oss << "\n";
    }
    oss << "\n";
  }
  return oss.str();
}

/*******************************************************************************
 * @brief Updates the `rmq` object of the algorithm
 * @param i uint current row when updating RMQ[i,j,q]
 * @param j uint current colum when updating RMQ[i,j,q]
 ******************************************************************************/
void LLCS2_SA_RMQ::rmqUpdate(uint i, uint j) {
  for (auto &[b, rmq] : RMQp) {
    rmq[i][j][0] = Mp[b][i][j];
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
  } /*b*/
}

/*******************************************************************************
 * @brief Determines the relevant area to look in for extending a LCS when left
 * and right map constraints are active
 * @param i uint current row when calculating M[i,j]
 * @param j uint current colum when calculating M[i,j]
 * @param gap Pair for calculating the right window size
 * @return [s0,e0, s1, e1] defines the area [s0:e1, s1:e1]
 ******************************************************************************/
LLCS2_SA_RMQ::WindowBounds LLCS2_SA_RMQ::getWindowOffSet(
    const uint i,
    const uint j,
    LLCS2_SIG_Algorithm::Pair gap) {
  auto [l, u] = gap;
  const uint s0 = i > u + 1 ? i - u - 1 : 0;
  const uint e0 = i > l + 1 ? i - l - 1 : 0;
  const uint s1 = j > u + 1 ? j - u - 1 : 0;
  const uint e1 = j > l + 1 ? j - l - 1 : 0;
  return {s0, e0, s1, e1};
}

} // end of namespace
