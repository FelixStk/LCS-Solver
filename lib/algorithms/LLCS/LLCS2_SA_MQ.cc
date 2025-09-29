/******************************************************************************
 * @file LLCS2_SA_MQ.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Implementation of an algorithm for LCS2_SR with the MaxQueue
 * @details Time Complexity: \f$\mathcal O (n*m*\sigma^2)\f$
 *****************************************************************************/
#include "algorithms/LLCS/LLCS2_SA_MQ.h"

#include <algorithm>
#include <memory>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::llcs {
/*******************************************************************************
 * Constructor for LLCS2_SA_MQ
 * @param vec vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 ******************************************************************************/
LLCS2_SA_MQ::LLCS2_SA_MQ(const StrPtrVector &vec, const ConstraintMap &map)
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
      n(s.size() == 2 ? s[1].size() : 0) { }

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS2_SA_MQ::isValid() const {
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
std::string_view LLCS2_SA_MQ::getDescription() const {
  return {R"(
LLCS Algorithm for additional left and right gap constraints. Uses a MaxQueue2D
/** Pseudocode (Adamson et. al, LCS wit Gap constraints page 19-20) ************
 *  LLCS2_SR_MQ(s1,s2,gc):
 *      // Initialize the relevant data structures
 *      let li = length(si) and li<=lj for all i<j
 *      init M, M_b with Zeros[lk + 1, lk + 1] for all b in Sigma
 *      init gapLR[b,a] = (max {left[b].l, right[a].l}, min {left[b].u, right[a].u}) for all a,b in Sigma
 *      init D_ba = MaxQueue2D for b,a in Sigma
 *
 *      // Fill M_p using dynamic programming and a window sliding approach. M_p is computed in O(l1*l2)
 *      for i in [0..l1]
 *          update D_ba for a,b in Sigma (lineSetups)
 *          for j in [0..l2]
 *              update D_ba for a,b in Sigma (computingSetup)
 *              set max = 0
 *              for a,b in Sigma x Sigma
 *                  let (l,u) = gapLR[b][a], I_ab = [i-u-1:i-l-1] and J_ab = [j-u-1:j-l-1]
 *                  if s1[i] == s2[j]
 *                      use D_ba to retrieve max, the maximum of M[I_ab, J_ab]
 *              if s1[i] == s2[j],
 *                  M[i][j] = M[a][i][j] = max + 1 where a = s1[i]
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
std::unique_ptr<BaseSolution> LLCS2_SA_MQ::query() {
  reset(ResetLevel::Full);
  if (!isValid())
    return std::make_unique<solutions::EmptySolution>();
  doPreprocessing();
  uint const llcs = getMaxMatrixElement(M);
  return std::make_unique<solutions::UnsignedSolution>(llcs);
}

/*******************************************************************************
 * @brief Executes the main dp loop of the algorithm
 * @details Calculates the gapIntersection Map (first key for left, second key
 * for the right map symbols), afterward is sets up MaxQueues to do finally the
 * dp algorithm to update the 2D matrices Mp[a] and M
 ******************************************************************************/
void LLCS2_SA_MQ::doPreprocessing() {
  M = Matrix(m + 1, Row(n + 1, 0));
  auto gapLR = calcGapIntersection();
  for (const auto &b : gapLR | std::views::keys) {
    Mp[b] = Matrix(m + 1, Row(n + 1, 0));
  }
  auto D = getMaxQueue(gapLR);

  // Fill Matrix with dynamic programming approach
  for (uint i = 1; i <= m; ++i) {
    doLineSetups(i, D, gapLR);
    for (uint j = 1; j <= n; ++j) {
      // ComputingSetup and get max extension length max
      uint max = 0;
      for (auto &[b, Db] : D) {
        for (auto &[a, Qba] : Db) {
          auto gap = gapLR[b][a];
          if (gap.first > gap.second)
            continue;
          const uint s0 = gap.second + 1;
          const uint e0 = gap.first + 1;
          const uint s1 = gap.second + 1;
          const uint e1 = gap.first + 1;
          Qba.doComputingSetup(i, j, e0, e1); // Inserts M[i-e0][i-e1] into Qs
          if (Phi(i, j) && a == v[i - 1]) {
            const uint res = Qba.query(i, j, {s0, e0, s1, e1});
            if (trackKeyPairs) track({i, j}, res + 1);
            max = std::max<uint>(max, res);
          }
        } /*a*/
      } /*b*/

      if (Phi(i, j)) {
        M[i][j] = max + 1;
        const Symbol a = v[i - 1];
        Mp[a][i][j] = max + 1;
      }

    } /*j*/
  } /*i*/
  setState(State::Preprocessed);
}

/*******************************************************************************
 * Resets the state of the LLCS2_SA_MQ object to its initial condition.
 * @note Side Effect: Erases the elements from M and frees memory
 * @param lvl BaseAlgorithm::ResetLevel to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_SA_MQ::reset(BaseAlgorithm::ResetLevel lvl) {
  for (auto &matrix : Mp | std::views::values) {
    matrix.clear();
    matrix.shrink_to_fit();
  }
  Mp.clear();
  M.clear();
  M.shrink_to_fit();
  setState(State::Constructed);
}
/*******************************************************************************
 * @brief Getter for the dp-Matrix
 * @return std::vector<std::vector<uint> M
 *******************************************************************************/
const LLCS2_Algorithm::Matrix &LLCS2_SA_MQ::getMatrix() const {
  return M;
}

/*******************************************************************************
 * isExtensible
 * @param a Pair (a1,a2) a one based position: s1[a1],s2[a2]
 * @param b Pair (b1,b2) a one based position: s1[b1],s2[b2]
 * @param llcsOfA Uint equal to LLCS(s1[1:a1],s2[1:a2])
 * @return Whether an LCS from position a can be extended to position b
 ******************************************************************************/
bool LLCS2_SA_MQ::isExtensible(Pair a, Pair b, uint llcsOfA) const {
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
LLCS2_Algorithm::Window LLCS2_SA_MQ::getPrevRange(const Pair &pair, uint llcs) const {
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
 * @brief Searches a matrix for its largest value
 * @param matrix std::vector<std::vector<uint>> represents a matrix
 * @return uint the larges value stored in the matrix
 ******************************************************************************/
LLCS2_SA_MQ::uint LLCS2_SA_MQ::getMaxMatrixElement(const Matrix &matrix) {
  auto max = matrix[0][0];
  for (const auto &row : matrix)
    for (const auto &elem : row)
      max = std::max(elem, max);
  return max;
}

/*******************************************************************************
 * @brief Creates a std::string to describe the state of the algorithm
 * @return std::string containing:
 * - the dp matrices M_b (for all b in Sigma)
 * - the dp matrix M
 ******************************************************************************/
std::string LLCS2_SA_MQ::DebugString() const {
  std::ostringstream oss;
  if (M.empty())
    oss << "M is empty \n";
  else {
    oss << "M is: \n";
    for (const auto &row : M) {
      for (const auto &val : row) {
        oss << val << " ";
      }
      oss << "\n";
    }
  }
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
  }
  return oss.str();
}

/*******************************************************************************
 * @brief Checks for matching symbols in the strings
 * @note i and j are one based indices for positions in s_1 and s_2
 * @param i uint represents a one based position in the first problem string
 * @param j uint represents a one based position in the second problem string
 * @return true iff the problems strings have matching symbols at position i & j
 ******************************************************************************/
bool LLCS2_SA_MQ::Phi(uint i, uint j) const {
  if (i == 0 || i > v.size()) return false;
  if (j == 0 || j > w.size()) return false;
  return v[i - 1] == w[j - 1];
}

/*******************************************************************************
 * @brief Calculates gap constraints via combining the maps left[b] and right[a]
 * @return SymbSymbTupleMap ( char x char -> std::pair<uint,uint> )
 ******************************************************************************/
LLCS2_SA_MQ::SymbSymbTupleMap LLCS2_SA_MQ::calcGapIntersection() const {
  std::unordered_map<Symbol, std::unordered_map<Symbol, std::pair<uint, uint>>> gapLR;
  for (const auto &[b, gapL] : left) {
    gapLR.emplace(b, std::unordered_map<Symbol, std::pair<uint, uint>>());
    for (const auto &[a, gapR] : right) {
      auto gap = std::make_pair(
          std::max<uint>(gapL.first, gapR.first),
          std::min<uint>(gapL.second, gapR.second)
      );
      std::unordered_map<Symbol, std::pair<uint, uint>> const x;
      gapLR.at(b).emplace(a, gap);
    }
  }
  return gapLR;
}

/*******************************************************************************
 * @brief Initializes MaxQueues for the algorithm
 * @note For D['b']['a'], 'b' is to the left of a gap and 'a' to its right
 *  (similar to gapLR)
 * @param gapLR SymbSymbTupleMap contains the gc infos
 * @return SymbSymbQueueMap (char x char -> MaxQueue2D)
 ******************************************************************************/
LLCS2_SA_MQ::SymbSymbQueueMap LLCS2_SA_MQ::getMaxQueue(SymbSymbTupleMap &gapLR) {
  using structures::MaxQueue2D;
  std::unordered_map<Symbol, std::unordered_map<Symbol, MaxQueue2D<uint >>> D;
  for (const auto &b : left | std::views::keys) {
    D.emplace(b, std::unordered_map<Symbol, MaxQueue2D<uint>>());
  }
  for (const auto &[b, right] : gapLR) {
    for (const auto &[a, gap] : right) {
      const uint yBuffer = std::max<uint>(1, gap.second - gap.first);
      D.at(b).emplace(a, MaxQueue2D<uint>(Mp[b], yBuffer, 0));
    }
  }
  return D;
}

/*******************************************************************************
 * @brief Does the doLineSetup for all MaxQueue2D[b][a] where b,a in Sigma
 * @param currentLine uint i while looping rows in the dp matrix M[i][j]
 * @param D SymbSymbQueueMap D[b][a] is a MaxQueue2D for b being to the
 * left of a gap and a to the right of a gap defined by gapLR[b][a]
 * @param gapLR SymbSymbTupleMap gapLR['b']['a'] is a gc-pair for a gap between
 * 'b' and 'a'
 ******************************************************************************/
void LLCS2_SA_MQ::doLineSetups(
    uint currentLine,
    LLCS2_SA_MQ::SymbSymbQueueMap &D,
    SymbSymbTupleMap &gapLR) {
  for (auto &[b, Db] : D) {
    for (auto &[a, Qba] : Db) {
      auto gap = gapLR[b][a];
      if (gap.first > gap.second)
        continue;
      const uint s0 = gap.second + 1;
      const uint e0 = gap.first + 1;
      Qba.doLineSetup(currentLine, s0, e0);
    }
  }
}

} // end of namespace