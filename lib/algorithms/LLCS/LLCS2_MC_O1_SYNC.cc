/*******************************************************************************
 * @file LLCS2_MC_O1_SYNC.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Implementation of an algorithm for LLCS2_MC_O1_SYNC
 * @details Time Complexity: O(n*m)
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_MC_O1_SYNC.h"

#include <algorithm>
#include <limits>
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
#include "util/Cantor.h"
#include "util/RadixSort.h"

namespace lcs_solver::algorithms::llcs {
/*******************************************************************************
 * Constructor for LLCS2_MC_O1_SYNC
 * @param vec The vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 ******************************************************************************/
LLCS2_MC_O1_SYNC::LLCS2_MC_O1_SYNC(const StrPtrVector &vec, const ConstraintMap &map)
    : LLCS2_MC_Algorithm(
          AlgoType::LLCS2_MC_O1_SYNC, vec, map,
          LLCS2_MC_Algorithm::getGaps(map, ConstraintType::MC_O1C_SYNC)),
      v(s.size() == 2 ? s[0] : util::StringView()),
      w(s.size() == 2 ? s[1] : util::StringView()),
      m(s.size() == 2 ? s[0].size() : 0),
      n(s.size() == 2 ? s[1].size() : 0) {}

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS2_MC_O1_SYNC::isValid() const {
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
std::string_view LLCS2_MC_O1_SYNC::getDescription() const {
  return {R"(
LLCS Algorithm h distinct gc-constraint-pairs. The gaps are synchronized
/** Pseudocode (Adamson et. al, LCS wit Gap constraints page 14-16) ************
 *  LLCS2_MC_O1_SYNC(s1,s2,gc):
 *      // Initialize the relevant data structures
 *      let li = length(si) and li<=lj for all i<j
 *      let Phi(s1[i], s2[j]) = ith symbol of s1 is equal to jth symbol of s2
 *      RadixSort the gc-tuples C and reducing them to list Cs of h distinct ones (label[i]=j iff C[i]=Cs[j])
 *      let M_r = Zeros[lk + 1][lk + 1] for all r in [h]
 *      let D_rp = MaxQueue2D(const MatrixPtr = &M, xRange = [0..l1], yRange[0..l2], bufferAxisSize = u-l) for rp in [h]
 *
 *      // Fill M_p using dynamic programming and a window sliding approach. M_p is computed in O(l1*l2)
 *      for i in [0..l1]
 *          update D_rp for rp in [h]
 *          for j in [0..l2]
 *              update D_rp for rp in [h]
 *              for rp in [h]
 *                  let l,u = Cs[rp]
 *                  let s0 = u+1, e0=l+1, s1=u+1, e1=l+1
 *                  use D_rp to retrieve m_rp, the maximum of the M[I_rp][J_rp] where I_rp = [i-s0:i-e0] and J_rp = [j-s1:j-e1]
 *                  m_rp is set to be 0 when I or J are empty
 *              compute set max_r = max{m_rp| rp in [h] and label[m_rp+1]=r} for all r in [h]
 *              for r in [h]
 *                  if Phi(s1[i], s2[j]) == 1 then
 *                      set M[i,j] = m + 1
 *                  else
 *                      set M[i][j] = max {M[i-1][j] , M[i][j-1])}
 *      return max {M[r][x][y] | r in [1..h],l in [1..l1] and y in [1..l2]
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
std::unique_ptr<BaseSolution> LLCS2_MC_O1_SYNC::query() {
  reset(ResetLevel::Full);
  if (!isValid()) return std::make_unique<solutions::EmptySolution>();
  doPreprocessing();
  return std::make_unique<solutions::UnsignedSolution>(findMaxIn3DMatrix(M));
  // return std::make_unique<solutions::UnsignedSolution>(R[m][n]);
}

/*******************************************************************************
 * @brief Getter for the dp-Matrix
 * @return std::vector<std::vector<uint> M
 *******************************************************************************/
const LLCS2_Algorithm::Matrix &LLCS2_MC_O1_SYNC::getMatrix() const {
  return R;
}

/*******************************************************************************
 * Function used to check for matching symbols in the strings
 * @param i position in the first string
 * @param j position in the second string
 * @return true if the ith position of the 1st string is equal to the jth of the
 *  2nd string
 ******************************************************************************/
bool LLCS2_MC_O1_SYNC::Phi(const uint i, const uint j) const {
  if (i == 0 || i > v.size()) return false;
  if (j == 0 || j > w.size()) return false;
  return v[i - 1] == w[j - 1];
}

/*******************************************************************************
 * @brief Executes the main dp loop of the algorithm
 * @details Radix Sorts the gaps and groups the gaps into `Cp`. Then it sets up
 *  the matrices `R` and `M_r`and uses a `MaxQueue2D` to fill `R[i][j]` with the
 *  llcs of s_1[1:i], s_2[1:j]
 ******************************************************************************/
void LLCS2_MC_O1_SYNC::doPreprocessing() {
  using structures::MaxQueue2D;
  using util::Cantor;
  using util::RadixSort;

  // Radix-sort the gap constraint vector
  std::vector<std::pair<size_t, std::pair<size_t,size_t>>> Cp;// Vector for key value pairs
  for (const auto &gap : C) {
    Cp.emplace_back(Cantor::calc(gap), gap);
  }
  auto [Cs, label] = RadixSort<size_t, std::pair<size_t,size_t>>::sort(Cp, true);

  // Setup D[r] for r in [0:h-1]
  const uint h = Cs.size();
  R = Matrix2D(m + 1, Row(n + 1, 0));
  M.reserve(h);
  for (uint r = 0; r < h; ++r) {
    M.emplace_back(m + 1, Row(n + 1, 0));
  }

  // Setup D[r] for M[r]
  std::vector<MaxQueue2D<uint>> D;
  for (uint r = 0; r < h; r++) {
    const uint yBuffer = std::max<uint>(1, Cs[r].second.second - Cs[r].second.first);
    D.emplace_back(M[r], yBuffer, 0);// D[r] operatees on M[r]
  }

  for (uint i = 0; i <= m; ++i) {
    for (uint r = 0; r < h; r++) {
      const uint s0 = Cs[r].second.second + 1;
      const uint e0 = Cs[r].second.first + 1;
      D[r].doLineSetup(i, s0, e0);
    }
    for (uint j = 0; j <= n; ++j) {
      auto mp = std::vector<uint>(h, 0);
      for (uint r = 0; r < h; r++) {
        const uint s0 = Cs[r].second.second + 1;
        const uint e0 = Cs[r].second.first + 1;
        const uint s1 = Cs[r].second.second + 1;
        const uint e1 = Cs[r].second.first + 1;
        D[r].doComputingSetup(i, j, e0, e1);
        mp[r] = D[r].query(i, j, {s0, e0, s1, e1});
      }
      //uint llcs = *std::max_element(p.begin(), p.end());
      // Compute set max_r = max_{rp in [h], label[m_rp+1]=r} for all r in [h]
      for (uint r = 0; r < h; ++r) {
        if (Phi(i, j)) {
            uint maxr = 0;
            for (uint index_in_cs = 0; index_in_cs < h; ++index_in_cs) {
              const uint candidate_length = mp[index_in_cs];
              if (candidate_length > 0) {
                const uint candidate_index_in_cs = label[candidate_length-1];
                if (mp[candidate_index_in_cs] == candidate_length) {
                  maxr = std::max<uint>(maxr, candidate_length);
                }
              }
            }
            M[r][i][j] = maxr + 1;

            // Tracking of points M[r][i][j]
            const uint& len = M[r][i][j];
            const Pair idx = {i, j};
            const bool do_update = keyPairs[len].empty() || (keyPairs[len].back() != idx);
            if (trackKeyPairs && do_update) track(idx, len);
            R[i][j] = std::max<uint>(R[i][j], M[r][i][j]);
        } else if (i > 0 && j > 0) {
          R[i][j] = std::max<uint>(R[i - 1][j], R[i][j - 1]);
        }
      }
    }
  }
  setState(State::Preprocessed);
}

/*******************************************************************************
 * @brief Finds the largest value in a 3D Matrix
 * @param matrix std::vector<std::vector<std::vector<uint>>>
 * @return uint that is the largest value found in matrix
 ******************************************************************************/
LLCS2_MC_O1_SYNC::uint LLCS2_MC_O1_SYNC::findMaxIn3DMatrix(const Matrix3D &matrix) {
  if (matrix.empty())
    return std::numeric_limits<uint>::min();
  uint maxVal = matrix[0][0][0];
  for (const auto &mat2D : matrix) {
    for (const auto &vec1D : mat2D) {
      for (const auto &elem : vec1D) {
        maxVal = std::max<uint>(maxVal, elem);
      }
    }
  }
  return maxVal;
}

/*******************************************************************************
 * @brief Resets the state of the LLCS2_MC_O1_SYNC object to its initial state
 * @details This function clears the matrix `M` and sets the object's state to
 *  `Constructed`.
 * @param lvl The level of reset to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_MC_O1_SYNC::reset(BaseAlgorithm::ResetLevel lvl) {
  M.clear();
  setState(State::Constructed);
}

/*******************************************************************************
 * @brief Creates a std::string to describe the state of the algorithm
 * @return std::string containing:
 * - the dp matrix M_r for all r in [h]
 * - the dp matrix R (with R[i][j] = llcs of s_1[1:i] and s_2[1:j])
 * - the maximum found in M_r over all r in [h]
 ******************************************************************************/
std::string LLCS2_MC_O1_SYNC::DebugString() const {
  std::ostringstream oss;
  oss << BaseAlgorithm::toString(s) << "\n";
  if (M.empty()) {
    oss << "dp is empty\n";
  } else {
    uint r = 0;
    for (const auto &mat2D : M) {
      std::string matrix_name = "M[r == " + std::to_string(r++) + "]";
      oss << toString(mat2D, s, false, matrix_name, true) << "\n";
    }
  }
  // oss << toString(R, s, false, "R", true) << "\n";
  oss << "Maximum is: " << findMaxIn3DMatrix(M) << "\n";
  oss << toString(keyPairs, trackKeyPairs) << "\n";
  return oss.str();
}

}// namespace lcs_solver::algorithms::llcs