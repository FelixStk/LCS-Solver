/*******************************************************************************
 * @file LLCS_STD_FL.cc
 * @author Steinkopp:Felix
 * @brief Generalisation of the folklore LCs algorithm to k input strings
 * @details Time Complexity: \f$ \mathcal O \left( \prod_{k} len(string_k) \right)\f$
 *          Space Complexity: \f$ \mathcal O \left( \prod_{k} len(string_k) \right)\f$
 ******************************************************************************/

#include "algorithms/LLCS/LLCS_STD_FL.h"

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"

namespace {
using ::lcs_solver::algorithms::solutions::EmptySolution;
using ::lcs_solver::algorithms::solutions::UnsignedSolution;
}// namespace

namespace lcs_solver::algorithms::llcs {

/*******************************************************************************
 * Constructor for LLCS_STD_FL
 * @param vec The vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 ******************************************************************************/
LLCS_STD_FL::LLCS_STD_FL(
    const BaseAlgorithm::StrPtrVector &vec,
    const ConstraintMap &map
) : LLCS_Algorithm(AlgoType::LLCS_STD_FL, vec, map) {}

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS_STD_FL::isValid() const {
  bool const resIsSorted = isSorted();
  bool const resIsFilled = isFilled();
  bool const resEachValid = isEachConstraintIndividuallyValid();
  return resIsSorted && resIsFilled && resEachValid;
}

/*******************************************************************************
 * @brief Getter for the algorithms pseudocode
 * @return std::string_view of the description string
 ******************************************************************************/
std::string_view LLCS_STD_FL::getDescription() const {
  static constexpr std::string_view msg = R"DESC(Pseudocode: LCS Folklore Algorithm for n strings
> LLCS(s1,s2,...,sk):
>   let li = length(si) and li<=lj for all i<j
>   let dp[l1 + 1][l2 + 1]...[lk + 1] filled with zeros
>   for i1 = 1 to l1 do
>     for i2 = 1 to l2 do
>       ...
>         for ik = 1 to lk do
>           if s1[i1-1] == s2[i2-1] == ... == sk[ik-1]
>             dp[i1][i2]...[ik] = dp[i1-1][i2-1]...[ik-1] + 1
>           else
>             dp[i1][i2]...[ik] = max( dp[i1-1][i2]...[ik], dp[i1][i2-1]...[ik], ... , dp[i1][i2]...[ik-1] )
>   return dp[l1][l2]...[lk]
)DESC";
  return msg;
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
std::unique_ptr<BaseSolution> LLCS_STD_FL::query() {
  reset(ResetLevel::Full);
  if (s.empty() || !isValid()) {
    return std::make_unique<EmptySolution>();
  }
  if (s.size() == 1) {
    return std::make_unique<UnsignedSolution>(s[0].size());
  }
  doPreprocessing();
  auto lastElem = dp[dp.lastIndex];
  return std::make_unique<UnsignedSolution>(lastElem);
}

/*******************************************************************************
 * Resets the state of the LLCS_STD_FL object to its initial condition.
 * @note Side effect: Erases the elements from M and frees memory
 * @param lvl BaseAlgorithm::ResetLevel to be applied. Currently unused.
 ******************************************************************************/
void LLCS_STD_FL::reset(BaseAlgorithm::ResetLevel lvl) {
  dp.clear();
}

/*******************************************************************************
 * Getter for the dp-Matrix
 * @return structures::Matrix<Unsigned> dp
 *******************************************************************************/
const LLCS_STD_FL::Matrix &LLCS_STD_FL::getMatrix() {
  if (getState() != State::Preprocessed)
    doPreprocessing();
  return dp;
}

/*******************************************************************************
 * @brief Executes the folklore llcs algorithm generalized to k strings
 * @note This function uses the k dim Matrix<uint> dp (that is based on a std::vector<uint>)
 ******************************************************************************/
void LLCS_STD_FL::doPreprocessing() {
  // Create a |s_d| x |s_{d-1}| x ... x |s0| matrix (with d = #strings - 1 )
  std::vector<size_t> dim(s.size());
  std::ranges::generate(dim, [n = 0, this]() mutable { return s[n++].size() + 1; });
  dp = Matrix(dim);// filled with zeros

  // Folklore Algorithm
  for (uint pos = 0; pos < dp.size(); pos++) {
    auto idx = dp.vectorizeIndex(pos);
    if (matchingSymbol(idx, s))
      dp[pos] = dp.stepDownDiag(pos) + 1;
    else
      for (uint d = 0; d < idx.size(); ++d)
        dp[pos] = std::max(dp[pos], dp.stepDown(pos, d));
  }
  setState(State::Preprocessed);
}

/******************************************************************************
 * @brief Checks for matching symbols in the strings
 * @param idx std::vector<uint> such svv[k][idx[k]] represents the idx[k]-th
 * symbol in the k-th string of the problem
 * @param svv reference to the vector of share_ptr for the problem strings
 * @return svv[0][idx[1]-1] == svv[1][idx[1]-1] == ... == svv[k-1][idx[k-1]-1]
 *****************************************************************************/
bool LLCS_STD_FL::matchingSymbol(
    const std::vector<size_t> &idx,
    const StringViewVector &svv) {
  for (uint k = 1; k < idx.size(); ++k) {
    if (idx[k] < 1 || idx[k] > svv[k].size())
      return false;
    if (svv[0][idx[0] - 1] != svv[k][idx[k] - 1])
      return false;
  }
  return true;
}

/*******************************************************************************
 * @brief Creates a std::string to describe the state of the algorithm
 * @return std::string containing:
 * - the strings of the problem
 * - the dp matrix
 ******************************************************************************/
std::string LLCS_STD_FL::DebugString() const {
  std::ostringstream ss;
  for (uint i = 0; i < s.size(); ++i)
    ss << "s[" << i << "]: " << util::to_string(s[i]) << "\n";
  ss << dp.DebugString();
  return ss.str();
}

} // end of namespace lcs_solver::algorithms::llcs
