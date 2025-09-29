/******************************************************************************
 * @file LCS_STD_S.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Generalisation of the folklore LCs algorithm to k input strings
 * @details Time Complexity: \f$ \mathcal O \left( \prod_{k} len(string_k) \right)\f$
 *          Space Complexity: \f$ \mathcal O \left( \prod_{k} len(string_k) \right)\f$
 * @todo Compare Tests LCS_STD_S <=> LCS2_STD_S
 *****************************************************************************/

#include "algorithms/LCS/LCS_STD_S.h"
#include "util/Logger.hpp"
#include <cassert>
#include <sstream>

#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/Points.h"

#include <structures/MaxSegTreeND.h>

namespace lcs_solver::algorithms::lcs {

/*******************************************************************************
 * Constructor for LCS_STD_S
 * @param vec The vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 ******************************************************************************/
LCS_STD_S::LCS_STD_S(
    const StrPtrVector &vec,
    const ConstraintMap &map) : BaseAlgorithm(AlgoCategory::LLCS, AlgoType::LCS_Stack, vec, map),
                                algo(std::make_unique<LLCS_STD_FL>(vec, map)) {}

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true :
 * - the algorithm is given 2 string and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LCS_STD_S::isValid() const {
  bool const resIsSorted = isSorted();
  bool const resIsFilled = isFilled();
  bool const resEachValid = isEachConstraintIndividuallyValid();
  return resIsSorted && resIsFilled && resEachValid;
}

/*******************************************************************************
 * @brief Getter for the algorithms pseudocode
 * @return std::string_view of the description string
 ******************************************************************************/
//TODO: Update pseudocode
std::string_view LCS_STD_S::getDescription() const {
  static constexpr std::string_view msg = R"DESC(Pseudocode: LCS Folklore Algorithm for n strings
> LCS(s1,s2,...,sk):
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
>   ...
)DESC";
  return msg;
}

/*******************************************************************************
 * @brief Generates a string to describe the objects used in the algorithm
 * @return std::string with the strings and matrix in algo
 * @see LLCS_STD_FL::DebugString
 ******************************************************************************/
std::string LCS_STD_S::DebugString() const {
  return algo->DebugString();
}

/*******************************************************************************
 * @brief Executes LLCS_STD_FL and creates an appropriate LCS Solution
 * @return std::unique_ptr<BaseSolution> A unique pointer to a solution object:
 * - `EmptySolution` if the object is not valid or zero strings are provided
 * - `Points` equivalent to to s[1] if only one valid string is provided
 * - `LCS_STD_FL::StackSolution` to generate all possible LCS Sequences (
 ******************************************************************************/
std::unique_ptr<BaseSolution> LCS_STD_S::query() {
  reset(ResetLevel::Full);
  if (s.empty() || !isValid()) {
    return std::make_unique<solutions::EmptySolution>();
  }
  if (s.size() == 1) {
    std::vector<std::vector<uint>> mat(s[0].size() - 1);
    for (uint i = 0; i < mat.size(); i++) { mat[i].push_back(i); }
    return std::make_unique<solutions::Points>(getStrPtrVec(), mat);
  }
  algo->doPreprocessing();
  LLCS_STD_FL::Matrix mat = algo->getMatrix();
  return std::make_unique<StackSolution>(algo->getStrPtrVec(), mat);
}

/*******************************************************************************
 * @brief Executes the folklore llcs algorithm generalized to k strings
 * @see LLCS_STD_FL::doPreprocessing
 ******************************************************************************/
void LCS_STD_S::doPreprocessing() {
  algo->doPreprocessing();
  setState(State::Preprocessed);
}

/*******************************************************************************
 * Resets the state of the LLCS_STD_FL object to its initial condition.
 * @note Side effect: Erases the elements from M and frees memory
 * @param lvl BaseAlgorithm::ResetLevel to be applied. Currently unused.
 ******************************************************************************/
void LCS_STD_S::reset(BaseAlgorithm::ResetLevel lvl) {
  algo->reset(lvl);
  setState(algo->getState());
}

/*******************************************************************************
 * StackSolution Constructor
 * @param spv vector<shared_ptr<const String>> containing the problems' strings
 * @param matrix LLCS_STD_FL::Matrix filled DP-Matrix
 ******************************************************************************/
LCS_STD_S::StackSolution::StackSolution(
    const StrPtrVector &spv,
    const Matrix &matrix)
    : spv(spv), m(matrix) {}

/*******************************************************************************
 * StackSolution::begin
 * @return StackIterator for returning all possible LCS Subsequences
 ******************************************************************************/
solutions::BaseCollector::AnyIterator LCS_STD_S::StackSolution::begin() const {
  auto temp = AnyIterator(StackIterator(m, spv, false));
  return temp;
}

/*******************************************************************************
 * StackSolution::end
 * @return end of a StackIterator
 ******************************************************************************/
solutions::BaseCollector::AnyIterator LCS_STD_S::StackSolution::end() const {
  return AnyIterator(StackIterator(m, spv, true));
}

/*******************************************************************************
 * StackSolution::clone
 * @return new StackSolution based on stored Matrix and StringViewVec members
 ******************************************************************************/
BaseSolution *LCS_STD_S::StackSolution::clone() const {
  return new StackSolution(spv, m);
}

/*******************************************************************************
 * StackIterator Constructor
 * @param m
 * @param spv
 * @param end
 ******************************************************************************/
LCS_STD_S::StackSolution::StackIterator::StackIterator(
    const Matrix &m,
    const StrPtrVector &spv,
    const bool end)
    : spv(spv), dp(m), s(BaseAlgorithm::transformToSVV(spv)) {
  if (end || dp.empty() || dp[dp.lastIndex] == 0) {
    is_end = true;
    current = nullptr;
  } else {
    // util::Logger::Info() << describe("stack push", dp.linearIndex(dp.lastIndex));
    stack.emplace(dp.linearIndex(dp.lastIndex));
    sol = std::make_unique<Points>(spv);
    current = sol.get();
  }
}

/*******************************************************************************
 * @brief Changes the Points * current (or rather sol) to the next lcs
 ******************************************************************************/
void LCS_STD_S::StackSolution::StackIterator::advance() {
  while (!stack.empty()) {
    Matrix::uint pos = stack.top();
    stack.pop();
    // util::Logger::Info() << describe("stack pop", pos);
    if (dp[pos] == 1 && matchingSymbol(dp, pos, s)) {
      // util::Logger::Info() << "return" << DebugString();
      return;
    }
    // Add Positions onto stack by removing one symbol from a string (this is done for all k strings)
    for (size_t k = 0; k < dp.dim.size(); ++k) {
      if (dp.vectorizeIndex(pos, k) > 0) {// non trivial position
        const Matrix::uint smaller = dp.getStepDownIndex(pos, k);
        if (dp[pos] == dp[smaller]) {
          // util::Logger::Info() << describe("stack push", smaller);
          stack.push(smaller);
        }
      }
    }

    // Look for a path by checking the existence of a matching symbol path
    const uint diag = dp.stepDownDiag(pos);
    if (dp[diag] == dp[pos] - 1
        && matchingSymbol(dp, pos, s)
        && matchingSymbol(dp, diag, s)) {
      if (!sol->empty()) { // Remove
        uint lastSolPosition = dp.linearIndex(sol->back());
        while (!sol->empty() && dp[lastSolPosition] >= dp[pos]) {
          // util::Logger::Info() << describe("stack pop2", stack.top());
          stack.pop();
          lastSolPosition = dp.linearIndex(sol->back());
        }
      }
      if (sol->empty() || pos != diag) {
        util::Logger::Info() << describe("stack push", diag);
        stack.push(diag);
        auto temp = dp.vectorizeIndex(pos);
        sol->emplace_back();
      }
    } // extension possible
  } // while(!stack.empty()
}

/*******************************************************************************
 * clone
 * @return Clone that respects spv, dp, s, sol, stack, current and is_end
 ******************************************************************************/
std::unique_ptr<solutions::BaseIterator> LCS_STD_S::StackSolution::StackIterator::clone() const {
  auto clonePtr = std::make_unique<StackIterator>(dp, spv, is_end); // sets spv, dp, s, is_end
  clonePtr->sol = std::make_unique<Points>(*sol);
  clonePtr->current = clonePtr->sol.get(); // basically a const cast in the respective object
  clonePtr->stack = stack;
  return clonePtr;
}

/*******************************************************************************
 * DebugString
 * @return string about the current value of the iterator and the end flag
 ******************************************************************************/
std::string LCS_STD_S::StackSolution::StackIterator::DebugString() const {
  auto oss = std::ostringstream();
  oss << "current:" << sol->DebugString() << "\n"
      << "is_end:" << is_end << "\n";
  return BaseIterator::DebugString();
}

/*******************************************************************************
 * @brief Helper - Checks for matching symbols in the strings
 * @param mat matrix
 * @param pos Matrix::uint linearized (row-major) Index in the dp matrix
 * @param svv std::vector<StringView> containing the problems' strings
 * @return false if positions do not match. True in the cases:
 * - svv[pos[0]] == ... == svv[pos[pos.size()-1]] && !oneBased
 ******************************************************************************/
bool LCS_STD_S::StackSolution::StackIterator::matchingSymbol(
    const Matrix &mat,
    const Matrix::uint &pos,
    const StrViewVector &svv) {
  for (size_t k = 1; k < mat.dim.size(); ++k) {// k-1 comparisons
    Matrix::uint i = mat.vectorizeIndex(pos, k);
    if (i == 0 || i > svv[k].size()) {
      return false;// svv[k][i] is not in [1:ssv[k].size()]
    }
    if (svv[0][i - 1] != svv[k][i - 1])
      return false;
  }
  return true;
}

/*******************************************************************************
 * @brief Describes an action when doing something at a position in a matrix
 * @param msg string with a meaning full message
 * @param pos uint equal to some position in the dp matrix of a StackIterator
 * @return std::string containing msg, pos, and the vector form of pos
 ******************************************************************************/
std::string LCS_STD_S::StackSolution::StackIterator::describe(const std::string &msg, const Matrix::uint pos) const {
  std::ostringstream oss;
  std::vector<Matrix::uint> vec = dp.vectorizeIndex(pos); // readable index
  oss << msg << "(pos: " << pos << " ~ (";
  for (const auto & index : vec) {
    oss << index << " ";
  }
  oss << ")" << std::endl;
  return oss.str();
}

}// namespace lcs_solver::algorithms::lcs