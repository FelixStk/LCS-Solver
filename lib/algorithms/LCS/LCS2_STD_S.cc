/******************************************************************************
 * @file LCS2_STD_S.cc
 * @author Steinkopp:Felix
 * @version 0.1
 * @brief Implementation of an algorithm for LCS2_Classic with RangeTrees
 *****************************************************************************/

#include "algorithms/LCS/LCS2_STD_S.h"

#include <sstream>
#include <utility>

#include "algorithms/AlgoFactory.h"
#include "algorithms/LLCS/LLCS2_Algorithm.h"
#include "algorithms/LLCS/LLCS2_STD_FL.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/Points.h"

namespace lcs_solver::algorithms::lcs {

/**
 * Constructor for LCS2_STD_S
 * @param spv The vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 */
LCS2_STD_S::LCS2_STD_S(const StrPtrVector &spv, const ConstraintMap &map)
    : BaseAlgorithm(AlgoCategory::LLCS, AlgoType::LCS2_STD_S, spv, map),
      algo_(std::make_unique<llcs::LLCS2_STD_FL>(spv, map)) {
  algo_->setTracking(false);// The dp[i][j] == LLCS(s1[1:i],s2[1:j])
}

/**
 * isValid
 * @return true, iff algo is valid, inherits from LLCS_Algorithm
 */
bool LCS2_STD_S::isValid() const {
  if (algo_ == nullptr) {
    return false;// *algo did not inherit from llcs::LLCS2_Algorithm
  }
  bool const svv_is_of_length_two = algo_->getStrPtrVec().size() == 2;
  bool const algo_is_valid = algo_->isValid();
  return algo_is_valid && svv_is_of_length_two;
}

/**
 * Getter for the algorithms pseudocode
 * @return std::string_view of the description string
 */
std::string_view LCS2_STD_S::getDescription() const {
  static constexpr std::string_view msg = R"DESC(Pseudocode: Reconstruction of LCS with Rangetree
> Initialization:
>   Let s1 and s2 be the two input sequences
>   Let l1 = length(s1), l2 = length(s2)
>   Let dp be a dp matrix with : dp[i][j] == LLCS(s0[1:i], s1[1:j])
>   Let s be a stack of states, where a state has the form ((i,j), len)
>   Let p be a vector for points (pairs (i,j)) representing the embeddings
>
> Function next() -> vector<pair<uint,uint>>:
>   If firstCall then s.push( (l1,l2,0) )
>   While s is not empty:
>     (pair, len) = stack.pop()
>     If p.size() > len then p.resize(len)
>
>     If i==0 or j == 0 then:
>       return p
>     If (max = std::max(dp[i][j-1],dp[i-1][j]); max==dp[i][j]) then:
>       If dp[i-1][j] >= dp[i][j-1] then s.push((i-1, j, len))
>       If dp[i][j-1] >= dp[i-1][j] then s.push((i, j-1, len))
>     If (len == 0 && isMatched(pair)) || (len > 0 && isExtensible(pair, p.back(), len)) then:
>       p.push_back((i,j))
>       s.push ((i-1, j-1, len+1))
>   return nullptr and set is_end flag to true
)DESC";
  return msg;
}

/**
 * DebugString
 * @return  std::string containing the type and DebugSting of algo
 */
std::string LCS2_STD_S::DebugString() const {
  std::ostringstream oss;
  oss << "type: " << static_cast<int>(type) << "\n";
  oss << algo_->DebugString();
  return oss.str();
}

/**
 * @brief Executes algo and creates an appropriate LCS RangeTree Solution
 * @return std::unique_ptr<BaseSolution> A unique pointer to a solution object:
 * - `EmptySolution` if the object is not valid or zero strings are provided
 * - `Points` equivalent to to s[1] if only one valid string is provided
 * - `LCS_STD_FL::StackSolution` to generate all possible LCS Sequences (
 */
std::unique_ptr<BaseSolution> LCS2_STD_S::query() {
  reset(ResetLevel::Full);
  if (!isValid()) {
    return std::make_unique<solutions::EmptySolution>();
  }
  if (getState() != State::Preprocessed)
    algo_->doPreprocessing();
  return std::make_unique<StackSolution>(algo_);
}

/**
 * @brief Executes the llcs algorithm
 * @see LLCS_STD_FL::doPreprocessing
 */
void LCS2_STD_S::doPreprocessing() {
  algo_->doPreprocessing();
  setState(algo_->getState());
}

/**
 * Resets the algo that is used to get the lcs points in the StackSolution
 * @param lvl BaseAlgorithm::ResetLevel to be applied
 */
void LCS2_STD_S::reset(const ResetLevel lvl) {
  algo_->reset(lvl);
  setState(algo_->getState());
}

/**
 * StackSolution - Constructor
 * @param ptr LLCS2_Algorithm for getting the necessary reconstruction data
 */
LCS2_STD_S::StackSolution::StackSolution(const std::shared_ptr<llcs::LLCS2_Algorithm> &ptr)
    : algo_(ptr) {}

/**
 * StackSolution::begin
 * @return Iterator for returning all possible LCS Subsequences
 */
LCS2_STD_S::StackSolution::AnyIterator LCS2_STD_S::StackSolution::begin() const {
  return AnyIterator(Iterator(algo_, false));
}

/**
 * StackSolution::end
 * @return Iterator for returning all possible LCS Subsequences
 */
LCS2_STD_S::StackSolution::AnyIterator LCS2_STD_S::StackSolution::end() const {
  return AnyIterator(Iterator(algo_, true));
}

/**
 * StackSolution::clone
 * @return cloned StackSolution (can be created by calling the constructor)
 */
BaseSolution *LCS2_STD_S::StackSolution::clone() const {
  return new StackSolution(algo_);
}

/**
 * Constructor for Iterator
 * @details The default constructor is used to generate the end for an iterator
 */
LCS2_STD_S::StackSolution::Iterator::Iterator()
    : algo_(nullptr),
      current_(std::make_unique<solutions::Points>()),
      llcs_(0) {
  current = current_.get();
}

/**
 * Constructor for Iterator
 * @param algo weak_ptr<LLCS2_Algorithm> to access `LLCS2_Algorithm::isMatched`
 *  and `LLCS2_Algorithm::isExtensible`. Also used to generate the range trees
 * @param end Whether the iterator has finished processing
 */
LCS2_STD_S::StackSolution::Iterator::Iterator(
    const std::weak_ptr<llcs::LLCS2_Algorithm> &algo,
    const bool end)
    : algo_(std::move(algo.lock())),
      current_(std::make_unique<solutions::Points>(algo_->getStrPtrVec(), true)),
      dp_(algo_->getMatrix()),
      llcs_(dp_.empty() ? 0 : (dp_.back().empty() ? 0 : dp_.back().back())),
      processed_(llcs_ + 1) {
  const auto &s = algo_->getStringViewVec();
  const Pair last_point_in_dp_mat = std::make_pair(s[0].size(), s[1].size());
  const State state = std::make_pair(last_point_in_dp_mat, 0);
  stack_.push(state);
  is_end = end;
  current = current_.get();
  advance();
}

/**
 * advance() executes the pseudocode writen down in LCS2_STD_S::getDescription
 */
void LCS2_STD_S::StackSolution::Iterator::advance() {
  if (is_end)
    return;

  solutions::Points &points = *current_;
  while (!stack_.empty()) {
    const auto [pair, len] = stack_.top();
    const auto &[i, j] = pair;
    stack_.pop();

    // Base Case: Report solutions
    if (i == 0 || j == 0) {
      if (points.empty())
        continue;
      return;
    }
    // Correct Content when a new solution is started
    const bool is_pair_known = processed_[dp_[i][j]].contains(pair);
    if (algo_->isMatched(pair, true) && !is_pair_known) {
      if (points.size() > len) {
        for (uint k = 1; k < llcs_ - len; ++k) {
          processed_[k].clear();
        }
        points.resize(len);
      }
    }
    // Push onto the stack without extending points
    if (dp_[i][j] == std::max(dp_[i][j - 1], dp_[i - 1][j])) {
      if (dp_[i - 1][j] >= dp_[i][j - 1]) {
        stack_.emplace(std::make_pair(i - 1, j), len);
      }
      if (dp_[i][j - 1] >= dp_[i - 1][j]) {
        stack_.emplace(std::make_pair(i, j - 1), len);
      }
    }
    // Push onto the stack if extending points is possible
    const bool b1 = len == 0 && algo_->isMatched(pair, true);
    const bool b2 = len > 0 && algo_->isExtensible(pair, points.BackPair(), dp_[i][j]);
    if (b1 || b2) {
      if (!is_pair_known) {
        points.push_back(pair);
        processed_[llcs_ - len].insert(pair);
        stack_.emplace(std::make_pair(i - 1, j - 1), len + 1);
      }
    }
  }
  is_end = true;
  current_->clear();
}

/**
 * clone creates a deep copy of the iterator (including its state)
 * @return std::unique_ptr<solutions::BaseIterator>
 */
std::unique_ptr<solutions::BaseIterator> LCS2_STD_S::StackSolution::Iterator::clone() const {
  auto ptr = std::make_unique<Iterator>();
  // Clone fields from LCS2_STD_S
  const auto p = dynamic_cast<solutions::Points *>(current_->clone());
  if (p == nullptr) return nullptr;
  ptr->algo_ = this->algo_;
  ptr->current_ = std::unique_ptr<solutions::Points>(p);
  ptr->stack_ = stack_;
  // Clone set fields from BaseIterator
  ptr->current = ptr->current_.get();
  ptr->is_end = is_end;
  return ptr;
}

}// namespace lcs_solver::algorithms::lcs