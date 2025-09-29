/*******************************************************************************
 * @file LLCS2_STD_H.cc
 * @author Steinkopp:Felix
 * @brief Hirschberg's LLCS Algorithm for two Strings
 * @details Time Complexity: \f$ \mathcal O \left(|s_1| \cdot |s_2|)\f$
 *          Space Complexity: \f$ \mathcal O \left(|s_2|)\f$
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_STD_H.h"

#include <algorithm>
#include <generator>
#include <memory>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"

namespace lcs_solver::algorithms::llcs {

constexpr std::string_view LLCS2_STD_H::description =
    R"DESC(Pseudocode: Hirschberg LLCS Algorithm for n strings
> LLCS(s1,s2,...,sk):
>   let li = length(si) and li<=lj for all i in {0,1}
>   let mid = l1 / 2
>   calc top_row    = dp[mid] via the folklore dp algorithm. Restrict memory to the two rows needed for the calculation.
>   calc bottom_row = dp[mid+1] similar to the folklore dp algorithm, but beginning in the bottom right corner. Restrict
                      memory to the two rows needed for the calculation.
>   return max_k( top_row[k] + bottom_row[k+1] )
)DESC";

/*******************************************************************************
 * Constructor for LLCS2_STD_H
 * @param vec The vector of shared points to the constant strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 ******************************************************************************/
LLCS2_STD_H::LLCS2_STD_H(const StrPtrVector& vec,
                         const ConstraintMap& map)
    : LLCS2_Algorithm(AlgoType::LLCS2_STD_H, vec, map, false, true) , llcs(0) {
}

/*******************************************************************************
 * @brief Checks if the algorithm is compatible with the given problem params
 * @return true, iff all four conditions hold true:
 * - the algorithm is given two strings, and
 * - the strings were sorted by length (increasingly)
 * - No string was empty
 * - Each constraint in the constraint map is valid
 ******************************************************************************/
bool LLCS2_STD_H::isValid() const {
  bool const res_is_sorted = isSorted();
  bool const res_is_filled = isFilled();
  bool const res_each_valid = isEachConstraintIndividuallyValid();
  return res_is_sorted && res_is_filled && res_each_valid && !trackKeyPairs;
}

/*******************************************************************************
 * @brief Getter for the algorithms pseudocode
 * @return std::string_view of the description string
 ******************************************************************************/
std::string_view LLCS2_STD_H::getDescription() const {
  return description;
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
std::unique_ptr<BaseSolution> LLCS2_STD_H::query() {
  reset(ResetLevel::Full);
  if (s.empty() || !isValid()) {
    return std::make_unique<solutions::EmptySolution>();
  }
  if (s.size() == 1) {
    return std::make_unique<solutions::UnsignedSolution>(s[0].size());
  }
  doPreprocessing();
  // uint lastElem = dp.empty() ? 0 : dp.back().back();
  return std::make_unique<solutions::UnsignedSolution>(llcs);
}

/*******************************************************************************
 * Resets the state of the LLCS2_STD_H object to its initial condition.
 * @note Side effect: Erases the elements from M and frees memory
 * @param lvl BaseAlgorithm::ResetLevel to be applied. Currently unused.
 ******************************************************************************/
void LLCS2_STD_H::reset(const ResetLevel lvl) {
  LLCS2_Algorithm::reset(lvl);
  llcs = 0;
}

std::string ToStringRows(const std::vector<uint>& top, const std::vector<uint>& bottom) {
  std::ostringstream oss;
  oss << "row1: ";
  for (const uint i : std::ranges::subrange(top.begin() + 0, top.end() - 0 )) {
    oss << i << " ";
  }
  oss << "\n";
  oss << "row2: ";
  for (const uint i : std::ranges::subrange(bottom.begin() + 0, bottom.end() - 0 )) {
    oss << i << " ";
  }
  oss << "\n";
  return oss.str();
}


/**
 * @brieft Two-dimensional Pair Generator in row-major order
 * @details Iterates a rectangle from top-left to bottom-right in row-major
 *          order. Example: `{a,b}, {a,b+1}, ..., {a,d}, ..., {c,b}, ... {c,d}`
 * @param top_left Pair `(a,b)` the upper-left corner of the area to iterate
 * @param bottom_right Pair `(c,d)` the bottom-right corner of the area
 * @return Sequence of row-major indices as a `std::generator`
 *
 * @pre `top_left.first > 0` && `top_left.second > 0`
 * @pre `top_left.first <= bottom_right.first`
 * @pre `top_left.second <= bottom_right.second`
 * @warning Uses `assert()` to enforce preconditions.
 */
std::generator<std::pair<uint, uint>> row_major_rect_inc(
    const std::pair<uint, uint> top_left, const std::pair<uint, uint> bottom_right) {
  assert(top_left.first > 0 && "Unsafe: top_left.first must be > 0");
  assert(top_left.second > 0 && "Unsafe: top_left.second must be > 0");
  assert(top_left.first <= bottom_right.first && "Invalid rectangle range");
  assert(top_left.second <= bottom_right.second && "Invalid rectangle range");
  for (uint i = top_left.first; i <= bottom_right.first; ++i)
    for (uint j = top_left.second; j <= bottom_right.second; ++j)
      co_yield {i, j};
}

/**
 * @brieft Two-dimensional Pair Generator in reverse row-major order
 * @details Iterates a rectangle from bottom-right to the top-left in row-major
 *          order. Example: `{c,d}, {c,d-1}, ..., {c,b}, ..., {a,d}, ... {a,b}`
 * @param top_left Pair `(a,b)` the upper-left corner of the area to iterate
 * @param bottom_right Pair `(c,d)` the bottom-right corner of the area
 * @return Sequence of row-major indices as a `std::generator`
 *
 * @pre `top_left.first > 0` && `top_left.second > 0`
 * @pre `top_left.first <= bottom_right.first`
 * @pre `top_left.second <= bottom_right.second`
 * @warning Uses `assert()` to enforce preconditions.
 */
std::generator<std::pair<uint, uint>> row_major_rect_dec(
  const std::pair<uint, uint> top_left, const std::pair<uint, uint> bottom_right) {
  assert(top_left.first > 0 && "Unsafe: top_left.first must be > 0");
  assert(top_left.second > 0 && "Unsafe: top_left.second must be > 0");
  assert(top_left.first <= bottom_right.first && "Invalid rectangle range");
  assert(top_left.second <= bottom_right.second && "Invalid rectangle range");
  for (uint i = bottom_right.first; i >= top_left.first; --i)
    for (uint j = bottom_right.second; j >= top_left.second; --j)
      co_yield {i, j};
}

/*******************************************************************************
 * @brief Executes the folklore llcs algorithm generalized to k strings
 ******************************************************************************/
void LLCS2_STD_H::doPreprocessing() {
  const uint l1 = s.size() < 2 ? 0 : s[0].size();
  const uint l2 = s.size() < 2 ? 0 : s[1].size();
  const uint mid = (l1 + 1) / 2;  // theoretically, there are l1+1 rows

  auto top_row = Row(l2 + 2, 0); // dp content is left and right zero-padded
  auto bottom_row = Row(l2 + 2, 0);

  // Step 1: Calc top_row from the top to mid
  auto current = Row(l2 + 2, 0); // next[0] and next[l2+1] stay zero
  std::ranges::for_each(
    row_major_rect_inc({1,1}, {mid, l2}),
    [&](const auto& index_pair) {
      const auto & [i, j] = index_pair;
      current[j] = isMatched(index_pair, true)
        ? top_row[j - 1] + 1
        : std::max(top_row[j], current[j-1]);
      if (j == l2) {
        std::swap(top_row, current); // set full row
      }
    }
  );

  // Step 2: Calc bottom_row from the bottom to mid+1
  std::ranges::for_each(
    row_major_rect_dec({mid+1,1}, {l1, l2}),
    [&](const auto& index_pair) {
      const auto & [i, j] = index_pair;
      current[j] = isMatched(index_pair, true)
        ? bottom_row[j + 1] + 1
        : std::max(bottom_row[j], current[j+1]);
      if (j == 1) {
        std::swap(bottom_row, current);
      }
    }
  );
  // std::cout << ToStringRows(top_row, bottom_row) << std::endl;

  // Step 3: Find LLCS: = max_k( top_row[k] + bottom_row[k+1] )
  auto sum = std::views::zip_transform(
    [](const uint a, const uint b) {return a+b;},
    std::views::counted(top_row.begin(), static_cast<int>(l2+1)), // Do not skip the empty string
    std::views::counted(bottom_row.begin()+1, static_cast<int>(l2+1))); // Connect cells diagonally

  llcs = std::ranges::max(sum);
  setState(State::Preprocessed);
}

/*******************************************************************************
 * @brief Creates a std::string to describe the state of the algorithm
 * @return std::string containing:
 * - the strings of the problem
 * - the dp matrix
 ******************************************************************************/
std::string LLCS2_STD_H::DebugString() const {
  std::stringstream oss;
  oss << BaseAlgorithm::toString(s) << "\n";
  oss << "llcs = " << llcs << "\n";
  return oss.str();
}

/*******************************************************************************
 * Getter for the dp-Matrix
 * @return structures::Matrix<Unsigned> dp
 * @note
 ******************************************************************************/
const LLCS2_STD_H::Matrix& LLCS2_STD_H::getMatrix() const {
  static Matrix empty_matrix;
  return empty_matrix;
}

/**
 * isExtensible
 * @param a Start Index (to be extended south-east to b)
 * @param b End Index (b1,b2)
 * @param llcs_of_a Uint equal to LLCS(s1[1:a1],s2[1:a2])
 * @return Whether an LCS from position a can be extended to position b
 */
bool LLCS2_STD_H::isExtensible(const Pair a,
                               const Pair b,
                               uint llcs_of_a) const {
  const bool match = isMatched(a) && isMatched(b);
  const bool gap_x = b.first > a.first && b.first - a.first >= 1;
  const bool gap_y = b.second > a.second && b.second - a.second >= 1;
  return match && gap_x && gap_y;
}

/**
 * getNextWindow
 * @param pair Index (x,y) for specifying the position s[1:x], s[1:y]
 * @param llcs The length of the longest common subsequence of s[1:x] and s[1:y]
 * @return {[1:x-1], [1:y-1]}
 */
LLCS2_Algorithm::Window LLCS2_STD_H::getPrevRange(const Pair& pair,
                                                  uint llcs) const {
  return {{1, pair.first - 1}, {1, pair.second - 1}};
}

}  // end of namespace lcs_solver::algorithms::llcs
