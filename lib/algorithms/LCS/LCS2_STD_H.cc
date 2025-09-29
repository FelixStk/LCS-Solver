/******************************************************************************
 * @file LCS2_STD_H.cc
 * @author Steinkopp:Felix
 * @version 0.2
 * @brief Implementation of LCS Hirschberg algorithm with return Type LCS iter
 *****************************************************************************/

#include "algorithms/LCS/LCS2_STD_H.h"

#include <cassert>
#include <format>
#include <generator>
#include <list>
#include <ranges>
#include <sstream>
#include <optional>
#include <utility>

#include "algorithms/AlgoFactory.h"
#include "algorithms/LLCS/LLCS2_Algorithm.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/Points.h"
#include "util/Logger.hpp"

#define LCS_SOLVER_ENABLE_HIRSCH_LOGGING 0

namespace lcs_solver::algorithms::lcs {

using util::StringView; // = std::string_view
using Row = std::vector<uint>;
using Pair = LCS2_STD_H::Pair; // = std::pair<size_t, size_t>
using PairVec = LCS2_STD_H::PairVec; // std::vector<Pair>
using InitList = std::initializer_list<Pair>;
using RowPair = std::pair<Row, Row>;
using StrViewPair = std::pair<StringView, StringView>;
using RowIter = Row::iterator;
using OptIdx = std::optional<uint>;
using OptIdxPair = std::pair<OptIdx, OptIdx>;
using RowIterPair = std::pair<RowIter, RowIter>;

//==============================================================================
// Declare logging functions and macros
//==============================================================================
void log_pairs(std::string_view msg, const PairVec& vec, int depth);
void log_sv(std::string_view msg, StringView s, StringView t, int depth = 0);
void log_filter_msg(std::string_view msg, int lvl, int depth = 0);

template <typename T>
void log_vector(std::string_view msg, const std::vector<T>& vec, int depth = 0);

#if LCS_SOLVER_ENABLE_HIRSCH_LOGGING
#define LOG_CALL(s, t, d) log_sv(("call"), (s), (t), (d))
#define LOG_D0_MSG(msg, d) log_filter_msg((msg), (0), (d))
#define LOG_PAIRS(msg, vec, d) log_pairs((msg), (vec), (d))
#define LOG_ROW(msg, vec, d) log_vector((msg), (vec), (d))
#define LOG_SV(msg, s, t, d) log_sv((msg), (s), (t), (d))
#else
#define LOG_CALL(s, t, d)
#define LOG_D0_MSG(msg, d)
#define LOG_PAIRS(msg, vec, d)
#define LOG_ROW(msg, vec, d)
#define LOG_SV(msg, s, t, d)
#endif

//==============================================================================
// Declare Auxiliary Functions
//==============================================================================
bool IsMatch(const Pair& pair, StringView s, StringView t);

std::generator<Pair> row_major_rect_inc(Pair top_left, Pair bottom_right);
std::generator<std::pair<uint, uint>> row_major_rect_dec(Pair top_left, Pair bottom_right);

LCS2_STD_H::Generator gen_empty();
LCS2_STD_H::Generator base_cases(StringView s, StringView t, bool swapped, int depth);
LCS2_STD_H::Generator middle_generator(StringView s, StringView t, const RowPair& top_pair, const RowPair& bottom_pair, int depth);

Row dp_hirsch_top(StringView s, StringView t);
Row dp_hirsch_bottom(StringView s, StringView t);
RowPair dp_extended_top(StringView s, StringView t);
RowPair dp_extended_bottom(StringView s, StringView t);

StrViewPair GetUpperLeft(StringView s, StringView t, const PairVec& middle);
StrViewPair GetLowerRight(StringView s, StringView t, const PairVec& middle);

PairVec Merge(StringView s, StringView t, const PairVec& l, const PairVec& m, const PairVec&r, bool swap_l, bool swap_m, bool swap_r);

std::generator<Pair> GenTopPairs(StringView s,StringView t, uint row, uint init_col, uint n);
std::generator<Pair> GenBtmPairs(StringView s,StringView t, uint row, uint init_col, uint n);

std::generator<OptIdxPair> GenerateKs(const Row& dp_top, const Row& dp_btm, const Row& pos_top, const Row& pos_btm, uint max_value);

//==============================================================================
// Main Algorithm
//==============================================================================
/**
 * @brief Implementation of LCS Hirschberg algorithm using c++20 coroutines
 * @param s First String
 * @param t Second String
 * @param depth Current depth in the recursion. Only used for debug output.
 * @return Generator with of type std::vector<std::pair<uint, uint>>
 */
LCS2_STD_H::Generator hirschberg(StringView s, StringView t, const int depth = 0) {
  LOG_CALL(s, t, depth);
  // Maintain the invariant s.length() >= t.length()
  const bool sv_swap = s.length() < t.length();
  if (sv_swap) {
    std::swap(s, t);
  }


  // Step 1: Handle Base Cases
  if (t.length() < 2) {
    for (auto&& pair_vector : base_cases(s, t, sv_swap, depth)) {
      co_yield pair_vector;
    }
    co_return;
  }

  // Step 2: Do Classy Dynamic Programming
  const RowPair top_row_pair = dp_extended_top(s, t);
  const RowPair bottom_row_pair = dp_extended_bottom(s, t);
  LOG_ROW("top_dp", top_row_pair.first, depth);
  LOG_ROW("top_n ", top_row_pair.second, depth);
  LOG_ROW("btm_dp", bottom_row_pair.first, depth);
  LOG_ROW("btm_n ", bottom_row_pair.second, depth);

  // Step 3: Recursion
  for (auto&& m_vector : middle_generator(s, t, top_row_pair, bottom_row_pair, depth)) {
    auto&& [s_init_view, t_init_view] = GetUpperLeft(s, t, m_vector);
    auto&& [s_tail_view, t_tail_view] = GetLowerRight(s, t, m_vector);
    LOG_PAIRS("middle_generator", m_vector, depth);
    LOG_SV("init", s_init_view, t_init_view, depth);
    LOG_SV("tail", s_tail_view, t_tail_view, depth);
    for (PairVec&& l_vector : hirschberg(s_init_view, t_init_view, depth+1)) {
      for (PairVec&& r_vector : hirschberg(s_tail_view, t_tail_view, depth+1)) {
        LOG_PAIRS("l_vector", l_vector, depth);
        LOG_PAIRS("m_vector", m_vector, depth);
        LOG_PAIRS("r_vector", r_vector, depth);
        // const bool l_swap = s_init_view.length() < t_init_view.length();
        // const bool r_swap = s_tail_view.length() < t_tail_view.length();
        PairVec result = Merge(s, t, l_vector, m_vector, r_vector, sv_swap, sv_swap, sv_swap);
        LOG_PAIRS("co_yield", result, depth);
        LOG_D0_MSG("-------------------------------------------------- \n", depth);
        co_yield result;
      }
    }
  } // for each in middle_generator
  LOG_D0_MSG("FINISHED \n\n ", depth);

}

//==============================================================================
// Implementation of Logging Helpers
//==============================================================================
void log_pairs(const std::string_view msg, const PairVec& vec, const int depth) {
  std::cout << std::string(depth * 2, ' ') << msg << ": {";
  bool first = true;
  for (auto elem : vec) {
    if (!first) std::cout << ",| ";
    std::cout << std::format("{} {}", elem.first, elem.second);
    first = false;
  }
  std::cout << "}" << std::endl;
}

void log_sv(const std::string_view msg, const StringView s, const StringView t, const int depth) {
  std::cout << std::string(depth * 2, ' ') << msg << " ";
  std::cout << (!s.empty() ? util::to_string(s) : "\"\"") << " ";
  std::cout << (!t.empty() ? util::to_string(t) : "\"\"") << std::endl;
}

void log_filter_msg(const std::string_view msg, const int lvl, const int depth) {
  if (lvl == depth) {
    std::cout << std::string(depth * 2, ' ') << msg << std::endl;
  }
}

template <typename T>
void log_vector(const std::string_view msg, const std::vector<T>& vec, const int depth) {
  std::cout << std::string(depth * 2, ' ') << msg << " ";
  for (const auto& elem : vec) {
    std::cout << std::format("{}", elem) << " ";
  }
  std::cout << std::endl;
}

//==============================================================================
// Implementation of Index Helpers
//==============================================================================

/**
 * @brief Tests if a one-based pair is a matching position in two strings
 * @param pair One-based Index Position (x,y)
 * @param s StringView for the first string
 * @param t StringView for the second string
 * @return bool s[pair.first -1 ] != t[pair.second -1]
 */
bool IsMatch(const Pair& pair, const StringView s, const StringView t) {
  // assert(s.length() >= t.length() && "Logic: Swapped Propagation");
  return s[pair.first -1 ] == t[pair.second -1];
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
std::generator<Pair> row_major_rect_inc(const Pair top_left, const Pair bottom_right) {
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

//==============================================================================
// Implementation Auxiliary Generators
//==============================================================================
std::generator<Pair> GenTopPairs(const StringView s, const StringView t, const uint row, const uint init_col, const uint n) {
  uint i = row, found = 0;
  while (found < n) {
    const Pair pair{i, init_col};
    if (IsMatch(pair, s, t)) {
      ++found;
      co_yield pair;
    }
    assert(i > 0);
    --i;
  }
}

std::generator<Pair> GenBtmPairs(const StringView s, const StringView t, const uint row, const uint init_col, const uint n) {
  uint i = row, found = 0;
  while (found < n) {
    const Pair pair{i, init_col};
    if (IsMatch(pair, s, t)) {
      ++found;
      co_yield pair;
    }
    assert(i <= s.length());
    ++i;
  }
}

std::generator<OptIdxPair> GenerateKs(const Row& dp_top,
                                      const Row& dp_btm,
                                      const Row& pos_top,
                                      const Row& pos_btm,
                                      uint max_value) {
  // Simple Cases
  if (pos_top.empty() && pos_btm.empty())
    co_return;
  if (pos_top.empty()) {
    for (const auto& elem_btm : pos_btm)
      co_yield {std::nullopt, elem_btm};
    co_return;
  }
  if (pos_btm.empty()) {
    for (const auto& elem_top : pos_top)
      co_yield {elem_top, std::nullopt};
    co_return;
  }

  // Now we have !pos_top.empty() && !pos_btm.empty()
  // First, let's consider the cases where the top half is not used
  if (auto btm_iter = pos_btm.begin(); dp_btm[*btm_iter] == max_value ) {
    while (btm_iter != pos_btm.end() && dp_btm[*btm_iter] == max_value) {
      co_yield {std::nullopt, *btm_iter};
      ++btm_iter;
    }
  }

  // Second, let's consider the cases where the top and btm halves are used
  auto valid_arrangement = [](auto top_pos_iter, const auto btm_pos_iter) {
    return *top_pos_iter < *btm_pos_iter;
  };
  auto valid_max_sum = [&](auto top_iter,  auto btm_iter) -> bool {
    if (btm_iter == pos_btm.end())
      return false;
    return dp_top[*top_iter] + dp_btm[*btm_iter] == max_value;
  };
  auto has_matching_next = [](const auto& pos_iter, const auto& pos, const auto& dp) {
    auto next = std::next(pos_iter);
    if (next == pos.end())
      return false;
    return dp[*pos_iter] == dp[*next];
  };

  ptrdiff_t diff = 0;
  for (auto top_iter = pos_top.begin(); top_iter != pos_top.end(); ++top_iter) {
    auto btm_iter = std::next(pos_btm.begin(), diff);
    while (btm_iter != pos_btm.end() && !valid_arrangement(top_iter, btm_iter))
      ++btm_iter;
    while (valid_max_sum(top_iter, btm_iter)) {
      co_yield {*top_iter, *btm_iter};
      ++btm_iter;
    }
    if (!has_matching_next(top_iter, pos_top, dp_top))
      diff = std::distance(pos_btm.begin(), btm_iter); // continue

    // Third, let's consider the cases where the btm half is not used
    if (dp_top[*top_iter] == max_value) {
      co_yield {*top_iter, std::nullopt};
    }
  }
}

/**
 * Empty Generator
 * @return Generator that yields std::vector<std::pair<uint, uint>>&&
 */
LCS2_STD_H::Generator gen_empty() {
  co_return;
}

/**
 * Handles the base case generation in the main algorithm
 * @param s Longer StringView in the problem
 * @param t Shorter StringView in the problem
 * @param swapped True, iff s and t were swapped to sort them s.t. |s| >= |t|
 * @param depth Current depth in the recursion. Only used for debug output.
 * @return Generator with of type std::vector<std::pair<uint, uint>>
 * - If s.empty and t.empty, one yield: {PairVec{}}
 * - Else if !swapped, multiple yields: `{(k,1) : s[k-1]=t[0]}`
 *        if swapped,  multiple yields: `{(1,k) : s[k-1]=t[0]}`
 * @pre t.length() < 2
 * @pre t.length() < s.length()
 */
LCS2_STD_H::Generator base_cases(const StringView s, const StringView t, const bool swapped, const int depth) {
  // Already swapped in caller. Remember the invariance s.length() >= t.length()
  if (t.empty()) {
    LOG_PAIRS("co_yield", PairVec{}, depth);
    co_yield {};  // Empty Vector for concatenation in the later recursion
    co_return;
  }
  if (t.length() == 1) {
    auto matches_pos = std::views::iota(0u, s.length())
      | std::views::filter([&s, &t](auto i){return s[i]==t[0];})
      | std::views::transform([](auto i){return i+1;});
    const auto b = matches_pos.begin();
    const auto e = matches_pos.end();
    if (b == e) {
      LOG_PAIRS("co_yield", PairVec{}, depth);
      co_yield {};
    } else {
      for (auto it = b; it != e; ++it) {
        auto base_sol = swapped ? PairVec{Pair(1,*it)} : PairVec{Pair(*it,1)};
        LOG_PAIRS("co_yield", base_sol, depth);
        co_yield base_sol;
      }
    }
    co_return;
  }
  assert(false && "base case management should lead to co_return");
}

LCS2_STD_H::Generator middle_generator(StringView s, StringView t, const RowPair& dp_top_pair, const RowPair& dp_bottom_pair, int depth) {
  const uint l1 = s.size();
  const uint l2 = t.size();
  const uint mid_row = (l1 + 1)/2;
  auto && [dp_top, n_top] = dp_top_pair;
  auto && [dp_btm, n_btm] = dp_bottom_pair;
  assert(dp_top.size() == l2 + 2 && "Logic: dp rows should have zero padding");
  assert(dp_btm.size() == l2 + 2 && "Logic: dp rows should have zero padding");
  assert(n_top.size() == l2 + 2 && "Logic: dp rows should have zero padding");
  assert(n_btm.size() == l2 + 2 && "Logic: dp rows should have zero padding");

  // Calculate maximum between connected upper and lower half dp table
  auto sum = std::ranges::to<std::vector<uint>>(std::views::zip_transform(
    [](const uint a, const uint b) {return a+b;},     // Add cells diagonally
    std::views::counted(dp_top.begin()+0, static_cast<int>(l2+1)),  // Do not skip empty represent (dp_top[0]=0)
    std::views::counted(dp_btm.begin()+1, static_cast<int>(l2+1))   // Do not skip empty represent (dp_top[l2+1]=0)
  ));
  const auto first_max_iter = std::ranges::max_element(sum);
  const auto last_max_iter = std::ranges::max_element(sum|std::views::reverse);
  const uint max_value = *first_max_iter;
  const uint first_max_pos = std::distance(sum.begin(), first_max_iter); // O(1) for std::vector
  const uint last_max_pos = sum.size() - 1 - std::distance(sum.rbegin(), last_max_iter); // zero-based right to left idx in `sum`

  auto top_row_pairs = std::ranges::to<std::vector<std::pair<uint,uint>>>(
    std::views::iota(1u) | std::views::take(l2)
      | std::views::filter([&sum, &max_value](const int j) {return sum[j-1] == max_value;})
      | std::views::filter([&s, &t, &mid_row](const uint j) {return s[mid_row-1] == t[j-1];})
      | std::views::transform([&mid_row](const uint j){return std::make_pair(mid_row, j); })
  );
  auto bottom_row_pairs = std::ranges::to<std::vector<std::pair<uint,uint>>>(
    std::views::iota(1u) | std::views::take(l2)
      | std::views::filter([&sum, &max_value](const int j) {return sum[j-1] == max_value;})
      | std::views::filter([&s, &t, &mid_row](const uint j) {return s[mid_row] == t[j];})
      | std::views::transform([&mid_row](const uint j){return std::make_pair(mid_row+1, j+1); })
  );

  const auto top_n_pos = std::ranges::to<std::vector<uint>>(
    std::views::iota(first_max_pos, last_max_pos + 1)
      | std::views::filter([&n_top](const int j) {return n_top[j] > 0;})
  );

  const auto btm_n_pos = std::ranges::to<std::vector<uint>>(
    std::views::iota(first_max_pos, last_max_pos + 1)
      | std::views::transform([](int i){return i + 1;}) // align index with dp_btm
      | std::views::filter([&n_btm](int j) { return n_btm[j] > 0; })
  );

  // Base Case
  if (top_n_pos.empty() && btm_n_pos.empty()) { // This logic is better for debugging than looking in dp
    co_yield {};
    co_return;
  }

  for (const auto&& [k_top_opt, k_btm_opt] : GenerateKs(dp_top,dp_btm, top_n_pos, btm_n_pos, max_value)) {
    if (!k_top_opt.has_value() && !k_btm_opt.has_value()) {
      co_yield {};
    }
    if (k_top_opt.has_value() && !k_btm_opt.has_value()) {
      for (auto&& top_pair : GenTopPairs(s,t, mid_row, *k_top_opt, n_top[*k_top_opt])) {
        co_yield {top_pair};
      }
    }
    if (!k_top_opt.has_value() && k_btm_opt.has_value()) {
      for (auto&& btm_pair : GenBtmPairs(s,t, mid_row+1, *k_btm_opt,n_btm[*k_btm_opt])) {
        co_yield {btm_pair};
      }
    }
    if (k_top_opt.has_value() && k_btm_opt.has_value()) {
      for (auto&& top_pair : GenTopPairs(s,t, mid_row, *k_top_opt, n_top[*k_top_opt])) {
        for (auto&& btm_pair : GenBtmPairs(s,t, mid_row+1, *k_btm_opt,n_btm[*k_btm_opt])) {
          co_yield {top_pair, btm_pair};
        }
      }
    }
  } // loop k pair
}

//==============================================================================
// Implementation DP Logic
//==============================================================================

Row dp_hirsch_top(StringView s, StringView t) {
  const uint l1 = s.length();
  const uint l2 = t.length();
  const uint mid = (l1 + 1) / 2;
  Row dp(l2 + 2, 0), current(l2 + 2, 0);

  // Classic llcs dp. Zero Padding: current[0] and current[l2+1] stay zero
  std::ranges::for_each(
    row_major_rect_inc({1, 1}, {mid, l2}),
    [&](const auto& index_pair) {
      const auto& [i, j] = index_pair;
      current[j] = IsMatch(index_pair, s, t)
        ? dp[j - 1] + 1
        : std::max(dp[j - 1], current[j - 1]);
      if (j == l2)
        std::swap(dp, current);  // set full row
  });
  return dp;
}

Row dp_hirsch_bottom(StringView s, StringView t) {
  const uint l1 = s.length();
  const uint l2 = t.length();
  const uint mid = (l1 + 1) / 2;
  Row dp(l2 + 2, 0), current(l2 + 2, 0);

  // Reverse llcs dp. Zero Padding: current[0] and current[l2+1] stay zero
  std::ranges::for_each(
    row_major_rect_dec({mid + 1, 1}, {l1, l2}),
    [&](const auto& index_pair) {
      const auto& [i, j] = index_pair;
      current[j] = IsMatch(index_pair, s, t)
        ? dp[j + 1] + 1
        : std::max(dp[j], current[j + 1]);
      if (j == 1)
        std::swap(dp, current);
  });
  return dp;
}

/**
 * @brief Extension of the dp lcs hirschberg algorithm to allow reconstruction.
 * @details Goes through the upper part of the classic lcs dp matrix until the
 *          (s.length()+1)/2-th row and returns it. Linear Space Complexity.
 * @param s StringView for the first string
 * @param t StringView for the second string
 * @return RowPair {row, n}. Such that (with s and t one-based):
 * - row[i] = dp[mid][j] with dp the classic lcs dp matrix and mid = (|s|+1)/2
 * - n[j] = #{(k,i)| k in [1:mid] and s[k]==t[j] and dp[mid][j]==dp[k][j]}
 *

 */
RowPair dp_extended_top(StringView s, StringView t) {
  const uint l1 = s.length();
  const uint l2 = t.length();
  const uint mid = (l1 + 1) / 2;
  Row dp(l2 + 2, 0), dp_curr(l2 + 2, 0);
  Row n(l2 + 2, 0);

  // Classic llcs dp. Zero Padding: current[0] and current[l2+1] stay zero
  std::ranges::for_each(
    row_major_rect_inc({1, 1}, {mid, l2}),
    [&](const auto& index_pair) {
      const auto& [i, j] = index_pair;
      const bool match = IsMatch(index_pair, s, t);
      dp_curr[j] = match
        ? dp[j - 1] + 1
        : std::max(dp[j], dp_curr[j - 1]);
      const bool contributes = dp[j] == dp_curr[j];
      n[j] = contributes
        ? (match ? n[j] + 1 : n[j])
        : match ? 1 : 0;
      if (j == l2) {
        std::swap(dp, dp_curr);
      }
  });
  return {dp,n};
}

/**
 * @brief Extension of the dp lcs hirschberg algorithm to allow reconstruction.
 * @details Goes through the lower part of the classic lcs dp matrix until the
 *          (s.length()+1)/2+1 th row and returns it. Linear Space Complexity.
 * @param s StringView for the first string
 * @param t StringView for the second string
 * @complexity
 * @return RowPair {row, n}. Such that (with s and t one-based):
 * - row[j] = dp[mid+1][j] with dp the reversed lcs dp matrix and mid=(|s|+1)/2
 * - n[j] = #{(k,j)| k in [mid+1:|t|] and s[k]==t[j] and dp[mid+1][j]==dp[k][j]}
 */
RowPair dp_extended_bottom(StringView s, StringView t) {
  const uint l1 = s.length();
  const uint l2 = t.length();
  const uint mid = (l1 + 1) / 2;
  Row dp(l2 + 2, 0), dp_curr(l2 + 2, 0);
  Row n(l2 + 2, 0);

  // Reverse llcs dp. Zero Padding: current[0] and current[l2+1] stay zero
  std::ranges::for_each(
    row_major_rect_dec({mid + 1, 1}, {l1, l2}),
    [&](const auto& index_pair) {
      const auto& [i, j] = index_pair;
      const bool match = IsMatch(index_pair, s, t);
      dp_curr[j] = match
        ? dp[j + 1] + 1
        : std::max(dp[j], dp_curr[j + 1]);
      const bool contributes = dp[j] == dp_curr[j];
      n[j] = contributes
        ? (match ? n[j] + 1 : n[j])
        : match ? 1 : 0;
      if (j == 1)
        std::swap(dp, dp_curr);
  });
  return {dp, n};
}

//==============================================================================
// Implementation Auxiliary Functions
//==============================================================================

StrViewPair GetUpperLeft(const StringView s, const StringView t, const PairVec& middle) {
  const uint l1 = s.length();
  const uint l2 = t.length();
  const uint mid_row = (l1 + 1) / 2;
  const uint mid_col = (l2 + 1) / 2;
  const StringView s_init_view = middle.empty()
    ? s.substr(0, mid_row)
    : s.substr(0, middle.front().first - 1);
  const StringView t_init_view = middle.empty()
    ? t.substr(0, mid_col)                          // include t[mid_col]
    : t.substr(0, middle.front().second - 1); // do not include last lcs position, convert to zero-based
  return {s_init_view, t_init_view};
}

StrViewPair GetLowerRight(StringView s, StringView t, const PairVec& middle) {
  const uint l1 = s.length();
  const uint l2 = t.length();
  const uint mid_row = (l1 + 1) / 2;
  const uint mid_col = (l2 + 1) / 2;
  // const uint s_offset = mid_row + 2;  // one-based index position (ignore `dp[mid_row:mid_row+1][]` )
  const uint s_offset = middle.empty()
    ? mid_row + 1 // dp[][mid_col] will be part of the first half
    : middle.back().first + 1;  // one-based index position (ignore `dp[][middle.back().second]`)
  const StringView s_tail_view = s_offset <= l1
    ? s.substr(s_offset - 1)
    : StringView();  // minus 1 bc zero-based
  const uint t_offset = middle.empty()
    ? mid_col + 1 // dp[][mid_col] will be part of the first half
    : middle.back().second + 1;  // one-based index position (ignore `dp[][middle.back().second]`)
  const StringView t_tail_view = t_offset <= l2
    ? t.substr(t_offset - 1)
    : StringView();
  return {s_tail_view, t_tail_view};
}

PairVec Merge(const StringView s, const StringView t, const PairVec& l, const PairVec& m, const PairVec& r, bool swap_l, bool swap_m, bool swap_r) {
  const uint l1 = s.length();
  const uint l2 = t.length();
  const uint mid_row = (l1 + 1) / 2;
  const uint mid_col = (l2 + 1) / 2;
  const uint s_offset = m.empty() // Number of positions before s_tail_view
    ? mid_row
    : m.back().first;
  const uint t_offset = m.empty() // Number of positions before t_tail_view
    ? mid_col
    : m.back().second;

  PairVec result;
  result.reserve(l.size() + m.size() + r.size());
  if (!swap_l) {
    std::ranges::move(l, std::back_inserter(result));
  }
  else {
    std::ranges::transform(l, std::back_inserter(result),
    [](const auto& p) {
      return std::pair{p.second, p.first};
    });
  }
  if (!swap_m) {
    std::ranges::move(m, std::back_inserter(result));
  }
  else {
    std::ranges::transform(m, std::back_inserter(result),
      [](const auto& p) {return std::pair{p.second, p.first};}
    );
  }
  if (!swap_r) {
    std::ranges::transform(r, std::back_inserter(result),
    [s_offset, t_offset](const auto& p) {
      return std::pair{p.first + s_offset  , p.second + t_offset};
    });
  }else {
    std::ranges::transform(r, std::back_inserter(result),
    [s_offset, t_offset](const auto& p) {
      return std::pair{p.second + t_offset, p.first + s_offset };
    });
  }
  return result;
}

//==============================================================================
// Implementation of algorithm class LCS2_STD_H
//==============================================================================

/**
 * Constructor for LCS2_STD_H
 * @param spv The vector of shared points to the strings of a problem
 * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
 */
LCS2_STD_H::LCS2_STD_H(const StrPtrVector &spv, const ConstraintMap &map)
    : BaseAlgorithm(AlgoCategory::LLCS, AlgoType::LCS2_STD_H, spv, map) { }

/**
 * isValid
 * @return true, iff LCS2_STD_H is constructed with two strings
 */
bool LCS2_STD_H::isValid() const {
  return getStrPtrVec().size() == 2;
}

/**
 * Getter for the pseudocode
 * @return std::string_view of the description string
 */
std::string_view LCS2_STD_H::getDescription() const {
  static constexpr std::string_view msg = R"DESC(Pseudocode: LCS Hirschberg
> Function hirsch(s1, s2) -> vector<pair<uint,uint>>:
>   Let length(si) and li<=lj for all i in {1,2}
>   Let mid = (l1 + 1) / 2
>   calc top_row    = dp[mid] via the folklore dp algorithm. Restrict memory to the two rows needed for the calculation.
>   calc bottom_row = dp[mid+1] similar to the folklore dp algorithm, but beginning in the bottom right corner. Restrict
>                     memory to the two rows needed for the calculation.
>   Let k = max_k( top_row[k] + bottom_row[k+1] ) = llcs(s1,s2)
>   Let left = hirsch( s1[1:mid], s2[1:k])
>   Let right = hirsch( s1[mid+1, l1], s2[k+1:l2])
>   For Each element (x,y) in right: set x to x+mid and set y to y + k
)DESC";
  return msg;
}

/**
 * DebugString
 * @return std::string containing the strings of the problem
 */
std::string LCS2_STD_H::DebugString() const {
  std::ostringstream oss;
  oss << toString(s) << "\n";
  return oss.str();
}

/**
 * @brief Executes algo and creates an appropriate LCS RangeTree Solution
 * @return std::unique_ptr<BaseSolution> A unique pointer to a solution object:
 * - `EmptySolution` if the object is not valid or zero strings are provided
 * - `LCS_STD_FL::GenSolution` to generate all possible LCS Sequences (
 */
std::unique_ptr<BaseSolution> LCS2_STD_H::query() {
  reset(ResetLevel::Full);
  if (!isValid()) {
    return std::make_unique<solutions::EmptySolution>();
  }
  if (getState() != State::Preprocessed)
    doPreprocessing();
  return std::make_unique<GenSolution>(getStrPtrVec());
}

/**
 * @brief Does the preprocessing of LCS2_STD_H
 * @details Because the algo is functionally impl., there isn't preprocessing.
 */
void LCS2_STD_H::doPreprocessing() {
  setState(State::Preprocessed);
}

/**
 * Resets the algo that is used to get the lcs points in the GenSolution
 * @param lvl BaseAlgorithm::ResetLevel to be applied
 */
void LCS2_STD_H::reset(const ResetLevel lvl) {
  setState(State::Constructed); //
}


//==============================================================================
// Implementation of Solution class GenSolution
//==============================================================================

/**
 * GenSolution - Constructor
 * @param spv Vector of pointers for the problem's strings
 */
LCS2_STD_H::GenSolution::GenSolution(util::StrPtrVector spv)
    : v_(std::move(spv)) {}

/**
 * GenSolution::begin
 * @return Iterator for returning all possible LCS Subsequences
 */
LCS2_STD_H::GenSolution::AnyIterator LCS2_STD_H::GenSolution::begin() const {
  auto gen = hirschberg(*v_[0],*v_[1]);
  return AnyIterator(GenIterator(v_, std::move(gen)));
}

/**
 * GenSolution::end
 * @return Iterator for returning all possible LCS Subsequences
 */
LCS2_STD_H::GenSolution::AnyIterator LCS2_STD_H::GenSolution::end() const {
  return AnyIterator(GenIterator());
}

/**
 * GenSolution::clone
 * @details Because coroutine states are not copyable, a new generator is
 *          created from the strings `s_` and `t_`
 * @return cloned GenSolution (can be created by calling the constructor)
 */
BaseSolution *LCS2_STD_H::GenSolution::clone() const {
  return new GenSolution(v_);
}

//==============================================================================
// Implementation of Solution Iterator class GenIterator
//==============================================================================

util::StrPtrVector LCS2_STD_H::GenSolution::GenIterator::empty_spv_ = {};

/**
 * @brief Default Constructor for GenIterator
 * @details The default constructor is used to generate the end for an iterator
 */
LCS2_STD_H::GenSolution::GenIterator::GenIterator()
    : generator_(gen_empty()),
      current_(std::make_unique<solutions::Points>()),
      iter_(generator_.begin()),
      spv_(empty_spv_)
{
  current = current_.get(); // Set pointer to the current val in BaseIterator
  is_end = true;            // Set end flag
}

/**
 * @brief Constructor for GenIterator
 * @param spv Vector of pointers for the problem's strings
 * @param gen Generator for embeddings
 */
LCS2_STD_H::GenSolution::GenIterator::GenIterator(
    const util::StrPtrVector& spv, Generator&& gen)
      : generator_(std::move(gen)),
        current_(std::make_unique<solutions::Points>(spv)),
        iter_(generator_.begin()),
        spv_(spv)
{
  current = current_.get(); // Set pointer to the current val in BaseIterator
  is_end = iter_ == std::default_sentinel;  // Set end flag
  advance();
}

/**
 * @brief Goes through the wrapped generator and updates the iterator
 * @details The value type of the `GenIterator` is defined in `BaseIterator` as
 *          `solutions::Points`, which an embeddings as a vector of tuples
 */
void LCS2_STD_H::GenSolution::GenIterator::advance() {
  if (iter_ == std::default_sentinel) {
    is_end = true;
    current_->clear();
    return;
  }
  current_->clear();
  current_->set_embedding(*iter_);
  ++iter_;
  ++n_;
}

/**
 * @brief Creates a copy of an Iterator at the current position
 * @details Coroutines do not support copying. So this method basically starts
 *          the hirschberg algorithm from the beginning and then advances the
 *          iterator, until the iterators stored in the instance match.
 * @return std::unique_ptr<solutions::BaseIterator>
 */
std::unique_ptr<solutions::BaseIterator> LCS2_STD_H::GenSolution::GenIterator::clone() const {
  util::Logger::Warning() << "Used LCS2_STD_H::GenSolution::GenIterator::clone()";
  auto gen = hirschberg(*spv_[0],*spv_[1]);
  auto ptr = std::make_unique<GenIterator>(spv_, std::move(gen));
  for (util::uint i = 0; i < n_; ++i) {
    ++ptr->iter_;
  }
  return ptr;
}

}  // namespace lcs_solver::algorithms::lcs