/*******************************************************************************
 * @file StringMaker.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Utility functions for `std::vector<std::shared_ptr<String>>` objects
 ******************************************************************************/

#include "util/StrPtrVector.h"

#include <algorithm>
#include <set>

#include "util/InOutHelper.h"

namespace lcs_solver::util {

StrPtrVector StrPtrVectorFrom(const std::vector<std::string> &vec) {
  StrPtrVector spv;
  for (const auto &s : vec) {
    String temp = to_String<Symbol>(s);
    spv.emplace_back(std::make_shared<String>(temp));
  }
  return spv;
}

/*******************************************************************************
 * @brief Construct a vector of shared String pointers from a list of string
 *literals.
 * @details Converts each `std::string` in the initializer list into your
 * internal `String` representation (via `to_String<Symbol>`), wraps it in a
 * `std::shared_ptr`, and collects them in a vector.
 * @param list  An initializer list of `std::string` values to convert.
 * @return A `StrPtrVector` (i.e. `std::vector<std::shared_ptr<String>>`)
 *         containing one shared pointer per input string.
 * @note Each call to `to_String<Symbol>` produces a new `String` instance;
 *       the caller takes ownership via `shared_ptr`.
 * @complexity O(N x M) where N = number of strings, M = average string length.
 ******************************************************************************/
StrPtrVector StrPtrVectorFrom(const std::initializer_list<std::string> list) {
  std::vector<StrPtr> spv;
  for (const std::string &s: list) {
    String temp = to_String<Symbol>(s);
    spv.emplace_back(std::make_shared<String>(temp));
  }
  return spv;
}

/*******************************************************************************
 * @brief Build a vector of StringView objects from shared String pointers.
 * @details For each shared pointer in `spv`, if it is non-null, constructs a
 * `StringView` that refers to the pointed-to `String`; otherwise, appends an
 * empty `StringView`
 * @param spv  A vector of `std::shared_ptr<const String>`
 * @return A `StringViewVector` (i.e. `std::vector<StringView>`) of the same size
 * @note Null pointers in `spv` become default-constructed (empty) `StringView`s
 * @complexity O(N) where N = `spv.size()`.
 ******************************************************************************/
StringViewVector StringViewVecFrom(const StrPtrVector &spv) {
  if (spv.empty()) {
    return {};
  }
  std::vector<StringView> result;
  result.reserve(spv.size());
  for (const auto &ptr : spv) {
    if (ptr) {
      result.emplace_back(*ptr);
    } else {
      result.emplace_back();
    }
  }
  return result;
}

/*******************************************************************************
 * @brief Compute the length of the shortest string in a vector of String pointers.
 * @details Iterates through all non-null pointers in `spv` and returns the
 * minimum string length. If `spv` is empty, returns 0.
 * @param spv  A vector of `std::shared_ptr<const String>`.
 * @return The minimum `String::length()` among all elements, or 0 if `spv` is empty.
 * @complexity O(N) where N = `spv.size()`.
 ******************************************************************************/
size_t CalcMinStrLen(const StrPtrVector &spv) {
  uint min = spv.empty() ? 0 : spv.front()->length();
  for (const auto &p : spv) {
    min = std::min(min, p->length());
  }
  return min;
}

/*******************************************************************************
 * @brief Compute the length of the longest string in a vector of String pointers.
 * @details Iterates through all non-null pointers in `spv` and returns the
 * maximum string length. If all pointers are null or the vector is empty,
 * returns 0.
 * @param spv  A vector of `std::shared_ptr<const String>`.
 * @return The maximum `String::length()` among all elements, or 0 if none are non-null.
 * @complexity O(N) where N = `spv.size()`.
 ******************************************************************************/
size_t CalcMaxStrLen(const StrPtrVector &spv) {
  uint max = 0;
  for (const auto &p : spv) {
    max = std::max(max, p->length());
  }
  return max;
}

std::vector<Symbol> CalcAlphabet(const std::span<const String> &s) {
  std::set<Symbol> alphabet;
  for (const auto &p : s) {
    for (const auto &c : p) {
      alphabet.insert(c);
    }
  }
  return {alphabet.begin(), alphabet.end()};
}

/*******************************************************************************
 * @brief Extract the sorted set of unique symbols (characters) across all
 *strings.
 * @details Scans each `String` in `spv`, inserts each character into a
 * `std::set<Symbol>` to deduplicate and sort, then returns the result as a
 * `std::vector<Symbol>`.
 * @param spv  A vector of `std::shared_ptr<const String>`.
 * @return A sorted `std::vector<Symbol>` containing each distinct character
 *found.
 * @complexity O(N * M * log U) where
 *   - N = number of strings,
 *   - M = average string length,
 *   - U = number of unique symbols.
 ******************************************************************************/
std::vector<Symbol> CalcAlphabet(const StrPtrVector &spv) {
  std::set<Symbol> alphabet;
  for (const auto &p : spv) {
    for (const auto &c : *p) {
      alphabet.insert(c);
    }
  }
  return {alphabet.begin(), alphabet.end()};
}

/*******************************************************************************
 * StrPtrVector factory function that uses streams
 * @param msg Message streamed to os in the first line when calling the function
 * @param default_value Sets the default value of the vector to be generated
 * @param in in stream (default is std::cin)
 * @param os out steam (default is std::cout)
 * @return StrPtrVector read from os
 ******************************************************************************/
StrPtrVector ReadStrPtrVec(const std::string_view msg,
                           const StrPtrVector* default_value,
                           std::istream& in,
                           std::ostream& os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    // os << (msg.empty() ? "Begin with the prompts for the strings.\n" : msg);
    os << (msg.empty() ? "" : msg);
  }

  constexpr std::string_view prompt = "Enter the number of strings";
  const uint n = ReadUnsigned(prompt, 2, in, os);
  return ReadStrPtrVec(n, msg, default_value, in, os);
}

/*******************************************************************************
 * StrPtrVector factory function that uses streams
 * @param length Size of the StrPtrVector to read in
 * @param msg Message streamed to os in the first line when calling the function
 * @param default_value Sets the default value of the vector to be generated
 * @param in in stream (default is std::cin)
 * @param os out steam (default is std::cout)
 * @return StrPtrVector read from os
 ******************************************************************************/
StrPtrVector ReadStrPtrVec(const uint length,
                           std::string_view msg,
                           const StrPtrVector* default_value,
                           std::istream& in,
                           std::ostream& os) {
  StrPtrVector spv;
  spv.reserve(length);
  for (uint i = 0; i < length; ++i) {
    StringView default_view;  // is default constructable (refers to empty str)
    if (default_value && i < default_value->size()) {
      default_view = *(*default_value)[i];
    }
    const String current = ReadUtilString(i, default_view, in, os);
    spv.push_back(std::make_shared<const String>(current));
  }
  return spv;
}

/*******************************************************************************
 * Helper: Calculates permutations for sorting original_spv_ and applies it
 * @details In a problem's constructor original_spv_ is set and not sorted. In
 *  this function sorted_spv_ is sorted by applying the following permutation:
 *  sorted_spv_[perm_[original_index]]:= original_spv_[original_index]
 *  So perm_ maps original indices to sorted indices and perm_r_ maps sorted
 *  indices to original indices. For sorting, std::stable_sort is used.
 ******************************************************************************/
std::pair<std::vector<size_t>, std::vector<size_t>> GenerateSortPerm(const StrPtrVector &spv) {
  std::pair<std::vector<size_t>, std::vector<size_t>> result;
  std::vector<size_t>& perm = result.first;
  std::vector<size_t>& perm_r = result.second;

  const size_t n = spv.size();
  std::vector<std::pair<size_t, size_t>> length_and_spv_idx;
  length_and_spv_idx.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    length_and_spv_idx.emplace_back(spv[i]->size(), i);
  }
  // Sort by length (using pre-computed lengths) and string content
  std::ranges::stable_sort(length_and_spv_idx, [spv](const auto &a, const auto &b) {
    if (a.first != b.first) return a.first < b.first;
    return *spv[a.second] < *spv[b.second];
  });

  // Generate permutation vectors
  perm.resize(n);
  perm_r.resize(n);
  for (size_t sorted_index = 0; sorted_index < length_and_spv_idx.size(); ++sorted_index) {
    const size_t original_index = length_and_spv_idx[sorted_index].second;
    perm[original_index] = sorted_index;
    perm_r[sorted_index] = original_index;
  }

  return result;
}

}// namespace lcs_solver::util