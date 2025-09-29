#ifndef LCS_SOLVER_UTIL_STRPTRVECTOR_H_
#define LCS_SOLVER_UTIL_STRPTRVECTOR_H_

#include <cassert>
#include <concepts>
#include <filesystem>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <vector>

#include "util/CommonTypes.h"

namespace lcs_solver::util {

// // Aliases for clarity
using StrPtr = std::shared_ptr<const ::lcs_solver::util::String>;
using StrPtrVector = std::vector<StrPtr>;
using StringViewVector = std::vector<::lcs_solver::util::StringView>;

// Definition of concepts
template<typename Container>
concept EmplaceBackableContainer = requires(Container c, typename Container::value_type v) {
  c.emplace_back(v);
};

template<typename Container>
concept InsertableContainer = requires(Container c, typename Container::value_type v) {
  c.insert(c.end(), v);
};

/*******************************************************************************
 * StrPtrVectorFrom(Strs... strs) is Factory Functions
 * @tparam Strs Variadic template parameter pack where each type must be
 *              convertible to std::string
 * @param strs trs Variable number of string-convertible arguments that will be
 *             converted to shared pointers of String objects
 * @return StrPtrVector A vector containing shared pointers to String objects
 *         created from the input arguments
 ******************************************************************************/
template<typename... Strs>
  requires(std::convertible_to<Strs, std::string> && ...)
StrPtrVector StrPtrVectorFrom(Strs... strs) {
  StrPtrVector spv;
  (spv.emplace_back(std::make_shared<String>(to_String<Symbol>(std::string(strs)))), ...);
  return spv;
}

/*******************************************************************************
 * Converts a StrPtrVector to another container type
 * @tparam Container Container The target container type to convert to
 * @param spv The source StrPtrVector to convert from
 * @return Container A new container of the specified type containing elements
 * converted from the input vector
 * @details This function takes a StrPtrVector (vector of string pointers) and
 * converts it to another container type. The target container must be default-
 * initializable and support either emplace_back or insert operations. Each
 * string pointer in the input vector is converted to the appropriate type for
 * the target container.
 *
 * Type conversions:
 * - If the target container's value_type is util::String, string pointers are dereferenced
 * - If the target container's value_type is std::string, string pointers are dereferenced and converted via util::to_string
 * - For other types, string pointers are dereferenced and static_cast to the target type
 * - Null pointers are handled by creating default-constructed values
 *
 * The function uses reserve() on the target container if available to optimize memory allocation.
 ******************************************************************************/
template<typename Container>
  requires std::default_initializable<Container> && (EmplaceBackableContainer<Container> || InsertableContainer<Container>)
Container StrPtrVectorTo(const StrPtrVector &spv) {
  using T = typename Container::value_type;
  Container result;

  if constexpr (requires(Container c) { c.reserve(0); }) {
    result.reserve(spv.size());
  }

  for (const StrPtr &ptr : spv) {
    T value;
    if constexpr (std::is_same_v<T, util::String>) {
      value = ptr ? *ptr : util::String();
    } else if constexpr (std::is_same_v<T, std::string>) {
      value = ptr ? util::to_string(*ptr) : std::string();
    } else {
      value = ptr ? static_cast<T>(*ptr) : T();
    }

    if constexpr (EmplaceBackableContainer<Container>) {
      result.emplace_back(std::move(value));
    } else {
      result.insert(result.end(), std::move(value));
    }
  }
  return result;
}

StrPtrVector StrPtrVectorFrom(const std::vector<std::string> &vec);
StrPtrVector StrPtrVectorFrom(std::initializer_list<std::string> list);
StringViewVector StringViewVecFrom(const StrPtrVector &);

size_t CalcMinStrLen(const StrPtrVector &spv);
size_t CalcMaxStrLen(const StrPtrVector &spv);
std::vector<Symbol> CalcAlphabet(const std::span<const String> &s);
std::vector<Symbol> CalcAlphabet(const StrPtrVector &spv);
std::pair<std::vector<size_t>, std::vector<size_t>> GenerateSortPerm(StrPtrVector &spv);

/*******************************************************************************
 * @brief Applies a permutation to a container.
 *
 * This function rearranges the elements of the input container according to the
 * specified permutation vector. For each index `i` in the input, the element
 * `input[i]` is moved to position `perm[i]` in the result.
 *
 * @tparam Container A random-access container type supporting `size()`,
 *         `operator[]`, and copy assignment of elements.
 * @param input The input container whose elements are to be permuted.
 * @param perm A vector of indices specifying the target positions for each
 *        element of the input. Must have the same size as the input container,
 *        and all indices must be valid (less than input.size()).
 * @return A new container of the same type as the input, with elements
 *         rearranged according to the permutation.
 * @note
 * - The permutation must be a rearrangement of indices in `[0, input.size())`
 * - This function uses `assert` to do some range indices checks
 *
 * @example
 * @code
 * std::vector<int> input = {10, 20, 30, 40};
 * std::vector<std::size_t> perm = {2, 0, 3, 1};
 * std::vector<int> result = ApplyPermutation(input, perm); // result is {20, 40, 10, 30}
 * @endcode
 ******************************************************************************/
template <typename Container>
Container ApplyPermutation(const Container& input, const std::vector<std::size_t>& perm) {
  assert(input.size() == perm.size() && "Permutation and container must have the same size.");

  Container result(input.size());
  for (std::size_t i = 0; i < perm.size(); ++i) {
    assert(perm[i] < input.size() && "Permutation index out of bounds.");
    result[perm[i]] = input[i];
  }
  return result;
}

StrPtrVector ReadStrPtrVec(
  std::string_view msg = "",
  const StrPtrVector * default_value = nullptr,
  std::istream &in=std::cin,
  std::ostream &os=std::cout
);

StrPtrVector ReadStrPtrVec(
  uint length,
  std::string_view msg = "",
  const StrPtrVector * default_value = nullptr,
  std::istream &in=std::cin,
  std::ostream &os=std::cout
);

}// namespace lcs_solver::util
#endif//LCS_SOLVER_UTIL_STRPTRVECTOR_H_
