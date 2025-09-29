/**
 * @file RadixSort.h
 * @author Felix Steinkopp
 * @version 1.1
 * @brief Implementation of RadixSort
 */

#ifndef LCS_SOLVER_UTIL_RADIXSORT_H_
#define LCS_SOLVER_UTIL_RADIXSORT_H_

#include <algorithm>
#include <bit>
#include <concepts>
#include <numeric>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace lcs_solver::util {

template<typename T>
concept BitwiseOperable = requires(T a, T b) {
  { a >> 1 } -> std::same_as<T>;
  { a & b } -> std::same_as<T>;
};

/*******************************************************************************
 * Implementation of Radix Sort (with binary digits)
 * @tparam T Typename of the keys
 * @tparam U Typename of the values
 ******************************************************************************/
template<typename T, typename U = std::monostate> requires BitwiseOperable<T>
struct RadixSort {
  using Pair = std::pair<T, U>;
  using KeyValSpan = std::span<const Pair>;
  using KeyValVec = std::vector<Pair>;
  using ReturnPair = std::pair<KeyValVec, std::vector<size_t>>;

 public:
  /*****************************************************************************
   * sort ( key value pairs )
   * @param arr std::span<T,U> where T,U are key and node type
   * @param rmDoubles if true delete duplicates from the sorted result
   * @return Pair: (vec of sorted Pairs, surjection from arr to vec)
   * @note label[i]=j means arr[i] = vec[j], where max(j) <= max(i)
   ****************************************************************************/
  static ReturnPair sort(const KeyValSpan &arr, bool rmDoubles = false) {
    if (arr.empty()) return {{}, {}};
    if (arr.size() == 1) return {{arr[0]}, {0}};

    // Setup variable
    auto vec = KeyValVec(arr.begin(), arr.end());
    auto label = std::vector<size_t>(vec.size());
    auto iterations = getWidth(arr);
    std::iota(std::begin(label),
              std::end(label),
              0); // Fill with [0:l.size()-1]

    // CountSort from the right most position to the left
    for (int i = 0; i < iterations; ++i) {
      countSort(vec, label, i);
    }

    // Calculate mapping from unsorted idx to sorted idx
    std::vector<size_t> outputLabel(label.size());
    if (!rmDoubles) {
      for (size_t i = 0; i < label.size(); ++i) {
        outputLabel[label[i]] = i;
      }
    }
    else{
      // Duplicate Removal (in linear time)
      auto p = [](const Pair &a, const Pair &b) {
        return a.first == b.first;
      };
      auto first = vec.begin(), last = vec.end(), result = first;
      size_t removals = 0;
      size_t pos = 0; // index in vec without removal affecting it
      outputLabel[label[0]] = 0;
      if (!vec.empty()) {
        while (++first != last) {
          ++pos;
          if (!p(*result, *first) && ++result != first) {
            // New district key value pair! Move it forward in vec
            *result = std::move(*first);
          } else if (p(*result, *first) && result != first) {
            // Key Value Pair did not change.
            ++removals;
          }
          outputLabel[label[pos]] = pos - removals;
        }
        ++result;
      }
      vec.erase(result, vec.end());
    }
    label = std::move(outputLabel); // label[i] == j <=> arr[i] == vec[j]


    return {vec, label};
  }

  /*****************************************************************************
   * sort (overload when values are not given)
   * @param arr std::span of keys
   * @param rmDoubles whether to delete duplicates from the sorted result
   * @return Pair: (vec of sorted keys, surjection from arr to vec)
   * @note lable[i]=j means arr[i] = vec[j], where max(j) <= max(i)
   ****************************************************************************/
  static ReturnPair sort(std::span<const T> arr, bool rmDoubles = false) {
    auto data = std::vector<Pair>();
    for (auto key : arr) {
      data.push_back({key, {}});
    }
    return sort(data,rmDoubles);
  }

 private:
  ///< Private constructor to prevent instantiation
  RadixSort() = default;

  /*****************************************************************************
   * getWidth
   * @param arr std::span<T,U> where T,U are key and node type
   * @return the largest bit width of a key in arr
   ****************************************************************************/
  static int getWidth(const KeyValSpan &arr) {
    auto p = std::max_element(arr.begin(),
                              arr.end(),
                              [](const Pair &a, const Pair &b) {
                                return std::bit_width<T>(a.first)
                                    < std::bit_width<T>(b.first);
                              });
    return std::bit_width<T>(p->first);
  }

  /*****************************************************************************
   * countSort
   * @param vec std::vector of key value pairs to sort
   * @param label the permutation that describes the process of sorting
   * @param bit the current bit index used to sort.
   ****************************************************************************/
  static void countSort(KeyValVec &vec, std::vector<size_t> &label, int bit) {
    std::vector<Pair> output(vec.size());
    std::vector<size_t> outputLabel(vec.size());
    size_t count[2] = {0}; // two buckets for binary radix sort

    // Count occurrences of 0 and 1 at position 'bit'
    for (const auto &item : vec) {
      count[(item.first >> bit) & 1]++;
    }

    // Calculate positions
    count[1] += count[0];

    // Build the output array
    for (int i = vec.size() - 1; i >= 0; i--) {
      int bitVal = (vec[i].first >> bit) & 1;
      output[--count[bitVal]] = vec[i];
      outputLabel[count[bitVal]] = label[i];
    }

    vec = std::move(output);
    label = std::move(outputLabel);
  }
};

} // end of namespace

#endif //LCS_SOLVER_UTIL_RADIXSORT_H_