#ifndef LCS_SOLVER_UTIL_METAPOW2LOG2_H_
#define LCS_SOLVER_UTIL_METAPOW2LOG2_H_

#include <cstddef>
#include <cassert>
#include <array>
#include <vector>

namespace lcs_solver::util {

template<size_t LEN>
struct Meta_Log2_table {
  static consteval auto gen(){
    std::array<size_t, LEN> a{};
    if constexpr (LEN > 0) a[0] = 0;
    if constexpr (LEN > 1) a[1] = 1;
    for (size_t i = 2; i < LEN; ++i)
      a[i] = a[(i-1) / 2] + 1;
    return a;
  }

  static constexpr std::array<size_t, LEN> arr = gen(); ///< arr = [log2(1), log2(2), ..., log2(LEN)]

  size_t operator[](const size_t i) const{
    assert(i-1 <= LEN); // if i==0 evaluate to false
    return arr[i-1];
  }
};

template<typename T>
struct Log2_table {
  explicit Log2_table(size_t len){
    arr.reserve(len);
    if (len > 0) arr.push_back(0);
    if (len > 1) arr.push_back(1);
    for (size_t i = 2; i < len; ++i)
      arr.push_back(arr[(i-1) / 2] + 1); // a[i] = a[(i-1) / 2] + 1;
  }

  std::vector<T> arr;

  T operator[](const size_t i) const{
    assert(i-1 <= arr.size()); // if i==0 evaluates to false
    return arr[i-1];
  }
};

template<size_t LEN>
struct Meta_Pow2_table {
  static consteval auto gen() {
    std::array<size_t, LEN> a{};
    if constexpr (LEN > 0) a[0] = 1;
    for (size_t i = 1; i < LEN; ++i)
      a[i] = (a[i - 1] << 1);
    return a;
  }

  static constexpr size_t N = std::min<size_t>(LEN, sizeof(size_t) * 8);
  static constexpr std::array<size_t, N> arr = gen(); ///< arr = [2**0, 2**1, ..., 2**LEN]

  size_t operator[](const size_t i) const{
    assert(i <= LEN);
    return arr[i];
  }
};

template<typename T>
struct Pow2_table {
  explicit Pow2_table(size_t len) {
    size_t n = std::min<size_t>(len, sizeof(T) * 8);
    arr.reserve(n);
    if (len > 0) arr.push_back(1);
    for (size_t i = 1; i < len; ++i)
      arr.push_back(arr.back() << 1); //a[i] = (a[i - 1] << 1);
  }

   std::vector<T> arr; ///< arr = [2**0, 2**1, ..., 2**LEN]

  T operator[](const size_t i) const{
    assert(i <= arr.size());
    return arr[i];
  }
};


}
#endif //LCS_SOLVER_UTIL_METAPOW2LOG2_H_