#ifndef LCS_SOLVER_UTIL_CANTOR_H_
#define LCS_SOLVER_UTIL_CANTOR_H_

#include <cassert>
#include <concepts>
#include <limits>
#include <type_traits>
#include <utility>

namespace lcs_solver::util {

template <typename T>
concept UnsignedInt = std::is_unsigned_v<T>;

struct Cantor {
 public:
  template<UnsignedInt T>
  static T calc(const std::pair<T,T>& pair) {
    T a = pair.first, b = pair.second;
    T sum = a + b;
    assert(sum >= a);
    T product = sum * (sum + 1);
    assert(product >= sum);
    T result = product / 2 + b;
    assert(result >= b);
    return result;
  }
};

} // namespace lcs_solver::util

#endif //LCS_SOLVER_UTIL_CANTOR_H_
