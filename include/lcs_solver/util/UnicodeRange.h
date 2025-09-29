#ifndef LCS_SOLVER_UTIL_UNICODE_RANGE_H_
#define LCS_SOLVER_UTIL_UNICODE_RANGE_H_

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "util/CommonTypes.h"

namespace lcs_solver::util {
class UnicodeRange {
 public:
  using uint = uint32_t;
  using Pair = std::pair<uint, uint>;
  const uint32_t invalid_code = 0x0000;

  enum class Alphabet {
    kASCIIDigits,
    kBasicLatin,
    kBinary,
    kLatinLower,
    kLatinUpper,
    kLatinLowerUpper,
    kPrintable
  };

  explicit UnicodeRange(std::span<std::pair<std::size_t, std::size_t>> ranges);
  explicit UnicodeRange(Alphabet alph);
  [[nodiscard]] uint32_t operator[](std::size_t i) const;
  [[nodiscard]] size_t size() const;
  [[nodiscard]] std::vector<uint32_t> uintAlphabetVector() const;
  [[nodiscard]] std::vector<Symbol> symbAlphabetVector() const;

 private:
  [[nodiscard]] std::vector<std::tuple<size_t, size_t, size_t>> calcInfoVec() const;
  std::vector<Pair> ranges_;
  std::vector<std::tuple<size_t, size_t, size_t>> info_;
};
}// namespace lcs_solver::util
#endif// LCS_SOLVER_UTIL_UNICODE_RANGE_H_
