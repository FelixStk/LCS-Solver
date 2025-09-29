#include "util/UnicodeRange.h"

#include <algorithm>
#include <cassert>

namespace lcs_solver::util {

/*******************************************************************************
 * Constructor
 * @param ranges Closed Intervals for Unicode Ranges
 ******************************************************************************/
UnicodeRange::UnicodeRange(std::span<std::pair<std::size_t, std::size_t>> ranges)
    : ranges_(ranges.begin(), ranges.end()), info_(calcInfoVec()) {}

/*******************************************************************************
 * Constructor
 * @param alph UnicodeRange::Alphabet
 ******************************************************************************/
UnicodeRange::UnicodeRange(const Alphabet alph) {
  switch (alph) {
    case Alphabet::kPrintable:
      ranges_.emplace_back(0x0020, 0x007E);  // Basic Latin (ASCII printable characters)
      ranges_.emplace_back(0x00A0, 0x02FF);  // Latin-1 Supplement, Latin Extended, IPA Extensions
      ranges_.emplace_back(0x0370, 0x052F);  // Greek, Cyrillic
      ranges_.emplace_back(0x2000, 0x206F);  // General punctuation
      ranges_.emplace_back(0x2100, 0x218F);  // Letter like symbols, number forms
      ranges_.emplace_back(0x2460, 0x25FF);  // Enclosed alphanumerics, geometric shapes
      ranges_.emplace_back(0x2600, 0x27BF);  // Miscellaneous symbols, Dingbats
      ranges_.emplace_back(0x4E00, 0x9FFF);  // CJK Unified Ideographs
      ranges_.emplace_back(0x1F300, 0x1F6FF);// Emoji symbols
      break;
    case Alphabet::kBasicLatin:
      ranges_.emplace_back(0x0020, 0x007E);// Basic Latin (ASCII printable characters)
      break;
    case Alphabet::kASCIIDigits:
      ranges_.emplace_back(0x0030, 0x0039);// ASCII Digits
      break;
    case Alphabet::kLatinUpper:
      ranges_.emplace_back(0x0041, 0x005A);// ASCII Upper Alphabet
      break;
    case Alphabet::kLatinLower:
      ranges_.emplace_back(0x0061, 0x007A);// ASCII Lower Alphabet
      break;
    case Alphabet::kBinary:
      ranges_.emplace_back(0x0030, 0x0031);
      break;
    case Alphabet::kLatinLowerUpper:
      ranges_.emplace_back(0x0041, 0x005A);// ASCII Upper Alphabet
      ranges_.emplace_back(0x0061, 0x007A);// ASCII Lower Alphabet
      break;
    default:
      assert(true && "Alphabet in switch is missing");
  }
  info_ = std::move(calcInfoVec());
}

/*******************************************************************************
 * operator[]
 * @param i size_t zero-Based position
 * @return Unicode of the ith symbol in a UnicodeRange
 ******************************************************************************/
uint32_t UnicodeRange::operator[](std::size_t i) const {
  if (i >= size()) {
    return invalid_code;
  }

  auto iter = std::ranges::upper_bound(
      info_, i,
      std::less<>{},
      [](const auto &t) { return std::get<2>(t); });
  if (iter != info_.begin()) {
    --iter;
  }
  const auto &[idx, length, start] = *iter;
  return ranges_[idx].first + i - start;
}

/*******************************************************************************
 * @brief Calculates the number of unicodes given specified by an alphabet
 * @detail Uses the calculated info_ vector. So its complexity O(1).
 * @return size_t
 ******************************************************************************/
size_t UnicodeRange::size() const {
  // size_t length = 0;
  // for (const auto &range : ranges_) {
  //   const auto &[l, u] = range;
  //   length += u - l + 1;
  // }
  // return length;
  if (info_.empty()) {
    return 0;
  }
  const auto &[i, len, first] = info_.back();
  return first + len;
}

/*******************************************************************************
 * @brief Generator of a vector of Symbols to represent an Alphabet
 * @return std::vector<Symbol> containing all unicodes specified with ranges_
 ******************************************************************************/
std::vector<uint32_t> UnicodeRange::uintAlphabetVector() const {
  std::vector<uint32_t> result;
  result.reserve(size());
  for (const auto &range : ranges_) {
    const auto &[l, u] = range;
    for (uint32_t i = l; i <= u; ++i) {
      result.emplace_back(i);
    }
  }
  return result;
}

std::vector<Symbol> UnicodeRange::symbAlphabetVector() const {
  std::vector<Symbol> result;
  result.reserve(size());
  for (const auto &range : ranges_) {
    const auto &[l, u] = range;
    for (uint32_t i = l; i <= u; ++i) {
      result.emplace_back(static_cast<Symbol>(i));
    }
  }
  return result;
}

/*******************************************************************************
 * @brief Helper: Calculates a vector<tupes> for searching with std::lower_bound
 * @return std::vector<std::tuple<size_t, size_t, size_t>> info where info[i] =
 *  (i, length, start) such that length is the size of ranges_[i] and start is
 *  the index of the symbol ranges_[i].first with respect in the zero based
 *  alphabet given concatenating ranges_.
 ******************************************************************************/
std::vector<std::tuple<size_t, size_t, size_t>> UnicodeRange::calcInfoVec() const {
  size_t idx = 0;
  size_t length = 0;
  size_t total = 0;
  std::vector<std::tuple<size_t, size_t, size_t>> info;
  info.reserve(ranges_.size());
  for (const auto &range : ranges_) {
    length = range.second - range.first + 1;
    info.emplace_back(idx++, length, total);
    total += length;
  }
  return info;
}

}// namespace lcs_solver::util