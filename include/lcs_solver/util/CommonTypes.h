#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <cstddef>// size_t
#include <cstdint>
#include <span>
#include <string>
#include <string_view>

namespace lcs_solver::util {
// using uint = uint16_t;
using uint = std::size_t;

// using Symbol = char;
using Symbol = char32_t;
using String = std::basic_string<Symbol>;
using StringView = std::basic_string_view<Symbol>;

std::string to_string(const char &c, bool asHex = false);
std::string to_string(const std::string &s, bool asHex = false);
std::string to_string(std::string_view v, bool asHex = false);
std::string to_string(const char *p, bool asHex = false);
std::string to_string(const std::string *p, bool asHex = false);
std::string to_string(const std::string_view *p, bool asHex = false);

std::string to_string(const char32_t &c, bool asHex = false);
std::string to_string(const std::u32string &s, bool asHex = false);
std::string to_string(std::u32string_view v, bool asHex = false);
std::string to_string(const char32_t *p, bool asHex = false);
std::string to_string(const std::u32string *p, bool asHex = false);
std::string to_string(const std::u32string_view *p, bool asHex = false);

std::ostream& operator<<(std::ostream& os, const char32_t& c);
std::ostream& operator<<(std::ostream& os, const std::u32string& s);
std::ostream& operator<<(std::ostream& os, const std::u32string_view & sv);

std::istream& operator>>(std::istream& is, char32_t& c);
std::istream& operator>>(std::istream& is, std::u32string& s);

void ConfigureUTF8Console();

char32_t char32_from_utf8(std::string_view v);
std::u32string u32string_from_utf8(std::string_view v, bool onlyHead = false);
std::string to_utf8(char32_t c);
std::string to_utf8(std::u32string_view v);

// For casting to util::String
template<typename T>
std::basic_string<T> to_String(std::string &&s) {
  if constexpr (std::is_same<T, char>()) {
    return s;
  } else if constexpr (std::is_same<T, char32_t>()) {
    return u32string_from_utf8(std::string_view(s));
  }
  return std::basic_string<T>();
}

template<typename T>
std::basic_string<T> to_String(const std::string & s) {
  if constexpr (std::is_same<T, char>()) {
    return s;
  } else if constexpr (std::is_same<T, char32_t>()) {
    return u32string_from_utf8(std::string_view(s));
  }
  return std::basic_string<T>();
}

struct SymbolPerfectHash {
  size_t operator()(const Symbol &c) const;
};

struct SymbolEqual {
  bool operator()(const Symbol &lhs, const Symbol &rhs) const;
};

}// namespace lcs_solver::util

#endif// COMMON_TYPES_H