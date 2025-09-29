#include "util/CommonTypes.h"

#include <algorithm>
#include <format>
#include <locale>
#include <ranges>
#include <sstream>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <windows.h>
#endif

namespace lcs_solver::util {
/*******************************************************************************
 * @brief Converts char to std::string
 * @param c char
 * @param asHex If true converts to an uint, otherwise to utf8 std::string
 * @return std::string describing a char
 ******************************************************************************/
std::string to_string(const char &c, const bool asHex) {
  if (asHex) {
    std::ostringstream oss;
    oss << std::hex << static_cast<int>(c);
    return oss.str();
  }
  return std::format("{}", c);
}

/*******************************************************************************
 * @brief Converts std::string to std::string
 * @param s std::string
 * @param asHex If true converts to an uint, otherwise to utf8 std::string
 * @return std::string
 ******************************************************************************/
std::string to_string(const std::string &s, const bool asHex) {
  const std::string_view v = s;
  return to_string(v, asHex);
}

/*******************************************************************************
 * @brief Converts util::string_view to std::string
 * @param v std::string_view
 * @param asHex If true converts to an uint, otherwise to utf8 std::string
 * @return std::string
 ******************************************************************************/
std::string to_string(const std::string_view v, const bool asHex) {
  if (asHex) {
    std::ostringstream oss;
    oss << "[";
    for (auto it = v.begin(); it != v.end(); ++it) {
      oss << to_string(*it, true);
      if (std::next(it) != v.end()) {
        oss << ", ";
      }
    }
    oss << "]";
    return oss.str();
  }
  return {v.data(), v.size()};
}

/*******************************************************************************
 * @brief Converts char pointer to meaningful std::string
 * @param p util::Symbol pointer
 * @param asHex If true converts to hex codes, otherwise to utf8 std::string
 * @return std::string describing the pointer p
 ******************************************************************************/
std::string to_string(const char *p, const bool asHex) {
  if (p == nullptr) {
    return "nullptr";
  }
  return to_string(*p, asHex);
}

/*******************************************************************************
 * @brief Converts std::string pointer to meaningful std::string
 * @param p util::String pointer
 * @param asHex If true converts to hex codes, otherwise to utf8 std::string
 * @return std::string describing the pointer p
 ******************************************************************************/
std::string to_string(const std::string *p, const bool asHex) {
  if (p == nullptr) {
    return "nullptr";
  }
  return to_string(*p, asHex);
}

/*******************************************************************************
 * @brief Converts std::string_view pointer to meaningful std::string
 * @param p std::string_view pointer
 * @param asHex If true converts to hex codes, otherwise to utf8 std::string
 * @return std::string describing the pointer p
 ******************************************************************************/
std::string to_string(const std::string_view *p, const bool asHex) {
  if (p == nullptr) {
    return "nullptr";
  }
  return to_string(*p, asHex);
}

//== 32bit char section ========================================================
/*******************************************************************************
 * @brief Converts char32_t to std::string
 * @param c char
 * @param asHex If true converts to hex codes, otherwise to utf8 std::string
 * @return std::string describing a char32_t
 ******************************************************************************/
std::string to_string(const char32_t &c, const bool asHex) {
  if (asHex) {
    std::ostringstream oss;
    oss << std::hex << static_cast<int>(c);
    return oss.str();
  }
  return to_utf8(c);
}

/*******************************************************************************
 * @brief Converts std::u32string to std::string
 * @param s std::u32string
 * @param asHex If true converts to hex codes, otherwise to utf8 std::string
 * @return std::string describing a std::u32string
 ******************************************************************************/
std::string to_string(const std::u32string &s, const bool asHex) {
  std::u32string_view v = s;
  return to_string(v, asHex);
}

/*******************************************************************************
 * @brief Converts std::u32string_view to std::string
 * @param v std::u32string_view
 * @param asHex If true converts to hex codes, otherwise to utf8 std::string
 * @return std::string describing a util::Symbol
 ******************************************************************************/
std::string to_string(std::u32string_view v, const bool asHex) {
  if (asHex) {
    std::ostringstream oss;
    oss << "[";
    for (auto it = v.begin(); it != v.end(); ++it) {
      oss << to_string(*it, true);
      if (std::next(it) != v.end()) {
        oss << ", ";
      }
    }
    oss << "]";
    return oss.str();
  }
  return to_utf8(v);
}

/*******************************************************************************
 * @brief Converts char32_t pointer to meaningful std::string
 * @param p char32_t pointer
 * @param asHex If true converts to hex codes, otherwise to utf8 std::string
 * @return std::string describing the pointer p
 ******************************************************************************/
std::string to_string(const char32_t *p, const bool asHex) {
  if (p == nullptr) {
    return "nullptr";
  }
  return to_string(*p, asHex);
}

/*******************************************************************************
 * @brief Converts std::string pointer to meaningful std::string
 * @param p std::u32string pointer
 * @param asHex If true converts to hex codes, otherwise to utf8 std::string
 * @return std::string describing the pointer p
 ******************************************************************************/
std::string to_string(const std::u32string *p, const bool asHex) {
  if (p == nullptr) {
    return "nullptr";
  }
  return to_string(*p, asHex);
}

/*******************************************************************************
 * @brief Converts std::u32string_view pointer to meaningful std::string
 * @param p std::u32string_view pointer
 * @param asHex If true converts to hex codes, otherwise to utf8 std::string
 * @return std::string describing the pointer p
 ******************************************************************************/
std::string to_string(const std::u32string_view *p, const bool asHex) {
  if (p == nullptr) {
    return "nullptr";
  }
  return to_string(*p, asHex);
}

/*******************************************************************************
 * @brief Stream insertion operator for UTF-32 character
 *
 * Allows sending a UTF-32 character (char32_t) directly to an output stream.
 * Converts the character to string representation using util::to_string().
 *
 * @param os Output stream to insert the character into
 * @param c UTF-32 character to insert
 * @return Reference to the output stream for chaining
 ******************************************************************************/
std::ostream& operator<<(std::ostream& os, const char32_t& c) {
 return os << util::to_string(c);
}

/*******************************************************************************
 * @brief Stream insertion operator for UTF-32 string
 *
 * Allows sending a UTF-32 string (std::u32string) directly to an output stream.
 * Converts the string to a printable representation using util::to_string().
 *
 * @param os Output stream to insert the string into
 * @param s UTF-32 string to insert
 * @return Reference to the output stream for chaining
 ******************************************************************************/
std::ostream& operator<<(std::ostream& os, const std::u32string& s) {
 return os << util::to_string(s);
}

/*******************************************************************************
 * @brief Stream insertion operator for UTF-32 string view
 *
 * Allows sending a UTF-32 string view (std::u32string_view) directly to an output stream.
 * Converts the string view to a printable representation using util::to_string().
 *
 * @param os Output stream to insert the string view into
 * @param sv UTF-32 string view to insert
 * @return Reference to the output stream for chaining
 ******************************************************************************/
std::ostream& operator<<(std::ostream& os, const std::u32string_view& sv) {
 return os << util::to_string(sv);
}

/*******************************************************************************
 * @brief Stream extraction operator for UTF-32 character
 *
 * Allows reading a UTF-32 character (char32_t) directly from an input stream.
 * Reads a string representation and converts it to char32_t.
 *
 * @param is Input stream to extract the character from
 * @param c UTF-32 character to store the extracted value
 * @return Reference to the input stream for chaining
 ******************************************************************************/
std::istream& operator>>(std::istream& is, char32_t& c) {
  std::string temp;
  is >> temp;
  c = util::char32_from_utf8(temp);
  return is;
}

/*******************************************************************************
 * @brief Stream extraction operator for UTF-32 string
 *
 * Allows reading a UTF-32 string (std::u32string) directly from an input stream.
 * Reads a string representation and converts it to std::u32string.
 *
 * @param is Input stream to extract the string from
 * @param s UTF-32 string to store the extracted value
 * @return Reference to the input stream for chaining
 ******************************************************************************/
std::istream& operator>>(std::istream& is, std::u32string& s) {
  std::string temp;
  is >> temp;
  s = util::char32_from_utf8(temp);
  return is;
}

/*******************************************************************************
 * @brief Converts a UTF-8 encoded string to a single UTF-32 character
 *
 * @details This function takes a UTF-8 encoded std::string_view and converts it
 * to a single UTF-32 character (char32_t). If the input string represents
 * multiple characters, only the first one is returned. If the input is empty or
 * contains invalid UTF-8, appropriate fallback values are returned.
 *
 * @param v UTF-8 encoded string view to convert
 * @return The first character as UTF-32 (char32_t), 0 if input is empty,
 *         or the Unicode replacement character (U+FFFD) if conversion fails
 ******************************************************************************/
char32_t char32_from_utf8(std::string_view v) {
  if (v.empty()) return 0;
  const std::u32string result = u32string_from_utf8(v, true);
  return result.empty() ? 0xFFFD : result[0];
}

/*******************************************************************************
 * @brief Converts a UTF-8 encoded string to a UTF-32 string
 *
 * This function takes a UTF-8 encoded string view and converts it to a UTF-32
 * string. It properly handles multibyte UTF-8 sequences (up to 4 bytes per
 * character) and includes error handling for invalid UTF-8 sequences.
 *
 * @param v UTF-8 encoded string view to convert
 * @param onlyHead If true, only the first character is converted and returned
 * @return UTF-32 string (std::u32string) containing the converted characters.
 *         Invalid UTF-8 sequences are replaced with the Unicode replacement
 *         character (U+FFFD)
 *
 * @details The function processes the input byte by byte, identifying UTF-8 sequences:
 *   - Single byte (0xxxxxxx): ASCII characters
 *   - Two bytes (110xxxxx 10xxxxxx): Characters in the range U+0080 to U+07FF
 *   - Three bytes (1110xxxx 10xxxxxx 10xxxxxx): Characters in the range U+0800 to U+FFFF
 *   - Four bytes (11110xxx 10xxxxxx 10xxxxxx 10xxxxxx): Characters in the range U+10000 to U+10FFFF
 *
 * If an invalid UTF-8 sequence is encountered, the Unicode replacement character is
 * inserted into the result and processing continues from the next byte.
 ******************************************************************************/
std::u32string u32string_from_utf8(std::string_view v, bool onlyHead) {
  std::u32string result;
  result.reserve(v.size());// Reserve space (will be more than needed)

  size_t i = 0;
  while (i < v.size()) {
    constexpr char32_t invalid_code = 0xFFFD;

    // Single byte (ASCII)
    if (const auto b0 = static_cast<unsigned char>(v[i]); (b0 & 0x80) == 0) {
      result.push_back(static_cast<char32_t>(b0));
      i += 1;
    }
    // Two bytes (110xxxxx 10xxxxxx)
    else if ((b0 & 0xE0) == 0xC0) {
      if (i + 1 >= v.size()
          || (static_cast<unsigned char>(v[i + 1]) & 0xC0) != 0x80) {
        result.push_back(invalid_code);
        i += 1;
        continue;
      }

      const auto b1 = static_cast<unsigned char>(v[i + 1]);
      result.push_back(
          ((b0 & 0x1F) << 6) | (b1 & 0x3F));
      i += 2;
    }
    // Three bytes (1110xxxx 10xxxxxx 10xxxxxx)
    else if ((b0 & 0xF0) == 0xE0) {
      if (i + 2 >= v.size()
          || (static_cast<unsigned char>(v[i + 1]) & 0xC0) != 0x80
          || (static_cast<unsigned char>(v[i + 2]) & 0xC0) != 0x80) {
        result.push_back(invalid_code);
        i += 1;
        continue;
      }

      const auto b1 = static_cast<unsigned char>(v[i + 1]);
      const auto b2 = static_cast<unsigned char>(v[i + 2]);
      result.push_back(
          ((b0 & 0x0F) << 12) | ((b1 & 0x3F) << 6) | (b2 & 0x3F));
      i += 3;
    }
    // Four bytes (11110xxx 10xxxxxx 10xxxxxx 10xxxxxx)
    else if ((b0 & 0xF8) == 0xF0) {
      if (i + 3 >= v.size()
          || (static_cast<unsigned char>(v[i + 1]) & 0xC0) != 0x80
          || (static_cast<unsigned char>(v[i + 2]) & 0xC0) != 0x80
          || (static_cast<unsigned char>(v[i + 3]) & 0xC0) != 0x80) {
        result.push_back(invalid_code);
        i += 1;
        continue;
      }
      const auto b1 = static_cast<unsigned char>(v[i + 1]);
      const auto b2 = static_cast<unsigned char>(v[i + 2]);
      const auto b3 = static_cast<unsigned char>(v[i + 3]);
      result.push_back(
          ((b0 & 0x07) << 18) | ((b1 & 0x3F) << 12) | ((b2 & 0x3F) << 6) | (b3 & 0x3F));
      i += 4;
    }
    // Invalid leading byte
    else {
      result.push_back(invalid_code);
      i += 1;
    }
    if (onlyHead) {
      return result;
    }
  }
  return result;
}

/*******************************************************************************
 * @brief Helper function for converting a char32_t into a std::string
 * @param c char32_t must be unicode encoded (e.g. U'x' or 128522 for a smiley)
 * @return std::string
 ******************************************************************************/
std::string to_utf8(const char32_t c) {
  std::string result;
  if (c <= 0x7F) {
    result += static_cast<char>(c);
  } else if (c <= 0x7FF) {
    result += static_cast<char>(0xC0 | (c >> 6));
    result += static_cast<char>(0x80 | (c & 0x3F));
  } else if (c <= 0xFFFF) {
    result += static_cast<char>(0xE0 | (c >> 12));
    result += static_cast<char>(0x80 | ((c >> 6) & 0x3F));
    result += static_cast<char>(0x80 | (c & 0x3F));
  } else {
    result += static_cast<char>(0xF0 | (c >> 18));
    result += static_cast<char>(0x80 | ((c >> 12) & 0x3F));
    result += static_cast<char>(0x80 | ((c >> 6) & 0x3F));
    result += static_cast<char>(0x80 | (c & 0x3F));
  }
  return result;
}

/*******************************************************************************
 * @brief Helper function for converting a char32_t into a std::string
 * @param v std::u32string_view
 * @return std::string
 ******************************************************************************/
std::string to_utf8(const std::u32string_view v) {
  std::ostringstream oss;
  for (const auto &c : v) {
    oss << to_utf8(c);
  }
  return oss.str();
}

std::string to_String(std::string &&v) {
  return v;
}

/*******************************************************************************
 * @brief Configures the console environment to properly display UTF-8 text
 * @details This function performs platform-specific setup to ensure that the
 * console can correctly display UTF-8 encoded characters:
 * - On Windows: Sets the console output code page to UTF-8 (CP_UTF8)
 *   through the Windows API function SetConsoleOutputCP().
 * - On Linux/macOS: No specific configuration is required as these platforms
 *   support UTF-8 by default in most terminal environments.
 *
 * Additional configuration options are included as commented code:
 * - Optional buffering settings for improved performance
 * - Input mode configuration for Windows
 *
 * @note This function should be called at the beginning of the program
 *       before any UTF-8 output is sent to the console.
 ******************************************************************************/
void ConfigureUTF8Console() {
#ifdef _WIN32
  // Windows-specific UTF-8 console configuration
  SetConsoleOutputCP(CP_UTF8);
  // setvbuf(stdout, nullptr, _IOFBF, 1000);

  // Optional: Enable UTF-8 mode for stdin too if needed
  // _setmode(_fileno(stdin), _O_U8TEXT);
#else
  // Linux and macOS already use UTF-8 by default
  // Just set up buffering if desired
  // std::ios_base::sync_with_stdio(false);
  // setvbuf(stdout, nullptr, _IOFBF, 1000);
#endif
}

/*******************************************************************************
 * @brief Hash function to be used in combination with std::unordered_map
 * @param c util::Symbol to hash
 * @return size_t the result of perfectly hashing the symbol c
 ******************************************************************************/
size_t SymbolPerfectHash::operator()(const Symbol &c) const {
  // return std::hash<Unsigned>()(static_cast<Unsigned>(c));
  return static_cast<uint>(c);
}

/*******************************************************************************
 * @brief Compare function to be used in combination with std::unordered_map
 * @param lhs left operand to the equal operator
 * @param rhs right operand to the equal operator
 * @return bool set to whether the symbols lhs and rhs are equal, or not
 ******************************************************************************/
bool SymbolEqual::operator()(const Symbol &lhs, const Symbol &rhs) const {
  return lhs == rhs;
}

}// namespace lcs_solver::util