/*******************************************************************************
 * @file InOutHelper.cc
 * @author Steinkopp:Felix
 * @version 2.0
 * @brief Provides utility functions for input/output operations.
 ******************************************************************************/

#include "util/InOutHelper.h"

#include <sstream>
#include <algorithm>
#include <cctype>      // for std::isspace, std::isprint

namespace lcs_solver::util {

/*******************************************************************************
 * Reads a char and prompts with a default value
 * @param msg The message to display. If empty, a default message is used.
 * @param default_value The value returned if input is empty.
 * @param in The input stream (e.g., std::cin).
 * @param os The output stream (e.g., std::cout). Prompts are only displayed if
 * os is `std::cout`.
 * @return std::string The input string or default_value if the input is empty
 ******************************************************************************/
char ReadChar(const std::string_view msg,
              const char default_value,
              std::istream& in,
              std::ostream& os) {
   if (os.rdbuf() == std::cout.rdbuf()) {
     std::string def(1, default_value);
     switch (default_value) {
       case '\'': def = "\\'"; break;
       case '\"': def = "\\\""; break;
       case '?' : def = "\\?"; break;
       case '\\': def = "\\\\"; break;
       case '\a': def = "\\a"; break;
       case '\b': def = "\\b"; break;
       case '\f': def = "\\f"; break;
       case '\n': def = "\\n"; break;
       case '\r': def = "\\r"; break;
       case '\t': def = "\\t"; break;
       case '\v': def = "\\v"; break;
       case '\0': def = "\\0"; break;
       default:
         if (!std::isprint(static_cast<unsigned char>(default_value))) {
           std::ostringstream oss;
           oss << "\\x" << std::hex << std::uppercase
               << static_cast<int>(static_cast<unsigned char>(default_value));
           def = oss.str();
         }
         break;
     }
     os << (msg.empty() ? "Enter a char" : msg)
        << " (default=" << def << "):";
   }
   std::string input;
   std::getline(in, input);
   if (input.empty()) return default_value;
   return input[0];
}

/*******************************************************************************
 * Reads a standard string from the input stream with a prompt and default value
 * @param msg The message to display. If empty, a default message is used.
 * @param default_value The value returned if input is empty.
 * @param in The input stream (e.g., std::cin).
 * @param os The output stream (e.g., std::cout). Prompts are only displayed if
 * os is `std::cout`.
 * @return std::string The input string or default_value if the input is empty
 ******************************************************************************/
std::string ReadStdString(
    const std::string_view msg,
    std::string default_value,
    std::istream &in,
    std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os << (msg.empty() ? "Enter a string" : msg)
       << (default_value.empty() ? ":" : " (default=" + default_value + "):");
  }
  std::string input;
  std::getline(in, input);
  if (input.empty()) return default_value;
  return input;
}

/*******************************************************************************
 * @brief Reads a util::String from the input stream with a custom prompt
 * @param msg The prompt message. Uses a default message if empty
 * @param in The input stream
 * @param os The output stream. Prompts are only displayed if os is std::cout
 * @return String The converted util::String
 ******************************************************************************/
String ReadUtilString(
    const std::string_view msg,
    std::istream &in,
    std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os << (msg.empty() ? "Enter a string" : msg) << ":";
  }
  std::string input;
  std::getline(in, input);
  return util::to_String<Symbol>(std::move(input));
}

/*******************************************************************************
 * @brief Reads a util::String for a specific index (e.g., array element)
 *
 * @param n The index to display in the prompt (e.g., "s[3]:")
 * @param in The input stream
 * @param os The output stream. Prompts are only displayed if os is std::cout
 * @return String The converted util::String
 ******************************************************************************/
String ReadUtilString(const uint n, std::istream &in, std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os << "s[" << n << "]:";
  }
  std::string input;
  std::getline(in, input);
  return util::to_String<Symbol>(std::move(input));
}

/*******************************************************************************
 * @brief Reads a util::String for a specific index (e.g., array element)
 *
 * @param n The index to display in the prompt (e.g., "s[3]:")
 * @param default_view The default value to use when nothing is entered
 * @param in The input stream
 * @param os The output stream. Prompts are only displayed if os is std::cout
 * @return String The converted util::String
 ******************************************************************************/
String ReadUtilString(uint n, StringView default_view, std::istream &in,
                      std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os << "s[" << n << "]:";
  }
  std::string input;
  std::getline(in, input);
  if (input.empty()) return {default_view.data()};
  return util::to_String<Symbol>(std::move(input));
}

/*******************************************************************************
 * @brief Reads an unsigned integer from the input stream with a prompt
 * @param msg The prompt message. Uses a default message if empty
 * @param in The input stream
 * @param os The output stream. Prompts are only displayed if os is std::cout
 * @return `uint` The parsed unsigned integer. Returns 0 if input is empty
 * @throw std::invalid_argument If input is not a valid number
 ******************************************************************************/
uint ReadUnsigned(const std::string_view msg, std::istream &in, std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os << (msg.empty() ? "Enter a unsigned int" : msg) << ":";
  }
  std::string input;
  std::getline(in, input);
  if (input.empty()) return {};
  const uint number = std::stoul(input);
  return number;
}

/*******************************************************************************
 * @brief Reads an unsigned integer with a default value
 * @param msg The prompt message. Uses a default message if empty
 * @param default_value The value returned if input is empty
 * @param in The input stream
 * @param os The output stream. Prompts are only displayed if os is std::cout.
 * @return uint The parsed value or default_value if input is empty
 * @throw std::invalid_argument If input is not a valid number
 ******************************************************************************/
uint ReadUnsigned(
    const std::string_view msg,
    const uint default_value,
    std::istream &in,
    std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os  << (msg.empty() ? "Enter a unsigned int" : msg)
        << " (default=" << default_value << "):";
  }
  std::string input;
  std::getline(in, input);
  if (input.empty()) return default_value;

  const uint result = std::stoul(input);
  return result;
}

/*******************************************************************************
 * @brief Reads an unsigned integer within a specified range
 * @param var_name The name of the variable to display in the prompt
 * @param default_value The value returned if input is empty
 * @param range The valid range [min, max] for the input
 * @param in The input stream
 * @param os The output stream. Prompts are only displayed if os is std::cout.
 * @return uint The parsed value within the range
 * @throw std::runtime_error If input is outside the valid range
 * @throw std::invalid_argument If input is not a valid number
 ******************************************************************************/
uint ReadUnsigned(
    const std::string_view var_name,
    const uint default_value,
    const std::pair<uint, uint> &range,
    std::istream &in,
    std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os  << "Enter " << var_name << (var_name.empty() ? "" : " ")
        << "a number from the interval " << util::to_string(range)
        << " (default=" << default_value << "):";
  }
  std::string input;
  std::getline(in, input);
  if (input.empty()) return default_value;

  const uint result = std::stoul(input);
  if (result < range.first || result > range.second)
    throw std::runtime_error("ReadUnsigned: uint input is not in range " +
                             to_string(range));
  return result;
}

/*******************************************************************************
 * @brief Reads a subinterval [first, second] within a specified range
 *
 * @param var_name The name of the variable to display in the prompt
 * @param default_value The default interval returned if input is empty
 * @param range The valid range that the subinterval must lie within
 * @param in The input stream
 * @param os The output stream. Prompts are only displayed if os is std::cout
 * @return std::pair<uint, uint> The parsed interval
 * @throw std::runtime_error If input is formatted invalid or the inputted
 * interval is out of bounds
 ******************************************************************************/
std::pair<uint, uint> ReadSubinterval(
    const std::string_view var_name,
    std::pair<uint, uint> default_value,
    const std::pair<uint, uint> &range,
    std::istream &in,
    std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os  << "Enter " << var_name << (var_name.empty() ? "" : " ")
        << "a pair of two numbers in the sub-interval " << to_string(range)
        << " (default=" << to_string(default_value) << "):";
  }

  std::string input;
  std::getline(in, input);
  if (input.empty()) return default_value;

  uint first, second;
  if (std::istringstream iss(input); !(iss >> first >> second))
    throw std::runtime_error("ReadSubinterval: pair input had invalid format");
  if (first < range.first || second > range.second || first > second)
    throw std::runtime_error("ReadSubinterval: pair input was not in range " +
                             to_string(range));
  return {first, second};
}

/*******************************************************************************
 * @brief Converts a pair of unsigned integers to a string in the format "[a,b]".
 *
 * @param pair The pair to convert.
 * @return std::string The formatted string.
 ******************************************************************************/
std::string to_string(const std::pair<uint, uint> &pair) {
  return "[" + std::to_string(pair.first) + "," + std::to_string(pair.second) +
         "]";
}

std::string &Trim(std::string &input) {
  input.erase(input.begin(), std::ranges::find_if(input, [](const unsigned char ch) {
    return !std::isspace(ch);
  }));
  input.erase(std::find_if(input.rbegin(), input.rend(), [](const unsigned char ch) {
    return !std::isspace(ch);
  }).base(), input.end());
  return input;
}

/*******************************************************************************
 * @brief Checks if a string represents a valid non-negative integer.
 *
 * @param sv The string view to check.
 * @return bool True if all characters are digits, false otherwise.
 ******************************************************************************/
bool IsNumber(const std::string_view sv) {
  if (sv.empty()) return false;
  return sv.find_first_not_of("0123456789") == std::string::npos;
}

/*******************************************************************************
 * @brief Reads a boolean value (0 or 1) with a default.
 *
 * @param msg The prompt message. Uses a default message if empty.
 * @param default_value The value returned if input is empty.
 * @param in The input stream.
 * @param os The output stream. Prompts are only displayed if os is std::cout.
 * @return bool The parsed boolean (true for non-zero values).
 * @throw std::invalid_argument If input is not a valid number.
 ******************************************************************************/
bool ReadBool(
    const std::string_view msg,
    const bool default_value,
    std::istream &in,
    std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os  << (msg.empty() ? "Enter a bool" : msg)
        << "  (yes or no, default=" << (default_value ? "yes" : "no") << "):";
  }

  std::string input;
  std::getline(in, input);

  if (input.empty()) return default_value;
  Trim(input);
  std::ranges::transform(input, input.begin(), [](const unsigned char ch) {
    return std::tolower(ch);
  });
  if (IsNumber(input)) {
    const uint n = std::stoul(input);
    return n > 0;
  }
  if (input == "yes" || input == "y") {
    return true;
  }
  if (input == "no" || input == "n") {
    return false;
  }
  throw std::runtime_error("ReadBool: Input was not a valid boolean value (expected 0, 1, yes, no, y, or n).");
}

}  // namespace lcs_solver::util