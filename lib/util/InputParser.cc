/*******************************************************************************
 * @file InputParser.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief InputParser is a utility class for parsing command-line arguments
 ******************************************************************************/

#include "util/InputParser.h"

#include <algorithm>
#include <filesystem>

namespace lcs_solver::util {

/*******************************************************************************
 * @brief Constructs an InputParser object and parses the command-line arguments
 * @param argc The number of command-line arguments
 * @param argv An array of command-line argument strings
 * @details This constructor initializes the InputParser object by assigning the
 * command-line arguments to the tokens_ vector, skipping the program name
 ******************************************************************************/
InputParser::InputParser(const int &argc, char **argv) :
  argc_(argc), argv_(argv)
{
  tokens_.assign(argv + 1, argv + argc); // skip program name
}

/*******************************************************************************
 * @brief Retrieves the value associated with a specific command-line option.
 * @param option The command-line option to retrieve.
 * @return The value associated with the specified option, or an empty string
 * if the option is not found.
 * @details This function uses the tokens_ vector to find the value associated
 * with the given option. It returns a reference to the value if found, or an
 * empty string if not.
 ******************************************************************************/
const std::string& InputParser::GetOption(const std::string_view option) const {
  if (auto itr = std::ranges::find(tokens_, option);
      itr != tokens_.end() && ++itr != tokens_.end()) {
    return *itr;
  }
  static const std::string empty_string;
  return empty_string;
}

const std::string& InputParser::GetFistOptionOf(
    std::span<const std::string_view> options) const {
  for (const auto &option : options) {
    if (const std::string& ref = GetOption(option); !ref.empty()) {
      return ref;
    }
  }
  static const std::string empty_string;
  return empty_string;
}

/*******************************************************************************
 * Retrieves the pair of values associated with a specific command-line option
 * @param option The command-line option to retrieve
 * @return A pair containing the values associated with the specified option, or
 * a pair of empty strings if the option is not found
 * @details This function uses the tokens_ vector to find the pair of values
 * associated with the given option. It returns a pair of references to the
 * values if found, or a pair of empty strings if not
 ******************************************************************************/
std::pair<const std::string &, const std::string &> InputParser::GetPair(
    const std::string_view option) const {
  static const std::string empty_string;
  auto itr = std::ranges::find(tokens_, option);
  if (itr != tokens_.end()) {
    ++itr;
    if (itr != tokens_.end()) {
      const std::string &first = *itr;
      ++itr;
      if (itr != tokens_.end()) {
        const std::string &second = *itr;
        return {first, second};
      }
    }
  }
  return {empty_string, empty_string};
}

/*******************************************************************************
 * @brief Checks if a specific command-line option exists
 * @param option The command-line option to check
 * @return `true` if the option exists, `false` otherwise
 ******************************************************************************/
bool InputParser::HasOption(const std::string_view option) const {
  return std::ranges::find(tokens_, option) != tokens_.end();
}

/*******************************************************************************
 * @brief Checks if a specific command-line option exists
 * @param options The command-line options to check
 * @return `true` if one of the options exists, `false` otherwise
 ******************************************************************************/
bool InputParser::HasAnyOptionOf(
    std::span<const std::string_view> options) const {
  return std::ranges::any_of(
      options, [this](const std::string_view opt) { return HasOption(opt); });
}

/*******************************************************************************
 * @brief Getter for the program name
 * @return Name of the program
 ******************************************************************************/
std::string InputParser::GetProgramName() const {
  if (argc_ < 1)
    return "";
  const std::filesystem::path path(argv_[0]); // trim path to the executable
  return path.filename().string();
}

/*******************************************************************************
 * @brief Getter for argv
 * @param i Integer to specify the ith command line argument
 * @return A `std::string_view` of the ith command line argument with that the
 *         instance was constructed
 ******************************************************************************/
std::string_view InputParser::GetArgVector(int i) const {
  if (i < 0 || i >= argc_) return {};
  return argv_[i];
}

/*******************************************************************************
 * @brief Getter for argc
 * @return The number of command line arguments
 ******************************************************************************/
int InputParser::GetArgNumber() const {
  return argc_;
}

}  // namespace lcs_solver::util