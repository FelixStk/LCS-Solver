#ifndef LCS_SOLVER_UTIL_INPUTPARSER_HPP_
#define LCS_SOLVER_UTIL_INPUTPARSER_HPP_

#include <span>
#include <string>
#include <vector>
namespace lcs_solver::util {

/*******************************************************************************
 * @brief A utility class for parsing command-line arguments
 * @details This class provides functionality to parse command-line arguments,
 * retrieve specific options or check if an option exists
 ******************************************************************************/
class InputParser {
 public:
  InputParser(const int &argc, char **argv);
  [[nodiscard]] const std::string &GetOption(std::string_view option) const;
  [[nodiscard]] const std::string &GetFistOptionOf(std::span<const std::string_view> options) const;
  [[nodiscard]] std::pair<const std::string &, const std::string &> GetPair(std::string_view option) const;
  [[nodiscard]] bool HasOption(std::string_view option) const;
  [[nodiscard]] bool HasAnyOptionOf(std::span<const std::string_view> options) const;
  [[nodiscard]] std::string GetProgramName() const;
  [[nodiscard]] std::string_view GetArgVector(int i) const;
  [[nodiscard]] int GetArgNumber() const;

 private:
  const int & argc_;
  char ** argv_;
  std::vector<std::string> tokens_;
};

}  // namespace lcs_solver::util
#endif  // LCS_SOLVER_UTIL_INPUTPARSER_HPP_