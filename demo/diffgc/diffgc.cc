/*******************************************************************************
 * @file diffgc.cc
 * @author Felix Steinkopp
 * @version 1.0.3
 * @brief Implementation of the diffgc demonstration for the lcs_solver library.
 * @details The core diff logic is implemented in `CalcDiffVector`, which
 * follows a custom difference generation algorithm based on Longest Common
 * Subsequence results.
 ******************************************************************************/

#include "diffgc.h"

#include "FileAdapter.h"

#include <lcs_solver/LibVersion.h>
#include <lcs_solver/algorithms/AlgoFactory.h>
#include <lcs_solver/algorithms/LCS/LCS2_RT.h>
#include <lcs_solver/constraints/BaseConstraint.h>
#include <lcs_solver/constraints/ConstraintFactory.h>
#include <lcs_solver/constraints/ConstraintMap.h>
#include <lcs_solver/problems/ProblemFactory.h>
#include <lcs_solver/util/InOutHelper.h>
#include <lcs_solver/util/InputParser.h>
#include <lcs_solver/util/ParamGenerator.h>
#include <lcs_solver/util/UnicodeRange.h>

#include <chrono>
#include <fstream>

namespace diffgc {

/*******************************************************************************
 * @brief Runs the diff generator
 * @details Creates a problem instance, sets input strings, constraints, and
 *          algorithm, then runs the algorithm and displays the diff.
 *
 * @param[in] argc  Argument count from main()
 * @param[in] argv  Argument vector from main()
 ******************************************************************************/
void Run(const int argc, char* argv[]) {
  auto [file1, file2, parser] = HandleOptions(argc, argv);
  if (parser == nullptr) return;

  const bool gnu = parser->HasOption("--gnu");
  const auto problem = std::make_unique<BaseProblem>();
  const std::unique_ptr<FileAdapter> adapter = SetInputStrings(problem.get(), file1, file2, parser.get());

  HandleConstraint(problem.get(), parser.get(), adapter.get());
  HandleSave(problem.get(), parser.get(), adapter.get());
  SetAlgorithm(problem.get());

  // Run the algorithm and show diff information
  const lcs_solver::util::String lcs = RunAndGetLcs(problem.get());
  const DiffVector diff_vector = CalcDiffVector(lcs, problem->GetStrPtrVector());
  if (gnu) {
    PrintUnifiedDiff(diff_vector, adapter.get());
  }else {
    PrintDiffVector(diff_vector, adapter.get());
  }
}

/*******************************************************************************
 * @brief Sets input strings for a problem from files and creates a FileAdapter
 * @details Reads two input files through a FileAdapter, sets their contents
 *          as the problem's input strings, and returns the adapter instance.
 *
 * @param[in,out] ptr    Problem instance to configure with input strings
 * @param[in]     file1  Path to the first input file (original version)
 * @param[in]     file2  Path to the second input file (modified version)
 * @param[in]     parser Parses terminal arguments as flags or options
 *
 * @return Unique pointer to FileAdapter initialized with the input files
 * @warning Terminates the program with EXIT_FAILURE if file reading fails
 * @note The FileAdapter manages line-to-symbol mapping for diff translation
 ******************************************************************************/
std::unique_ptr<FileAdapter> SetInputStrings(
    BaseProblem *ptr,
    const Path &file1,
    const Path &file2,
    const InputParser *parser) {
  FileAdapter file_adapter{};
  if (parser->HasOption("-j")) {
    const Path config_json = parser->GetOption("-j");
    std::ifstream file_stream(config_json);
    if (!file_stream) {
      throw std::runtime_error("Could not open file " + config_json.string());
    }
    nlohmann::json j = nlohmann::json::parse(file_stream);
    nlohmann::json temp = j["adapter"];
    from_json(temp, file_adapter);
  }
  file_adapter.ChangeFiles(file1, file2, false); // keep mapping
  ptr->SetStrPtrVector(file_adapter.GetStrPtrVector());
  return std::make_unique<FileAdapter>(file_adapter);
}

/*******************************************************************************
 * @brief Configures constraints for the problem instance
 * @details Loads constraints either interactively or from a JSON config file.
 *          When saving, it prompts the user to constraint type and parameters.
 *
 * @param[in,out] ptr       Problem instance to configure
 * @param[in]     parser    Parses terminal arguments as flags or options
 * @param[in]     adapter   Provides a mapping from lines to Symbols
 *
 * @post Problem instance will have zero or one constraint added
 * @note When `save` is true, creates new constraint even if config exists
 ******************************************************************************/
void HandleConstraint(BaseProblem *ptr, const InputParser *parser, const FileAdapter* adapter) {
  using lcs_solver::constraints::BaseConstraint;
  using lcs_solver::constraints::ConstraintFactory;
  using lcs_solver::constraints::ConstraintType;

  // Read in a constraint
  std::shared_ptr<BaseConstraint> constraint;
  if (parser->HasOption("-i")) {
    const std::string prompt = "Would you like to add a constraint?";
    if (lcs_solver::util::ReadBool(prompt, false)) {
      const ConstraintType type = ConstraintFactory::ReadConstraintType();
      const auto& r_map = adapter->GetReverseMap(); // we want to display the actual strings to set a constraint
      constraint = ConstraintFactory::ReadConstraint(ptr->GetStrPtrVector(), type, &r_map);
      ptr->AddConstraint(constraint);
    }
  }else if (parser->HasOption("-j")) {
    const Path config_json = parser->GetOption("-j");
    std::ifstream file_stream(config_json);
    if (!file_stream) {
      throw std::runtime_error("Could not open file " + config_json.string());
    }
    nlohmann::json j = nlohmann::json::parse(file_stream);
    nlohmann::json temp = j.at("constraint");
    lcs_solver::constraints::from_json(j["constraint"], constraint);
    ptr->AddConstraint(constraint);
  }
}

/*******************************************************************************
 * @brief Serializes problem constraints to JSON configuration
 * @details Saves the first constraint from the problem's constraint map to
 *          a JSON file. Prompts for a file path if none is provided.
 *
 * @param[in] ptr     Problem instance containing constraints
 * @param[in] parser  Parses terminal arguments as flags or options
 * @param[in] adapter Provides a mapping from lines to Symbols
 *
 * @note Only saves the first constraint if multiple exists
 * @warning Silently ignores errors during JSON serialization
 ******************************************************************************/
void HandleSave(const BaseProblem *ptr, const InputParser * parser, const FileAdapter* adapter) {
  bool do_saving = parser->HasOption("-s");
  Path file_name;
  if (parser->HasOption("-i")) {
    const std::string prompt = "Would you like to save the constraint?";
    do_saving = lcs_solver::util::ReadBool(prompt, false);
    if (do_saving) {
      file_name = lcs_solver::util::ReadStdString("Enter the save filepath", "./res/constraint.json");
    }
  } else if (parser->HasOption("-s")) {
    file_name = Path(parser->GetOption("-s"));
    if (file_name.empty()) {
      throw std::runtime_error("Error: Invalid save filepath");
    }
  }
  if (do_saving) {
    nlohmann::json j, temp;

    // Add Constraint to JSON
    lcs_solver::constraints::BaseConstraint *constraint_ptr;
    if (const auto &map =ptr->GetConstraints(); map.empty()) {
      constraint_ptr = nullptr;
    } else {
      constraint_ptr = map.begin()->second.get();
    }
    lcs_solver::constraints::to_json(temp, constraint_ptr);
    j["constraint"] = temp;

    // Add Adapter to JSON
    temp.clear();
    to_json(temp, *adapter);
    j["adapter"] = temp;

    // Write to the specified file
    std::ofstream fs;
    fs.open(file_name, std::ofstream::out | std::ofstream::trunc);
    if (!fs) {
      throw std::runtime_error("Could not open file " + file_name.string() + " for writing");
    }
    fs << j.dump(4);
    fs.close();
  }
}

/*******************************************************************************
 * @brief Configures algorithm for the problem instance
 * @details Selects optimal LCS algorithm based on constraints and problem type
 *
 * @param[in,out] ptr Problem instance to configure
 *
 * @post Problem instance will have algorithm set and ready for execution
 * @note Prefers LCS2_RT algorithm with maximum queues and avoids segment trees
 ******************************************************************************/
void SetAlgorithm(BaseProblem *ptr) {
  using lcs_solver::algorithms::AlgoFactory;
  using lcs_solver::algorithms::AlgoType;
  using lcs_solver::constraints::GetLcsAlgoTypes;

  const std::pair<AlgoType, AlgoType> conf_pair = GetLcsAlgoTypes(
      ptr->GetConstraints(),                      // map
      true,                                       // prefer_max_queues
      true,                                       // avoid_segment_tree
      lcs_solver::algorithms::AlgoCategory::LCS2  // effectively uses LCS2_RT
  );
  BaseProblem::AlgorithmPtr algo_ptr = AlgoFactory::Create(
      conf_pair, ptr->GetStrPtrVector(), ptr->GetConstraints());
  ptr->SetAlgorithm(std::move(algo_ptr));
}

/*******************************************************************************
 * @brief Parses command-line arguments and options
 * @details Handles input files and configuration options for diff generation
 *
 * @param[in] argc  Argument count from main()
 * @param[in] argv  Argument vector from main()
 *
 * @return Tuple containing:
 *         - Path to original file
 *         - Path to modified file
 *         - Unique pointer to InputParser instance
 *
 * @note Terminates with usage info if insufficient arguments provided
 ******************************************************************************/
std::tuple<Path, Path, std::unique_ptr<dna::InputParser>> HandleOptions(const int argc, char *argv[]) {
  auto parser = std::make_unique<dna::InputParser>(argc, argv);
  Path file1, file2, config;
  if (parser->HasOption("-h") || parser->HasOption("--help")) {
    PrintUsage(argv[0]);
    return std::make_tuple("","",nullptr);
  }
  if (parser->HasOption("--version")) {
    PrintVersion(argv[0]);
    return std::make_tuple("","",nullptr);
  }
  if (argc >= 3) {
    file1 = argv[1];
    file2 = argv[2];
  } else {
    PrintUsage(argv[0]);
    return std::make_tuple("","",nullptr);
  }
  return std::make_tuple(file1, file2, std::move(parser));
}

std::string GetBinName(char name[]) {
  const std::filesystem::path path(name);
  return path.filename().string();
}

void PrintVersion(char name[]) {
  std::cout << GetBinName(name);
  std::cout << " Version 1.0.2 - (based on " << lcs_solver::Name() << ")";
  std::cout << std::endl;
}

void PrintUsage(char name[]) {
  std::cout << "Usage: " << GetBinName(name);
  std::cout << R"( file1 file2 [-j config.json] [-s] [--gnu]
Parameters:
  file1 file2    : Filepaths to the original and modified file
Options:
  -j config.json : A JSON file that contains zero or one constraint
  -s file.json   : Save a constraint to a JSON file
  -i             : Prompt for a constraint type and parameters
  --gnu          : Print in the unified output format
  --help or -h   : Display this usage message
  --version      : Display the version of the program
)" << std::endl;
}

/*******************************************************************************
 * @brief Executes LCS algorithm and extracts the lcs from an embedding
 * @details Runs the configured algorithm and converts the solution to a diff
 *          compatible format using RangeTreeSolution
 *
 * @param[in] ptr Problem instance with configured algorithm
 * @return Embedded LCS string from the first valid solution
 *
 * @throws std::runtime_error If the generation of a RangeTreeSolution failed
 * @warning Assumes single valid embedding exists in a solution
 ******************************************************************************/
lcs_solver::util::String RunAndGetLcs(const BaseProblem *ptr) {
  using lcs_solver::algorithms::lcs::LCS2_RT;
  if (ptr == nullptr) return {};
  auto lcs = lcs_solver::util::String();
  const std::unique_ptr<BaseProblem::BaseSolution> result = ptr->ExecuteAlgo();
  const auto sol = dynamic_cast<LCS2_RT::RangeTreeSolution *>(result.get());
  if (!sol) {
    throw std::runtime_error("Error: The algorithm did not return a RangeTreeSolution");
  }
  if (auto it = sol->begin(); it != sol->end()) {
    if (const auto *points_ptr = it.operator->(); points_ptr != nullptr) {
      const auto embedding = points_ptr->GenEmbedding(0, true);
      lcs = embedding.getEmbeddedStr();
    }
  }
  return lcs;
}

/*******************************************************************************
 * @brief Computes difference vector between original and modified versions
 * @details Implements custom diff algorithm using LCS results to identify:
 *          - Kept elements (common sequence)
 *          - Added elements (new in modified)
 *          - Removed elements (deleted from original)
 *
 * @param[in] lcs  Longest Common Subsequence result
 * @param[in] spv  String pointers to original and modified versions
 *
 * @return Vector of (symbol, change type) pairs representing full diff
 *
 * @throws std::runtime_error If input doesn't contain exactly two strings
 * @see GetDescription() for algorithm pseudocode
 ******************************************************************************/
DiffVector CalcDiffVector(
    const lcs_solver::util::String &lcs,
    const BaseProblem::StrPtrVec &spv) {
  if (spv.size() != 2) {
    throw std::runtime_error("Error: The number of strings must be 2");
  }
  // Initialize Variables and iterators
  std::vector<std::pair<FileAdapter::Symbol, ChangeType>> v;
  auto it1 = spv[0]->begin();
  auto it2 = spv[1]->begin();
  auto lcs_iter = lcs.begin();
  const auto end1 = spv[0]->end();
  const auto end2 = spv[1]->end();
  const auto lcs_end = lcs.end();

  while (it1 != end1 || it2 != end2) {
    if (lcs_iter == lcs_end) {  // dereferencing lcsIter could be problematic
      while (it1 != end1) {
        v.emplace_back(*it1, ChangeType::kRemove);
        ++it1;
      }
      while (it2 != end2) {
        v.emplace_back(*it2, ChangeType::kAdd);
        ++it2;
      }
      continue;
    }
    if (it1 == end1) {
      while (it2 != end2) {
        v.emplace_back(*it2, ChangeType::kAdd);
        ++it2;
      }
      continue;
    }
    if (it2 == end2) {
      while (it1 != end1) {
        v.emplace_back(*it1, ChangeType::kRemove);
        ++it1;
      }
      continue;
    }
    if (it1 != end1 && it2 != end2 && *it1 == *it2) {
      v.emplace_back(*it1, ChangeType::kKeep);
      ++it1;
      ++it2;
      ++lcs_iter;
      continue;
    }
    auto x = it1;
    while (x != end1 && *x != *lcs_iter) {
      ++x;
    }
    auto y = it2;
    while (y != end2 && *y != *lcs_iter) {
      ++y;
    }
    while (it1 != x) {
      v.emplace_back(*it1, ChangeType::kRemove);
      ++it1;
    }
    while (it2 != y) {
      v.emplace_back(*it2, ChangeType::kAdd);
      ++it2;
    }
  }
  return v;
}

/*******************************************************************************
 * @brief Displays diff vector in human-readable format
 * @details Prints symbols with change markers (+/-/=). Can translate symbols
 *          back to original file lines using FileAdapter mapping.
 *
 * @param[in] dv              Diff vector to display
 * @param[in] adapter         FileAdapter for symbol-to-line translation
 * @param[in] translate_back  Enable line number translation
 *
 * @throws std::runtime_error If translation requested without a valid adapter
 * @note Without translation, outputs raw symbol values
 ******************************************************************************/
void PrintDiffVector(const DiffVector &dv, const FileAdapter *adapter,
                     const bool translate_back) {
  if (translate_back && adapter == nullptr) {
    throw std::runtime_error("Error: adapter required for translation");
  }
  size_t i = 0;
  for (const auto &[symbl, change] : dv) {
    std::cout << static_cast<char>(change) << " | ";
    if (translate_back) {
      std::cout << adapter->GetLine(symbl);
    } else {
      std::cout << lcs_solver::util::to_string(symbl, true);
    }
    ++i;
    if (i != dv.size()) {
      std::cout << "\n";
    }
  }
}

/*******************************************************************************
 * @brief Prints a unified diff view to standard output.
 * @details Outputs differences between two versions in the standard unified diff
 *          format used by tools like `diff` and `git diff`, including optional
 *          context lines before and after each change hunk.
 *
 * @param[in] dv              The diff vector (from CalcDiffVector).
 * @param[in] adapter         Optional FileAdapter for symbol-to-line translation.
 * @param[in] translate_back  Whether to translate symbols using the adapter.
 *
 * @throws std::runtime_error if `translate_back` is true but `adapter` is null.
 *
 * @note The output format includes:
 *       - `@@ -start,len +start,len @@` headers for each hunk
 *       - Lines starting with:
 *           - `' '` for unchanged lines (context)
 *           - `'-'` for deletions
 *           - `'+'` for additions
 * @see https://www.gnu.org/software/diffutils/manual/html_node/Unified-Format.html
 ******************************************************************************/
void PrintUnifiedDiff(
    const DiffVector &dv,
    const FileAdapter *adapter,
    const bool translate_back) {
  if (translate_back && adapter == nullptr) {
    throw std::runtime_error("Error: adapter required for translation");
  }

  // File header: Get paths and last modification times
  if (adapter != nullptr) {
    const std::span<const Path> files = adapter->GetFiles();
    if (files.size() == 2) {
      auto print_timestamp = [](const std::filesystem::path &path) {
        using namespace std::chrono;
        const auto ftime = std::filesystem::last_write_time(path);
        const auto sys_time = clock_cast<system_clock>(ftime);
        zoned_time zt{current_zone(), sys_time};
        std::cout << format("{:%F %T %z}", zt);
      };

      std::cout << "--- " << files[0].string() << "\t";
      print_timestamp(files[0]);
      std::cout << "\n+++ " << files[1].string() << "\t";
      print_timestamp(files[1]);
      std::cout << "\n";
    }
  }

  // Hunk generation
  size_t orig_line = 1, mod_line = 1;
  size_t hunk_orig_start = 0, hunk_mod_start = 0;
  size_t hunk_orig_count = 0, hunk_mod_count = 0;
  std::vector<std::string> hunk_lines;

  auto flush_hunk = [&] {
    if (hunk_lines.empty()) return;
    std::cout << "@@ -" << hunk_orig_start << (hunk_orig_count != 1 ? "," + std::to_string(hunk_orig_count) : "")
              << " +" << hunk_mod_start << (hunk_mod_count != 1 ? "," + std::to_string(hunk_mod_count) : "")
              << " @@\n";
    for (const auto &line : hunk_lines) {
      std::cout << line << '\n';
    }
    hunk_lines.clear();
    hunk_orig_count = 0;
    hunk_mod_count = 0;
  };

  for (const auto &[sym, change] : dv) {
    std::string line = translate_back ? adapter->GetLine(sym)
                                      : lcs_solver::util::to_string(sym, true);
    switch (change) {
      case ChangeType::kKeep:
        flush_hunk();
        ++orig_line;
        ++mod_line;
        break;
      case ChangeType::kRemove:
        if (hunk_lines.empty()) {
          hunk_orig_start = orig_line;
          hunk_mod_start = mod_line;
        }
        hunk_lines.push_back("-" + line);
        ++hunk_orig_count;
        ++orig_line;
        break;
      case ChangeType::kAdd:
        if (hunk_lines.empty()) {
          hunk_orig_start = orig_line;
          hunk_mod_start = mod_line;
        }
        hunk_lines.push_back("+" + line);
        ++hunk_mod_count;
        ++mod_line;
        break;
    }
  }
  flush_hunk();
}


std::string_view GetDescription() {
  static constexpr std::string_view msg = R"DESC(Pseudocode:
> Diff(s1, s2, lcs):
>   let v be an empty array
>   let it1 = s1.begin(), it2 = s2.begin(), lcsIter = lcs.begin()
>   let end1 = s1.end(), end2 = s2.end(), lcsEnd = lcs.end()
>
>   while (it1 != s1.end() || it2 != s2.end())
>     if lcsIter == lcsEnd then
>       mark rest of s1 as removed and rest of s2 as added and return
>     if it1 == end1 then mark rest of s2 as added and return
>     if it2 == end2 then mark rest of s1 as added and return
>     if *it1 == *it2 increment it1, it2 and lcsIter and continue
>     let x = it1 and increment it until *x == *lcsIter
>     let y = it2 and increment it until *y == *lcsIter
>     while it1 != x
>       mark *it1 as removed and increment it1
>     while it2 != y
>       mark *it2 as removed and increment it2
)DESC";
  return msg;
}

/*******************************************************************************
 * @brief Generates Latin lowercase alphabet symbols
 * @return Vector containing Unicode symbols from U+0061 to U+007A ('a'-'z')
 ******************************************************************************/
std::vector<lcs_solver::util::Symbol> GetAlphabet() {
  using lcs_solver::util::UnicodeRange;
  const UnicodeRange range(UnicodeRange::Alphabet::kLatinLower);
  return range.symbAlphabetVector();
}

/*******************************************************************************
 * @brief Creates randomized LCS problem example with unique solution
 * @return AlgoParam structure containing:
 *         - Name: Generated example name
 *         - Strings: Original/modified string pair
 *         - Constraints: Monotonicity constraint (MC)
 *         - Solution: Reference solution (not for direct diff use)
 *
 * @note Uses fixed seed (19,035,066) for reproducible examples
 ******************************************************************************/
lcs_solver::util::AlgoParam GenerateRndExample() {
  auto vector = lcs_solver::util::ParamGenerator::genWithUniqSol(
      lcs_solver::constraints::ConstraintType::MC,
      19035066,      // std::random_device{}()  // seed
      1,             // Number of Problems to generate
      {              // v[i] = {min(s[i].size(), min(s[i].size()}
       {5, 10},      // bound
       {5, 10}},     //
      {2, 5},        // Pair (l,u) such that l <= llcs(s[0],s[1]) <= u
      GetAlphabet()  // v[0] is the only common symbol of s[0], s[1]
  );
  return vector.front();
}

}  // namespace diffgc