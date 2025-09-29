#ifndef LCS_SOLVER_INCLUDE_UTIL_FILE_ADAPTER_H_
#define LCS_SOLVER_INCLUDE_UTIL_FILE_ADAPTER_H_

#include "lcs_solver/util/CommonTypes.h"
#include "lcs_solver/util/StrPtrVector.h"

#include <filesystem>
#include <unordered_map>
#include <vector>

#include <nlohmann/json.hpp>

namespace diffgc {

/**
 * @class FileAdapter
 * @brief Translates between file contents and symbolic representation
 * @details Maintains bidirectional mapping between file lines and unique
 * symbols. Processes input files to create normalized String representations
 * used by LCS algorithms. Enables translation between diff results and original
 * file lines.
 *
 * @note Manages state to avoid redundant file processing. Symbols are assigned
 * sequentially starting from kFirstSymbol (0). Each unique line gets a unique
 * symbol.
 */
class FileAdapter {
 public:
  using Path = std::filesystem::path;
  using String = lcs_solver::util::String;
  using Symbol = lcs_solver::util::Symbol;
  using StrPtrVector = lcs_solver::util::StrPtrVector;

  FileAdapter(const std::vector<Path> files, std::unordered_map<Symbol, std::string> map);

  /**
   * @brief Construct from multiple files
   * @param files Vector of file paths to process
   * @post Internal state marked as changed, files will be processed on first
   * access
   */
  explicit FileAdapter(const std::vector<Path> &files);

  /**
   * @brief Construct from two files for diff comparison
   * @param file1 Path to the original version file
   * @param file2 Path to the modified version file
   */
  FileAdapter(const Path &file1, const Path &file2);

  /// Default constructor creates empty adapter
  FileAdapter() = default;

  /**
   * @brief Saves adapter instance to a file
   * @param path Path to the save file
   * @param adapter FileAdapter instance to be saved
   */
  static void ToJson(const Path& path, const FileAdapter* adapter);

  /**
   * @brief SavesI adapter instance to a file
   * @param path Path to the save file
   */
  void ToJson(const Path& path) const;

  /**
   * @brief Get current file paths
   * @return Read-only view of managed file paths
   */
  [[nodiscard]] std::span<const Path> GetFiles() const;

  /**
   * @brief Get processed strings for LCS computation
   * @return Vector of shared pointers to processed symbol strings
   * @throws std::runtime_error If any file cannot be opened
   * @note Processes files if needed (when changed_ is true)
   */
  StrPtrVector GetStrPtrVector();

  /**
   * @brief Get a symbol for a file line (const version)
   * @param line Exact line content to look up
   * @return Associated symbol or 0 if not found
   */
  Symbol GetSymbol(const std::string &line) const;

  /**
   * @brief Get original line content for a symbol
   * @param symbol Numeric symbol to translate
   * @return Original line content or empty string if not found
   */
  std::string GetLine(Symbol symbol) const;

  /**
   * @brief Reset adapter to empty state
   * @post Clears all file references, mappings, and processed data
   */
  void Clear();

  /**
   * @brief Update files to compare (pair version)
   * @param file1 New original file path
   * @param file2 New modified file path
   * @param clear_maps Whether to reset symbol mappings
   * @post Marks state as changed, files will be reprocessed
   */
  void ChangeFiles(
    const Path &file1,
    const Path &file2,
    bool clear_maps = true);

  /**
   * @brief Update files to compare (`std::vector` version)
   * @param files New set of file paths to process
   * @param clear_maps Whether to reset existing symbol mappings
   * @post Marks state as changed, files will be reprocessed
   */
  void ChangeFiles(const std::vector<Path> &files, bool clear_maps = true);

  /**
   * @brief Getter for line_to_symbol_map_
   * @return Map to translate symbols into std::strings
   */
  const std::unordered_map<std::string, Symbol>& GetMap() const;

  /**
   * @brief Getter for line_to_symbol_map_
   * @return Map to translate symbols into std::strings
   */
  const std::unordered_map<Symbol, std::string>& GetReverseMap() const;

 private:
  /// @brief Process files and update symbol mappings
  /// @throws std::runtime_error If any file cannot be opened
  void Process();

  /**
   * @brief Get or create symbol for a line
   * @param line Line content to register
   * @return Existing or newly created symbol for this line
   * @note Maintains bijection between lines and symbols
   */
  Symbol GetOrCreateSymbol(const std::string &line);

  /**
   * @brief Process single file into symbol string
   * @param filename Path to file to process
   * @return String of symbols representing file lines
   * @throws std::runtime_error If the file cannot be opened
   */
  String ProcessFile(const Path &filename);

  std::vector<std::filesystem::path> files_;
  static constexpr Symbol k_first_symbol = 0x61;
  Symbol next_symbol_ = k_first_symbol;
  std::unordered_map<std::string, Symbol> line_to_symbol_map_;
  std::unordered_map<Symbol, std::string> symbol_to_line_map_;
  StrPtrVector spv_;
  bool changed_ = false;
};

/**
 * @brief Serialize FileAdapter to JSON
 * @param[out] j Output JSON object
 * @param[in] p FileAdapter to serialize
 * @note Only serializes file paths, not symbol mappings
 */
void to_json(nlohmann::json &j, const FileAdapter &p);

/**
 * @brief Deserialize FileAdapter from JSON
 * @param[in] j Input JSON object
 * @param[out] p FileAdapter to populate
 * @post Resets symbol mappings and marks state as changed
 */
void from_json(const nlohmann::json &j, FileAdapter &p);

}  // namespace diffgc
#endif  // LCS_SOLVER_INCLUDE_UTIL_FILE_ADAPTER_H_
