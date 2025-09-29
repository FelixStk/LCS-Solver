/*******************************************************************************
 * @file FileAdapter.cc
 * @author Felix Steinkopp
 * @version 1.0.1
 * @brief Implements file-to-symbol translation for LCS-based diff generation
 * @details Provides bidirectional mapping between file lines and numeric
 *          symbols, enabling efficient LCS computation while maintaining
 *          traceability to original content. Processes files on-demand and
 *          maintains the processing state.
 ******************************************************************************/
#include "FileAdapter.h"

#include <fstream>
#include <stdexcept>

namespace diffgc {

/**
 * @brief Construct with multiple files and a symbol mapping
 * @details Moves files, moves `map` into `symbol_to_line_map_` (accessible via
 *          `GetReverseMap()`) and fills the appropriate normal map. Also sets
 *          `next_symbol` based on `k_fist_symbol` and the size of `map`
 * @param files Files to process (typically original and modified versions)
 * @param map Map from Symbols to standard strings
 * @post Marks adapter as needing processing (changed_ = true)
 */
FileAdapter::FileAdapter(std::vector<Path> files,
                         std::unordered_map<Symbol, std::string> map)
                           : files_(std::move(files)), symbol_to_line_map_(std::move(map))
{
  for (const auto& [key, value] : symbol_to_line_map_) {
    line_to_symbol_map_.emplace(value, key);
  }
  next_symbol_ = k_first_symbol + symbol_to_line_map_.size();
  changed_ = true;
}


/**
 * @brief Construct with multiple files
 * @param files Files to process (typically original and modified versions)
 * @post Marks adapter as needing processing (changed_ = true)
 */
FileAdapter::FileAdapter(const std::vector<Path> &files) : files_(files) {
  changed_ = true;
}

/**
 * @brief Construct with two files for comparison
 * @param file1 Original version file path
 * @param file2 Modified version file path
 */
FileAdapter::FileAdapter(const Path& file1, const Path& file2)
    : FileAdapter(std::vector{file1, file2}) {
  changed_ = true;
}

void FileAdapter::ToJson(const Path& path, const FileAdapter* adapter) {
  std::ofstream fs;
  fs.open(path, std::ofstream::out | std::ofstream::trunc);
  if (!fs) {
    throw std::runtime_error("Could not open file " + path.string() + " for writing");
  }
  nlohmann::json j;
  if (adapter == nullptr) {
    FileAdapter temp{};
    to_json(j, temp);
  } else {
    to_json(j, *adapter);
  }
  fs << j.dump(4);
  fs.close();
}

void FileAdapter::ToJson(const Path& path) const {
  ToJson(path, this);
}

/**
 * @brief Get current file paths
 * @return Read-only view of managed file paths
 */
std::span<const FileAdapter::Path> FileAdapter::GetFiles() const {
  return files_;
}

/**
 * @brief Process files if needed and return symbol strings
 * @return Vector of shared pointers to processed symbol strings
 * @throws std::runtime_error If any file cannot be opened
 * @details Processes files only if changes occurred since last access
 */
FileAdapter::StrPtrVector FileAdapter::GetStrPtrVector() {
  if (changed_) {
    Process();
  }
  return spv_;
}

/**
 * @brief Retrieve symbol for a line (const-safe)
 * @param line Exact line content to look up
 * @return Associated symbol or 0 (kFirstSymbol) if not found
 */
FileAdapter::Symbol FileAdapter::GetSymbol(const std::string& line) const {
  if (const auto it = line_to_symbol_map_.find(line);
      it != line_to_symbol_map_.end()) {
    return it->second;
  }
  return 0;
}

/**
 * @brief Get original line content for a symbol
 * @param symbol Numeric symbol to translate
 * @return Original line content or empty string if not found
 */
std::string FileAdapter::GetLine(const Symbol symbol) const {
  if (const auto it = symbol_to_line_map_.find(symbol);
      it != symbol_to_line_map_.end()) {
    return it->second;
  }
  return "";
}

/**
 * @brief Reset adapter to initial empty state
 * @post Clears all files, mappings, and resets symbol counter
 */
void FileAdapter::Clear() {
  files_.clear();
  line_to_symbol_map_.clear();
  symbol_to_line_map_.clear();
  next_symbol_ = k_first_symbol;
  spv_.clear();
  changed_ = false;
}

/**
 * @brief Internal file processing method
 * @throws std::runtime_error If any file cannot be opened
 * @details Converts all managed files to symbol strings using ProcessFile
 */
void FileAdapter::Process() {
  if (changed_) {
    spv_.clear();
    for (const auto &file : files_) {
      String s = ProcessFile(file.string());
      spv_.emplace_back(std::make_shared<String>(s));
    }
    changed_ = false;
  }
}

/**
 * @brief Update files to compare (pair version)
 * @param file1 New original file path
 * @param file2 New modified file path
 * @param clear_maps Whether to reset symbol mappings
 * @post Marks state as changed, files will be reprocessed
 */
void FileAdapter::ChangeFiles(
    const Path &file1,
    const Path &file2,
    const bool clear_maps) {
  if (clear_maps)
    Clear();
  files_ = {file1, file2};
  changed_ = true;
}

/**
 * @brief Update files to compare (`std::vector` version)
 * @param files New set of file paths to process
 * @param clear_maps Whether to reset existing symbol mappings
 * @post Marks state as changed, files will be reprocessed
 */
void FileAdapter::ChangeFiles(
    const std::vector<Path> &files,
    const bool clear_maps) {
  if (clear_maps)
    Clear();
  files_ = files;
  changed_ = true;
}

/**
 * @brief Get or create symbol for a line
 * @param line Line content to register
 * @return Existing or newly created symbol for this line
 * @note Maintains bijection between lines and symbols
 */
FileAdapter::Symbol FileAdapter::GetOrCreateSymbol(const std::string& line) {
  if (const auto it = line_to_symbol_map_.find(line);
      it != line_to_symbol_map_.end()) {
    return it->second;
  }
  const Symbol new_symbol = next_symbol_++;
  line_to_symbol_map_[line] = new_symbol;
  symbol_to_line_map_[new_symbol] = line;
  return new_symbol;
}

const std::unordered_map<FileAdapter::Symbol, std::string>&
FileAdapter::GetReverseMap() const {
  return symbol_to_line_map_;
}

const std::unordered_map<std::string, FileAdapter::Symbol>&
FileAdapter::GetMap() const {
  return line_to_symbol_map_;
}

/**
 * @brief Process a single file into a symbol sequence
 * @param filename File to process
 * @return Sequence of symbols representing file lines
 * @throws std::runtime_error For file open errors
 * @note Empty lines are preserved with unique symbols
 */
FileAdapter::String FileAdapter::ProcessFile(const Path &filename) {
  std::ifstream file(filename.string());
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + filename.string());
  }

  String result;
  std::string line;
  while (std::getline(file, line)) {
    // Remove any trailing whitespace or newline characters
    // line = line.substr(0, line.find_last_not_of(" \n\r\t") + 1);
    // Skip empty lines
    // if (!line.empty()) {
    //   const Symbol symbol = GetOrCreateSymbol(line);
    //   result.push_back(symbol);
    // }
    const Symbol symbol = GetOrCreateSymbol(line);
    result.push_back(symbol);
  }
  return result;
}

/**
 * @brief Serialize FileAdapter to JSON
 * @param[out] j Output JSON object
 * @param[in] p FileAdapter to serialize
 * @note Only serializes file paths, not symbol mappings
 */
void to_json(nlohmann::json& j, const FileAdapter& p) {
  j["Files"] = p.GetFiles();
  j["MapR"] = p.GetReverseMap();
}

/**
 * @brief Deserialize FileAdapter from JSON
 * @param[in] j Input JSON object
 * @param[out] p FileAdapter to populate
 * @post Resets symbol mappings and marks state as changed
 */
void from_json(const nlohmann::json& j, FileAdapter& p) {
  std::vector<std::filesystem::path> files;
  std::unordered_map<FileAdapter::Symbol, std::string> map_r;
  j.at("Files").get_to(files);
  j.at("MapR").get_to(map_r);
  p = FileAdapter(files, map_r);
}

}  // namespace diffgc