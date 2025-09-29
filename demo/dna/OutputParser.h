#ifndef DEMO_DNA_OUTPUTPARSER_H_
#define DEMO_DNA_OUTPUTPARSER_H_

#include <filesystem>
#include <iostream>
#include <map>
#include <string>

namespace dna {

/**
 * @class OutputParser
 * @brief Formats benchmark result according to TSV format
 * @see https://afproject.org/app/benchmark/genetree/swisstree/dataset/
 */
class OutputParser {
 public:
  using Path = std::filesystem::path;
  using PathPair = std::pair<Path, Path>;

  OutputParser() = default;
  void Add(const PathPair& pair, size_t llcs);

  void Print(std::ostream& os = std::cout);
  static void Print(PathPair, size_t llcs, std::ostream& os = std::cout);
  void Save(const std::filesystem::path& file_path);

 private:
  std::map<std::string, std::map<std::string, size_t>> map_;  ///< Stores LLCS
};

}  // namespace dna

#endif  // DEMO_DNA_OUTPUTPARSER_H_
