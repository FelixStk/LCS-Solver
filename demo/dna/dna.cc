/*******************************************************************************
 * @file dna.cc
 * @author Felix Steinkopp
 * @version 1.1.0
 * @brief Implementation of the dna demonstration for the lcs_solver library
 * @details Run (constraint) llcs algorithms on the swisstree dataset
 * @see https://afproject.org/app/benchmark/genetree/swisstree/dataset/
 ******************************************************************************/

#include "dna.h"

#include <lcs_solver/LibVersion.h>
#include <lcs_solver/algorithms/AlgoFactory.h>
#include <lcs_solver/algorithms/solutions/UnsignedSolution.h>
#include <lcs_solver/problems/ProblemFactory.h>
#include <lcs_solver/util/InputParser.h>
#include <nlohmann/json.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <regex>

namespace dna {

void Run(const int argc, char* argv[]) {
  // Setup data, problem params and other configure options
  const std::unique_ptr<InputParser> parser = HandleArgs(argc, argv);
  if (parser == nullptr)
    return;  // HandleArgs has already printed an appropriate msg
  const bool verbose = GetVerbose(parser.get(), false);
  const auto [save, save_file] = GetSaveInfo(parser.get());
  const ProblemType prob_t = GetProblemType(parser.get());
  const AlgoType algo_t = GetAlgorithmType(parser.get());
  const std::unique_ptr<GapVector> gc = GetGapVector(parser.get());
  auto [alph, l_gaps, r_gaps] = GetSigmaVectors(parser.get());
  const std::vector<PathPair> pairs = GetCombinations(parser.get());

  // Run LLCS Algorithm on all file combinations in pairs
  const std::unique_ptr<OutputParser> result = DoMainLoop(pairs,
                                                          prob_t,
                                                          algo_t,
                                                          gc.get(),
                                                          alph.get(),
                                                          l_gaps.get(),
                                                          r_gaps.get(),
                                                          verbose);

  // Print or save results
  if (!save) {
    result->Print();
  } else {
    result->Save(save_file);
  }
}

std::unique_ptr<OutputParser> DoMainLoop(const std::vector<PathPair>& comb_vec,
                                         const ProblemType& prob_t,
                                         const AlgoType& algo_t,
                                         const GapVector* gc,
                                         const SymbolVector* alph,
                                         const GapVector* l,
                                         const GapVector* r,
                                         const bool verbose) {
  using lcs_solver::algorithms::solutions::UnsignedSolution;
  using lcs_solver::problems::ProblemFactory;

  auto result = std::make_unique<OutputParser>();
  const size_t total = comb_vec.size();
  size_t counter = verbose ? ShowProgress(0, total) : 0;
  for (const auto& file_pair : comb_vec) {
    // Set up a problem
    std::array s = {From(file_pair.first), From(file_pair.second)};
    const auto [gc_span, alph_span, l_span, r_span] = GetCorrectedParams(s, gc, alph, l, r);
    const auto ptr = ProblemFactory::Create(prob_t, s, gc_span, alph_span, l_span, r_span, algo_t);

    // Execute Algorithm and save result
    const auto sol = ptr->ExecuteAlgo<UnsignedSolution>();
    result->Add(file_pair, s[0].length() + s[1].length() - 2*sol->GetNumber());
    //result->Add(file_pair, 0);

    // Update Progressbar
    counter = verbose ? ShowProgress(counter, total) : counter + 1;
  }
  return result;
}

/**
 * @brief Backs up the constraint vectors and generates appropriate spans
 * @param s Span of input strings
 * @param gc Predefined gap constraints (must
 * @param alph Alphabet used for the strings
 * @param l Predefined values of the left map values with respect to alph
 * @param r Predefined values of the right map values with respect to alph
 * @return tuple containing corrected span versions (of gc, alph, l and r)
 * @note The current implementation assumes that gc, alph, l and r do not change
 *       from problem to problem. It only corrects missing Symbols in alph, l
 *       and r by adding Symbols in s or tuples (0,max(s[0].size(),s[1].size()).
 *       If the gc tuple is not long enough, it will be extended with relaxed
 *       gap length constraint data (0,max(s[0].size(),s[1].size())
 */
std::tuple<GapSpan, SymbolSpan, GapSpan, GapSpan> GetCorrectedParams(
    const std::span<const String> s,
    const GapVector* gc,
    const SymbolVector* alph,
    const GapVector* l,
    const GapVector* r) {
  static GapVector gc_s;
  static SymbolVector alph_s;
  static GapVector l_s;
  static GapVector r_s;

  // Setup gc_ptr
  const uint min = std::min(s[0].length(), s[1].length());
  const uint max = std::max(s[0].length(), s[1].length());
  if (gc){
    gc_s.clear();
    for (const auto& gap : *gc) {
      gc_s.emplace_back(gap);
    }
    while (gc_s.size() < min-1) {
      gc_s.emplace_back(0, max);
    }
  }

  // Convert alph, l and r into accessible information
  using SigmaTupleMap = std::map<Symbol, std::pair<uint, uint>>;
  static SigmaTupleMap there_l;
  static SigmaTupleMap there_r;
  if (alph) {
    if (there_l.empty()) {
      for (uint i = 0; i < alph->size(); i++) {
        there_l.emplace(alph->at(i), l->at(i));
      }
    }
    if (there_r.empty()) {
      for (uint i = 0; i < alph->size(); i++) {
        there_r.emplace(alph->at(i), r->at(i));
      }
    }
  }

  // Setup alph_ptr
  const SymbolVector temp_alph = lcs_solver::util::CalcAlphabet(s);
  const bool alph_update = temp_alph.size() > alph_s.size();
  if (alph_update) {
    alph_s.clear();
    for (const auto& symb : temp_alph) {
      alph_s.emplace_back(symb);
    }
  }

  // Setup l_ptr
  if (alph_update && !there_l.empty()) {
    for (const auto& symb : temp_alph) {
      if (auto it = there_l.find(symb); it != there_l.end()) {
        l_s.emplace_back(it->second);
      } else {
        l_s.emplace_back(0, max);
      }
    }
  }

  // Setup r_ptr
  if (alph_update && !there_r.empty()) {
    for (const auto& symb : temp_alph) {
      if (auto it = there_r.find(symb); it != there_r.end()) {
        r_s.emplace_back(it->second);
      } else {
        r_s.emplace_back(0, max);
      }
    }
  }
  return {std::span(gc_s).subspan(0, min-1), alph_s, l_s, r_s};
}

std::unique_ptr<InputParser> HandleArgs(const int argc, char* argv[]) {
  auto ptr = std::make_unique<InputParser>(argc, argv);
  Path file = "";
  if (ptr->HasAnyOptionOf(std::span(kOptVersion))) {
    PrintVersion(argv[0]);
    return nullptr;
  }
  if (ptr->HasAnyOptionOf(std::span(kOptHelp))) {
    PrintUsage(argv[0], ptr.get());
    return nullptr;
  }
  // if (!ptr->HasAnyOptionOf(std::span(kOptAlgoType))) {
  //   PrintUsage(argv[0], ptr.get());
  //   return nullptr;
  // }
  // if (!ptr->HasAnyOptionOf(std::span(kOptProbType))) {
  //   PrintUsage(argv[0], ptr.get());
  //   return nullptr;
  // }
  return ptr;
}

std::string WrapText(const std::string& text,
                     size_t width,
                     size_t initial_indent = 0) {
  std::istringstream iss(text);
  std::ostringstream oss;

  size_t pos = initial_indent;
  std::string word;
  std::string line;
  std::getline(iss, line);
  bool continue_line = true;
  while (iss) {
    std::istringstream word_stream(line);
    while (word_stream >> word) {
      if (pos + word.size() <= width && continue_line) {
        oss << word << " ";
        pos += word.size() + 1;
      } else {
        oss << "\n" << std::string(initial_indent, ' ');
        oss << word << " ";
        pos = initial_indent + word.size() + 1;
        continue_line = true;
      }
    }

    if (iss) {
      std::getline(iss, line);
      if (!line.empty()) {
        continue_line = false;
      }
    }
  }
  return oss.str();
}

void PrintUsage(const Path& program_path, const InputParser* p) {
  std::ostringstream oss, temp;

  // Build usage line
  std::ostringstream usage_oss;
  usage_oss << program_path.filename().string() << " ";
  usage_oss << "[" << kOptProbType[0] << " ProblemType] ";
  usage_oss << "[" << kOptAlgoType[0] << " AlgoType] ";
  usage_oss << "[" << kOptFiles[0] << " file1.fasta file2.fasta] "
            << "[" << kOptVerbose[0] << "] "
            << "[" << kOptHelp[0] << "] "
            << "[" << kOptSave[0] << " result.tsv] "
            << "[" << kOptVersion[0] << "] ";
  usage_oss << "[" << kOptMcGap[0] << " gc.json] ";
  usage_oss << "[" << kOptAlph[0] << " alph.json "
            << "[" << kOptLeftGap[0] << " left.json | " << kOptRightGap[0]
            << " right.json | both]] "
            << "[" << kOptPath[0] << " path/to/swisstree]";

  std::string usage_str = usage_oss.str();
  std::cout << "Usage: " << WrapText(usage_str, 80, 7) << "\n";

  std::map<std::string_view, std::string> args = {
      {kOptProbType[0], "ProblemType"},
      {kOptAlgoType[0], "AlgoType"},
      {kOptMcGap[0], "gc.json"},
      {kOptAlph[0], "alph.json"},
      {kOptLeftGap[0], "left.json"},
      {kOptRightGap[0], "right.json"},
      {kOptFiles[0], "file1.fasta file2.fasta"},
      {kOptSave[0], "result.tsv"},
      {kOptVerbose[0], ""},
      {kOptHelp[0], ""},
      {kOptVersion[0], ""},
      {kOptPath[0], "path/to/swisstree"},
  };

  // ProblemType description
  temp.str("");
  temp << "Specifies the problem's type:\n";
  using lcs_solver::problems::ProblemFactory;
  for (const auto& type : ProblemFactory::kAvailable) {
    temp << ProblemFactory::GetName(type) << " ";
  }
  std::string problem_type_msg = temp.str();

  // AlgoType description
  temp.str("");
  temp << "Specifies the algorithm:\n";
  using lcs_solver::algorithms::AlgoFactory;
  for (const auto& type : AlgoFactory::GetLLCSAvailable()) {
    temp << "  " << AlgoFactory::GetName(type) << " ";
  }
  std::string algo_type_msg = temp.str();

  std::map<std::string_view, std::string> msg = {
      {kOptProbType[0], problem_type_msg},
      {kOptAlgoType[0], algo_type_msg},
      {kOptMcGap[0], "JSON file with gap constraints for MC algorithms."},
      {kOptAlph[0], "Enable symbol-dependent constraints. Requires -l/-r."},
      {kOptLeftGap[0], "Left gap constraint tuple (sorted by alph.json)."},
      {kOptRightGap[0], "Right gap constraints tuple (sorted by alph.json)."},
      {kOptFiles[0], "Input FASTA files."},
      {kOptSave[0], "Save results to TSV file."},
      {kOptVerbose[0], "Enable verbose output."},
      {kOptHelp[0], "Show this help message."},
      {kOptVersion[0], "Display version information."},
      {kOptPath[0], "Path to SwissTree dataset."},
  };

  // size_t max_len = 0;
  // for (const auto& [flag, placeholder] : args) {
  //   max_len = std::max(max_len, flag.size() + placeholder.size() + 3);
  // }

  std::cout << "Options:\n";
  for (const auto& [flag, placeholder] : args) {
    std::string opt_line = "  " + std::string(flag);
    if (!placeholder.empty()) opt_line += " " + placeholder;
    opt_line += " : ";
    const std::string& desc = msg.at(flag);

    std::string wrapped_desc = WrapText(desc, 80, opt_line.size());
    std::cout << opt_line << wrapped_desc << "\n";
  }

  oss << "\nBenchmark of " << lcs_solver::Name()
      << " based on the SwissTree Dataset:\n";
  std::cout << oss.str();
  PrintSwissTreeInfo(p);
}

void PrintSwissTreeInfo(const InputParser* p) {
  for (const auto& [key, pair] : GetSwissTreeData(p)) {
    std::cout << key << "\t";
    std::cout << std::left << std::setw(50) << std::setfill(' ');
    std::cout << pair.first << "\t";
    std::cout << pair.second.size() << std::endl;
  }
}
void PrintVersion(const Path& program_path) {
  std::cout << program_path.filename().string();
  std::cout << " Version 1.0.2 - (based on " << lcs_solver::Name() << ")";
}

bool GetVerbose(const InputParser* p, const bool default_value) {
  return p ? p->HasOption("-v") : default_value;
}

std::pair<bool, Path> GetSaveInfo(const InputParser* p) {
  if (!p) return {false, ""};
  bool save = p->HasOption("-s");
  Path file = p->GetOption("-s");
  return {save, file};
}

std::unique_ptr<GapVector> GetGapVector(
    const InputParser* p, const std::span<const std::string_view> options) {
  if (!p->HasAnyOptionOf(options)) {
    return nullptr;
  }
  const Path file = p->GetFistOptionOf(options);
  std::ifstream f(file);
  if (!f) {
    throw std::runtime_error("Could not open file: " + file.string());
  }
  GapVector gc = nlohmann::json::parse(f);
  return std::make_unique<GapVector>(gc);
}

AlgoType GetAlgorithmType(const InputParser* p) {
  if (!p->HasAnyOptionOf(std::span(kOptAlgoType))) {
    return AlgoType::LLCS2_STD_FL;
    // throw std::runtime_error("No algorithm type was set");
  }
  const std::string& type = p->GetFistOptionOf(kOptAlgoType);
  return lcs_solver::algorithms::AlgoFactory::GetMapNameToAlgoType().at(type);
}

std::shared_ptr<SymbolVector> GetAlphabet(const InputParser* p) {
  if (!p->HasAnyOptionOf(kOptAlph)) {
    return nullptr;
  }
  const Path file = p->GetFistOptionOf(std::span(kOptAlph));
  std::ifstream f(file);
  if (!f) {
    throw std::runtime_error("Could not open file: " + file.string());
  }
  const std::vector<std::string> vec = nlohmann::json::parse(f);
  SymbolVector alph;
  for (const std::string& s : vec) {
    String temp = lcs_solver::util::to_String<Symbol>(s);
    alph.emplace_back(temp.front());
  }
  return std::make_shared<SymbolVector>(alph);
}

ProblemType GetProblemType(const InputParser* p) {
  if (!p->HasOption("-t")) {
    return ProblemType::LCS_Classic;
    // return lcs_solver::problems::ProblemFactory::ReadProblemType();
  }
  const std::string& type = p->GetFistOptionOf(kOptProbType);
  return lcs_solver::problems::ProblemFactory::GetNameTypeMap().at(type);
}

std::tuple<std::shared_ptr<SymbolVector>,
           std::unique_ptr<GapVector>,
           std::unique_ptr<GapVector>>
GetSigmaVectors(const InputParser* p) {
  std::shared_ptr<SymbolVector> alph = GetAlphabet(p);
  std::unique_ptr<GapVector> left = GetGapVector(p, std::span(kOptLeftGap));
  std::unique_ptr<GapVector> right = GetGapVector(p, std::span(kOptRightGap));
  return std::make_tuple(std::move(alph), std::move(left), std::move(right));
}

std::vector<PathPair> GetCombinations(const InputParser* parser) {
  std::vector<PathPair> result;
  for (const auto& opt : kOptFiles) {
    if (parser->HasOption(opt)) {
      const auto [file1, file2] = parser->GetPair(opt);
      result.emplace_back(file1, file2);
      return result;
    }
  }
  // No Files Flag
  const InfoMap& map = GetSwissTreeData(parser);
  // for (const auto& group_pair: GetFamilyCombinations(map)) {
  //   for (const auto& file_pair: GetPathCombinations(map, group_pair)) {
  //     result.emplace_back(file_pair);
  //   }
  // }
  // return result;
  return GetPathCombinations(map);
}

const InfoMap& GetSwissTreeData(const InputParser* p) {
  constexpr std::string_view default_path = "./swisstree";
  std::string path;
  if (p->HasAnyOptionOf(kOptPath)) {
    path = p->GetFistOptionOf(kOptPath);
  } else {
    path = default_path;
  }
  std::vector<string> short_names = {"ST001",
                                     "ST002",
                                     "ST003",
                                     "ST004",
                                     "ST005",
                                     "ST007",
                                     "ST008",
                                     "ST009",
                                     "ST010",
                                     "ST011",
                                     "ST012"};

  std::vector<string> gene_families = {
      "Popeye domain-containing protein family",
      "NOX 'ancestral-type' subfamily NADPH oxidases",
      "V-type ATPase beta subunit",
      "Serine incorporator family",
      "SUMF family",
      "Ribosomal protein S10/S20",
      "Bambi family",
      "Asterix family",
      "Cited family",
      "Glycosyl hydrolase 14 family",
      "Ant transformer protein"};

  const std::regex re(R"((\d+))");  // get the number in the short_name
  int last = 1;
  std::vector<std::vector<Path>> files{11};
  auto files_iter = files.begin();
  for (const auto& entry : std::filesystem::directory_iterator(path)) {
    std::string filename = entry.path().filename().string();
    if (std::smatch match; std::regex_search(filename, match, re)) {
      if (const int i = std::stoi(match[1]); i != last) {
        ++files_iter;
        last = i;
      }
      files_iter->push_back(entry.path());
    }
  }

  static std::map<string, std::pair<string, std::vector<Path>>> map;
  for (auto [name, desc, file] :
       std::views::zip(short_names, gene_families, files)) {
    map[name] = {desc, file};
  }
  return map;
}

std::vector<GroupPair> GetFamilyCombinations(const InfoMap& map) {
  static std::vector<GroupPair> pairs;
  if (!pairs.empty()) {
    return pairs; // avoid calculating the pairs twice
  }
  for (const string& s1 : std::views::keys(map)) {
    for (const string& s2 : std::views::keys(map)) {
      if (s1 <= s2) {
        pairs.emplace_back(s1, s2);
      }
    }
  }
  return pairs;
}

void AddPathCombinations(std::vector<PathPair>& v,
                         const InfoMap& map,
                         const GroupPair& pair) {
  std::vector<PathPair> pairs;
  const auto& [key1, key2] = pair;
  const std::vector<Path>& vec1 = map.at(key1).second;
  const std::vector<Path>& vec2 = map.at(key2).second;
  for (const Path& path1 : vec1) {
    for (const Path& path2 : vec2) {
      if (path1 < path2) {
        v.emplace_back(path1, path2);
      }
    }
  }
}

std::vector<PathPair> GetPathCombinations(const InfoMap& map,
                                          const GroupPair& pair) {
  std::vector<PathPair> result;
  AddPathCombinations(result, map, pair);
  return result;
}

std::vector<PathPair> GetPathCombinations(const InfoMap& map) {
  std::vector<PathPair> result;
  for (const auto& group_pair : GetFamilyCombinations(map)) {
    AddPathCombinations(result, map, group_pair);
  }
  return result;
}

lcs_solver::util::String From(const Path& file_path) {
  std::ostringstream oss;
  std::ifstream file(file_path);
  if (!file) {
    throw std::runtime_error("Unable to open file: " + file_path.string());
  }
  std::string line;
  std::ostringstream result;
  std::getline(file, line);  // Skip the first line (e.g., ">ST012_013")

  while (std::getline(file, line)) {
    result << line;
  }
  std::string str = result.str();
  return lcs_solver::util::to_String<lcs_solver::util::Symbol>(str);
}

void OutputParser::Add(const PathPair& pair, size_t llcs) {
  const auto& [file1, file2] = pair;
  map_[file1.filename().string()].emplace(file2.filename().string(), llcs);
}

void OutputParser::Print(std::ostream& os) {
  std::ostringstream oss;
  for (const auto& [key1, map1] : map_) {
    for (const auto& [key2, llcs] : map1) {
      oss << key1 << "\t" << key2 << "\t" << llcs << '\n';
    }
  }
  os << oss.str();
}

void OutputParser::Print(PathPair file_pair,
                         const size_t llcs,
                         std::ostream& os) {
  const auto [file1, file2] = file_pair;
  os << file1.filename().string() << "\t" << file2.filename().string() << "\t"
     << llcs << '\n';
}

void OutputParser::Save(const std::filesystem::path& file_path) {
  std::ofstream file(file_path, std::ios::out | std::ios::trunc);
  if (!file) {
    throw std::runtime_error("Unable to open file: " + file_path.string());
  }
  this->Print(file);
  file.close();
  std::cout << "Saved output to " << file_path << std::endl;
}

size_t ShowProgress(size_t counter, const size_t total) {
  using Clock = std::chrono::steady_clock;
  static Clock::time_point start;
  static Clock::time_point last_print;

  if (counter == 0) {
    start = Clock::now();
    last_print = start;
  }

  const auto now = Clock::now();
  const auto min_delta = std::chrono::seconds(2);

  if (now - last_print >= min_delta || counter == total || counter == 0) {
    last_print = now;

    constexpr int width = 50;
    const float progress = static_cast<float>(counter) / static_cast<float>(total);
    const int pos = static_cast<int>(progress * width);

    std::ostringstream oss;
    oss << "[";
    for (int i = 0; i < width; ++i) {
      if (i < pos)
        oss << "=";
      else if (i == pos)
        oss << ">";
      else
        oss << " ";
    }
    oss << "] ";

    std::ostringstream count_stream;
    count_stream << counter << "/" << total;
    oss << std::setw(15) << std::left << count_stream.str();

    const std::chrono::duration<double> elapsed = now - start;
    oss << std::fixed << std::setprecision(3) << elapsed.count() << " sec";

    std::string output = oss.str();
    if (constexpr size_t console_width = 100; output.size() < console_width) {
      output.append(console_width - output.size(), ' ');
    }
    std::cout << "\r" << output << std::flush;
  }

  return ++counter;
}

}  // namespace dna