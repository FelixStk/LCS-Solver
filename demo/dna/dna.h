#ifndef DEMO_DNA_DNA_H_
#define DEMO_DNA_DNA_H_

#include "OutputParser.h"

#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "algorithms/AlgoType.h"
#include "lcs_solver/problems/ProblemType.h"
#include "lcs_solver/util/CommonTypes.h"
#include "util/InputParser.h"

namespace dna {
using std::string;
using uint          = lcs_solver::util::uint;
using GapVector     = std::vector<std::pair<uint, uint>>;
using GapSpan       = std::span<const std::pair<uint, uint>>;
using Symbol = lcs_solver::util::Symbol;
using SymbolVector  = std::vector<lcs_solver::util::Symbol>;
using SymbolSpan    = std::span<const lcs_solver::util::Symbol>;

using AlgoType      = lcs_solver::algorithms::AlgoType;
using Path          = std::filesystem::path;
using PathPair      = std::pair<Path, Path>;
using GroupPair     = std::pair<string, string>;
using InfoMap       = std::map<string, std::pair<string, std::vector<Path>>>;
using ProblemType   = lcs_solver::problems::ProblemType;
using InputParser   = lcs_solver::util::InputParser;
using String        = lcs_solver::util::String;

constexpr std::array<std::string_view, 2> kOptProbType = {"-p", "--problem"};
constexpr std::array<std::string_view, 2> kOptAlgoType = {"-t", "--type"};
constexpr std::array<std::string_view, 2> kOptFiles = {"-f", "--files"};
constexpr std::array<std::string_view, 2> kOptMcGap = {"-g", "--mc"};
constexpr std::array<std::string_view, 2> kOptAlph = {"-a", "--alph"};
constexpr std::array<std::string_view, 2> kOptLeftGap = {"-l", "--left"};
constexpr std::array<std::string_view, 2> kOptRightGap = {"-r", "--right"};
constexpr std::array<std::string_view, 2> kOptSave = {"-s", "--save"};
constexpr std::array<std::string_view, 2> kOptVerbose = {"-v", "--verbose"};
constexpr std::array<std::string_view, 2> kOptHelp = {"-h", "--help"};
constexpr std::array<std::string_view, 2> kOptVersion = {"-V", "--version"};
constexpr std::array<std::string_view, 2> kOptPath = {"-P", "--path"};

void Run(int argc, char* argv[]);
std::unique_ptr<OutputParser> DoMainLoop(const std::vector<PathPair> & comb_vec,
                                         const ProblemType& prob_t,
                                         const AlgoType& algo_t,
                                         const GapVector* gc,
                                         const SymbolVector* alph,
                                         const GapVector* l,
                                         const GapVector* r,
                                         bool verbose);
std::tuple<GapSpan, SymbolSpan, GapSpan, GapSpan> GetCorrectedParams(std::span<const String> s, const GapVector* gc,
                                         const SymbolVector* alph,
                                         const GapVector* l,
                                         const GapVector* r);

// Program Start and options
std::unique_ptr<InputParser> HandleArgs(int argc, char* argv[]);
void PrintUsage(const Path& program_path, const InputParser* p);
void PrintSwissTreeInfo(const InputParser* p);
void PrintVersion(const Path& program_path);

// Configure Main Loops
bool GetVerbose(const InputParser* p, bool default_value = false);
std::pair<bool, Path> GetSaveInfo(const InputParser* p);
ProblemType GetProblemType(const InputParser* p);
AlgoType GetAlgorithmType(const InputParser* p);

std::shared_ptr<SymbolVector> GetAlphabet(const InputParser* p);
std::unique_ptr<GapVector> GetGapVector(const InputParser* p, std::span<const std::string_view> options = std::span(kOptMcGap));
std::unique_ptr<GapVector> GetMapValues(const InputParser* p, std::span<const std::string_view> options);
std::tuple<std::shared_ptr<SymbolVector>, std::unique_ptr<GapVector>, std::unique_ptr<GapVector>> GetSigmaVectors(const InputParser* p);


std::vector<PathPair> GetCombinations(const InputParser* parser);

// Functions to create input sequence combinations
const InfoMap& GetSwissTreeData(const InputParser* p);
std::vector<GroupPair> GetFamilyCombinations(const InfoMap&);
void AddPathCombinations(std::vector<PathPair>& v, const InfoMap&, const GroupPair&);
std::vector<PathPair> GetPathCombinations(const InfoMap&, const GroupPair&);
std::vector<PathPair> GetPathCombinations(const InfoMap&);
std::vector<PathPair> GetAllPathCombinations(const InfoMap&);
lcs_solver::util::String From(const Path&);

size_t ShowProgress(size_t counter, size_t total);

}  // namespace dna

#endif  // DEMO_DNA_DNA_H_
