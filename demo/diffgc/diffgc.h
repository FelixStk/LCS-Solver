#ifndef DEMO_DIFFGC_UTILITY_H
#define DEMO_DIFFGC_UTILITY_H

#include "FileAdapter.h"

#include <lcs_solver/problems/BaseProblem.h>
#include <lcs_solver/util/AlgoTestParam.h>
#include <lcs_solver/util/CommonTypes.h>

#include <filesystem>
#include <memory>
#include <string>
#include <tuple>

#include "../dna/dna.h"

namespace diffgc {

enum class ChangeType : char{
  kRemove = '-',
  kKeep = '=',
  kAdd = '+'
};

using lcs_solver::util::InputParser;
using lcs_solver::problems::BaseProblem;

using DiffVector = std::vector<std::pair<FileAdapter::Symbol, ChangeType>>;
using Path = std::filesystem::path;
using ProblemPtr = std::unique_ptr<BaseProblem>;
using Symbol = lcs_solver::util::Symbol;
using String = lcs_solver::util::String;

void Run(int argc, char* argv[]);

std::tuple<Path, Path, std::unique_ptr<InputParser>> HandleOptions(int argc, char *argv[]);

std::unique_ptr<FileAdapter> SetInputStrings(BaseProblem * ptr, const Path &file1, const Path &file2, const InputParser * parser);
std::string GetBinName(char name[]);
void PrintVersion(char name[]);
void PrintUsage(char name[]);


void HandleConstraint(BaseProblem * ptr,const InputParser * parser, const FileAdapter* adapter);
void SetAlgorithm(BaseProblem * ptr);
void HandleSave(const BaseProblem * ptr, const InputParser * parser, const FileAdapter * adapter);

String RunAndGetLcs(const BaseProblem * ptr);
DiffVector CalcDiffVector(
  const lcs_solver::util::String& lcs,
  const BaseProblem::StrPtrVec& spv
);

void PrintDiffVector(const DiffVector& dv, const FileAdapter * adapter, bool translate_back=true);
void PrintUnifiedDiff(const DiffVector &dv, const FileAdapter *adapter,
                      bool translate_back = true);

std::string_view GetDescription();

lcs_solver::util::AlgoParam GenerateRndExample();
std::vector<lcs_solver::util::Symbol> GetAlphabet();

}// namespace diffgc
#endif //DEMO_DIFFGC_UTILITY_H
