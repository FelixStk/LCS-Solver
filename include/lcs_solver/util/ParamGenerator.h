#ifndef LCS_SOLVER_UTIL_PARAM_GENERATOR_H_
#define LCS_SOLVER_UTIL_PARAM_GENERATOR_H_

#include <random>
#include <string>
#include <utility>

#include "util/AlgoTestParam.h"
#include "algorithms/BaseSolution.h"
#include "algorithms/solutions/Vector3DSolution.h"
#include "algorithms/BaseAlgorithm.h"

namespace lcs_solver::util {

class ParamGenerator {
  using BaseSolution = algorithms::BaseSolution;
  using ConstraintType = constraints::ConstraintType;
  using ConstraintMap = algorithms::BaseAlgorithm::ConstraintMap;
  using StringVec = std::vector<String>;

 public:
  using uint = size_t;
  using Pair = std::pair<uint, uint>;
  using Point = std::vector<uint>;
  using Points = std::vector<Point>;

  // Generates AlgoParams such that in each string the lcs is unique
  static std::vector<AlgoParam> genWithUniqSol(
      ConstraintType t,
      unsigned int initSeed,
      size_t nParams,
      const std::vector<Pair> &l,
      const Pair &solLen,
      const std::vector<Symbol> &c,
      bool use_unsigned_sol = true,
      bool one_based_positions = true
  );

  // Generates AlgoParams with the constraint relaxed to the std lcs problem
  static std::vector<AlgoParam> genWithStdSol(
      ConstraintType t,
      unsigned int initSeed,
      size_t nParams,
      const std::vector<Pair> &l,
      const std::vector<Symbol> &alphabet
  );

  // Generates AlgoParams for MC constraints (uses LLCS2_MC for getting the sol)
  static std::vector<AlgoParam> genWithMCSol(
      ConstraintType t,
      unsigned int initSeed,
      size_t nParams,
      const std::vector<Pair> &l,
      Pair b,
      const std::vector<Symbol> &alphabet
  );

  // Generate AlgoParams without a solution
  static AlgoParam genWithoutSol(
      ConstraintType t,
      unsigned int initSeed,
      const std::vector<Pair> &l,
      const Pair &b,
      const std::vector<Symbol> &alphabet
  );

  static std::vector<AlgoParam> genWithoutSol(
      ConstraintType t,
      unsigned int initSeed,
      size_t nParams,
      const std::vector<Pair> &l,
      const Pair &b,
      const std::vector<Symbol> &alphabet
  );

  static ConstraintMap genRelaxedToStdMap(ConstraintType type, const StringVec &strings);

 private:
  // Methods added for genWithUniqSol
  static std::vector<size_t> genStrLengths(std::mt19937 &gen, const std::vector<Pair> &l);
  static size_t chooseLLCS(std::mt19937 &gen, const Pair &p, const std::vector<size_t> &len);
  static std::pair<StringVec, Points> genStrLLCS(std::mt19937 &gen, size_t llcs, const std::vector<size_t> &length, const std::vector<Symbol> &c, bool one_based_positions);
  static std::vector<Pair> getUniqueGaps(const StringVec &strings, size_t llcs, Symbol c);
  static void expandInterval(size_t value, std::pair<size_t, size_t> &interval);
  static void relaxToIncProperty(std::vector<Pair> &gaps);
  static void relaxToMC1CProperty(std::vector<Pair> &gaps);
  static void relaxToSYNCProperty(std::vector<Pair> &gaps);
  static bool subset(const Pair &a, const Pair &b);

  // Methods added for genWithStdSol
  static StringVec genRndStrings(std::mt19937 &gen, const std::vector<size_t> &length, const std::vector<Symbol> &alphabet);
  static ConstraintMap genRelaxedToStdMap(ConstraintType type, const StringVec &strings, const std::vector<Symbol> &alphabet);
  static std::vector<Pair> genRelaxedGaps(const StringVec &strings);
  static Pair genRelaxedPair(const StringVec &strings);
  static std::vector<Symbol> getSymbols(const StringVec &strings);
  static std::unique_ptr<BaseSolution> genStdFLSol(const StringVec &strings);

  // Methods added for genWithMCSol
  static std::unique_ptr<BaseSolution> genMCSol(const StringVec &strings, const std::vector<Pair> &gap);
  static ConstraintMap genRndMap(std::mt19937 &gen,ConstraintType type, const StringVec &strings, const Pair &b, const std::vector<Symbol> &alphabet);
  static std::vector<Pair> genRndGaps(std::mt19937 &gen, size_t n, const Pair &bounds);
  static Pair genRndPair(std::mt19937 &gen, const Pair &bounds);

  static void sortByLength(StringVec &vec);
};

}
#endif //LCS_SOLVER_UTIL_PARAM_GENERATOR_H_
