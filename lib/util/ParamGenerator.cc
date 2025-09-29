/*******************************************************************************
 * @file ParamGenerator.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Generates test parameters for testing
 ******************************************************************************/
#include "util/ParamGenerator.h"

#include <cassert>
#include <list>
#include <set>
#include <sstream>

#include "algorithms/LLCS/LLCS_STD_FL.h"
#include "algorithms/LLCS/LLCS2_MC.h"
#include "algorithms/solutions/UnsignedSolution.h"
#include "algorithms/solutions/Vector3DSolution.h"
#include "algorithms/solutions/EmptySolution.h"
#include "constraints/BaseConstraint.h"
#include "constraints/local/Constraint_MC.h"
#include "constraints/local/Constraint_MC_1C.h"
#include "constraints/local/Constraint_MC_INC.h"
#include "constraints/local/Constraint_MC_O1C.h"
#include "constraints/local/Constraint_MC_O1C_SYNC.h"
#include "constraints/local/Constraint_Sigma.h"
#include "constraints/local/Constraint_Sigma_R.h"
#include "constraints/local/Constraint_Sigma_L.h"
#include "constraints/global/Constraint_BR.h"
#include "util/Logger.hpp"

namespace{
using ::lcs_solver::algorithms::AlgoType;
using ::lcs_solver::algorithms::BaseAlgorithm;
using ::lcs_solver::algorithms::BaseSolution;
using ::lcs_solver::algorithms::solutions::EmptySolution;
using ::lcs_solver::algorithms::solutions::UnsignedSolution;
using ::lcs_solver::constraints::ConstraintType;
using ::lcs_solver::constraints::BaseConstraint;
using ::lcs_solver::constraints::local::Constraint_MC;
using ::lcs_solver::constraints::local::Constraint_MC_INC;
using ::lcs_solver::constraints::local::Constraint_MC_1C;
using ::lcs_solver::constraints::local::Constraint_MC_O1C;
using ::lcs_solver::constraints::local::Constraint_MC_O1C_SYNC;
using ::lcs_solver::constraints::local::Constraint_Sigma;
using ::lcs_solver::constraints::local::Constraint_Sigma_R;
using ::lcs_solver::constraints::local::Constraint_Sigma_L;
using ::lcs_solver::constraints::global::Constraint_BR;

using ConstraintMap = ::lcs_solver::algorithms::BaseAlgorithm::ConstraintMap;
using Pair = ::lcs_solver::util::ParamGenerator::Pair;
using RndDist = std::uniform_int_distribution<size_t>;
}


namespace lcs_solver::util {

//=== Public Methods ===========================================================

/*******************************************************************************
 * ParamGenerator::genWithUniqSol generates Parameters such that the strings
 * contain exactly one unique longest common subsequence, that is contained in
 * each of the exactly one time.
 * @param type ConstraintType to add into the ConstraintMap
 * @param initSeed unsigned int used to set up the random number generator
 * @param nParams size_t number of problems to create
 * @param l std::vector<Pair> such that string[k].size is in the interval l[k]
 * @param solLen Pair that specifies bounds for the randomly chosen llcs
 * @param c std::vector<Symbol> c[0] is a symbol shared across all strings,
 *        otherwise c[k] is a unique symbol for filling the kth string
 * @param use_unsigned_sol Flag is true if the solution is an UnsignedSolution
 *        otherwise the solution is of Type Points. (default = true)
 * @param one_based_positions Flag is true if the positions in the solution are
 *        one-based (default = true)
 * @return std::vector<AlgoParam>
 ******************************************************************************/
std::vector<AlgoParam> ParamGenerator::genWithUniqSol(
    ConstraintType type,
    unsigned int initSeed,
    size_t nParams,
    const std::vector<Pair> &l,
    const Pair &solLen,
    const std::vector<Symbol> &c,
    const bool use_unsigned_sol,
    const bool one_based_positions
) {
  std::vector<AlgoParam> params;
  params.reserve(nParams);
  for (size_t i = 0; i < nParams; ++i) {
    const size_t seed = initSeed + i;
    std::mt19937 gen(seed);

    auto lengths = genStrLengths(gen, l);
    std::ranges::sort(lengths);
    const auto llcs = chooseLLCS(gen, solLen, lengths);

    std::string name = "Rnd_" + std::to_string(seed);
    auto [strings, mat] = genStrLLCS(gen, llcs, lengths, c, one_based_positions);

    // Set up Constraints
    auto gapBounds = getUniqueGaps(strings, llcs, c[0]);
    auto map = BaseAlgorithm::ConstraintMap();
    switch (type) {
      case ConstraintType::Empty: break;
      case ConstraintType::MC: {
        map[type] = std::make_shared<Constraint_MC>(gapBounds);
        break;
      }
      case ConstraintType::MC_INC: {
        relaxToIncProperty(gapBounds);
        map[type] = std::make_shared<Constraint_MC_INC>(gapBounds);
        break;
      }
      case ConstraintType::MC_1C: {
        relaxToMC1CProperty(gapBounds);
        map[type] = std::make_shared<Constraint_MC_1C>(gapBounds);
        break;
      }
      case ConstraintType::MC_O1C: {
        map[type] = std::make_shared<Constraint_MC_O1C>(gapBounds);
        break;
      }
      case ConstraintType::MC_O1C_SYNC: {
        relaxToSYNCProperty(gapBounds);
        map[type] = std::make_shared<Constraint_MC_O1C_SYNC>(gapBounds);
        break;
      }
      case ConstraintType::SIGMA_R: {
        relaxToMC1CProperty(gapBounds);
        gapBounds.resize(1);
        std::vector<Symbol> s = {c[0]};
        for (size_t p = 1; p < c.size(); ++p) {
          s.push_back(c[p]);
          gapBounds.push_back(genRelaxedPair(strings));
        }
        //auto relaxed = std::vector<Pair>(s.size(), genRelaxedPair(strings));
        map[type] = std::make_shared<Constraint_Sigma_R>(s, gapBounds);
        break;
      }
      case ConstraintType::SIGMA_L: {
        relaxToMC1CProperty(gapBounds);
        gapBounds.resize(1);
        std::vector<Symbol> s = {c[0]};
        for (size_t p = 1; p < c.size(); ++p) {
          s.push_back(c[p]);
          gapBounds.push_back(genRelaxedPair(strings));
        }
        map[type] = std::make_shared<Constraint_Sigma_L>(s, gapBounds);
        break;
      }
      case ConstraintType::SIGMA: {
        relaxToMC1CProperty(gapBounds);
        gapBounds.resize(1);
        std::vector<Symbol> s = {c[0]};
        for (size_t p = 1; p < c.size(); ++p) {
          s.push_back(c[p]);
          gapBounds.push_back(genRelaxedPair(strings));
        }
        map[type] = std::make_shared<Constraint_Sigma>(s, gapBounds, gapBounds);
        break;
      }
      default: {
        Logger::Warning() << "genWithUniqSol: ConstraintType == "
                          << std::to_string(static_cast<unsigned char>(type))
                          << " is not implemented (nothing is added to the map)";
      }
    }

    // Setup StrPtrVec
    AlgoParam::StrPtrVector spv;
    for (auto& s : strings) {
      auto ptr = std::make_shared<const String>(s);
      spv.push_back(ptr);
    }

    // Setup Solutions
    BaseSolution *solution;
    if (use_unsigned_sol) {
      solution = new UnsignedSolution(llcs);
    }else {
      using algorithms::solutions::Vector3DSolution;
      std::vector points_vec = {
        Vector3DSolution::Points(spv, mat)
      };
      solution = new Vector3DSolution(points_vec);
    }
    params.emplace_back(name, spv, map, solution);
  }

  return params;
}

/*******************************************************************************
 * ParamGenerator::genWithStdSol
 * @param type ConstraintType to add into the ConstraintMap
 * @param initSeed unsigned int used to set up the random number generator
 * @param nParams size_t number of problems to create
 * @param l std::vector<Pair> such that string[k].size is in the interval l[k]
 * @param alphabet std::vector<Symbol> alphabet for the strings
 * @return std::vector<AlgoParam> with relaxed constraint such the constraint is
 * reduces to an equivalent of the standard lcs problem that can be solved with
 * the folklore algorithm
 ******************************************************************************/
std::vector<AlgoParam> ParamGenerator::genWithStdSol(
    ParamGenerator::ConstraintType type,
    const unsigned int initSeed,
    const size_t nParams,
    const std::vector<Pair> &l,
    const std::vector<Symbol> &alphabet
) {
  std::vector<AlgoParam> params;
  params.reserve(nParams);
  for (size_t i = 0; i < nParams; ++i) {
    const size_t seed = initSeed + i;
    std::mt19937 gen(seed);
    auto lengths = genStrLengths(gen, l);
    std::ranges::sort(lengths);
    std::string name = "Rnd_" + std::to_string(seed);;
    StringVec strings = genRndStrings(gen, lengths, alphabet);
    auto map = genRelaxedToStdMap(type, strings, alphabet);
    auto solution = genStdFLSol(strings);
    params.emplace_back(std::move(name),
                        std::move(strings),
                        std::move(map),
                        solution.release());
  }
  return params;
}

/*******************************************************************************
 * ParamGenerator::genWithMCSol
 * @param type ConstraintType to add into the ConstraintMap
 * @param initSeed unsigned int used to set up the random number generator
 * @param nParams size_t number of problems to create
 * @param l std::vector<Pair> such that string[k].size is in the interval l[k]
 * @param b Pair fixes the interval from with gaps are drawn
 * @param alphabet std::vector<Symbol> alphabet for the strings
 * @return std::vector<AlgoParam>
 ******************************************************************************/
std::vector<AlgoParam> ParamGenerator::genWithMCSol(
    ConstraintType type,
    unsigned int initSeed,
    size_t nParams,
    const std::vector<Pair> &l,
    Pair b,
    const std::vector<Symbol> &alphabet
) {
  std::vector<AlgoParam> params;
  params.reserve(nParams);
  for (size_t i = 0; i < nParams; ++i) {
    const size_t seed = initSeed + i;
    std::mt19937 gen(seed);
    auto lengths = genStrLengths(gen, l);
    std::ranges::sort(lengths);
    std::string name = "Rnd_" + std::to_string(seed);;
    StringVec strings = genRndStrings(gen, lengths, alphabet);
    auto map = genRelaxedToStdMap(type, strings, alphabet);
    std::unique_ptr<BaseSolution> solution;
    if (const auto *gc = dynamic_cast<Constraint_MC *>(map[type].get()))
      solution = genMCSol(strings, gc->GetGapVector());
    else
      solution = std::make_unique<EmptySolution>();
    params.emplace_back(std::move(name),
                        std::move(strings),
                        std::move(map),
                        solution.release());
  }
  return params;
}

/*******************************************************************************
 * ParamGenerator::genWithoutSol
 * @param t ConstraintType to add into the ConstraintMap
 * @param seed unsigned int used to set up the random number generator
 * @param l std::vector<Pair> such that string[k].size is in the interval l[k]
 * @param b Pair that represents the closed interval that gaps are drawn from
 * @param alphabet std::vector<Symbol> for drawing random Symbols
 * @return AlgoParam
 ******************************************************************************/
AlgoParam ParamGenerator::genWithoutSol(
    ConstraintType t,
    unsigned int seed,
    const std::vector<Pair> &l,
    const Pair &b,
    const std::vector<Symbol> &alphabet) {
  std::mt19937 gen(seed);
  auto lengths = genStrLengths(gen, l);
  std::ranges::sort(lengths);
  std::string name = "Rnd_" + std::to_string(seed);;
  StringVec strings = genRndStrings(gen, lengths, alphabet);
  auto map = genRndMap(gen, t, strings, b, alphabet);
  return {name, strings, map, nullptr};
}

/*******************************************************************************
 * ParamGenerator::genWithoutSol
 * @param t ConstraintType to add into the ConstraintMap
 * @param seed unsigned int used to set up the random number generator
 * @param nParams size_t number of problems to create
 * @param l std::vector<Pair> such that string[k].size is in the interval l[k]
 * @param b Pair that represents the closed interval that gaps are drawn from
 * @param alphabet std::vector<Symbol> for drawing random Symbols
 * @return std::vector<AlgoParam>
 ******************************************************************************/
std::vector<AlgoParam> ParamGenerator::genWithoutSol(
    ConstraintType t,
    unsigned int seed,
    size_t nParams,
    const std::vector<Pair> &l,
    const Pair &b,
    const std::vector<Symbol> &alphabet) {
  std::vector<AlgoParam> params;
  params.reserve(nParams);
  for (size_t i = 0; i < nParams; ++i) {
    params.push_back(
        genWithoutSol(
            t,
            seed+i,
            l,
            b,
            alphabet));
  }
  return params;
}

//==== Private Methods =========================================================

/*******************************************************************************
 * ParamGenerator::genStrLengths
 * @param gen std::mt19937 random number generator
 * @param l std::vector<Pair> l[k]= {min, max} bounds the length of str[k]
 * @return std::vector<size_t> res with l[k].first <= res[k] <= l[k].second
 ******************************************************************************/
std::vector<size_t> ParamGenerator::genStrLengths(
    std::mt19937 &gen,
    const std::vector<Pair> &l
) {
  std::vector<size_t> res;
  res.reserve(l.size());
  for (auto &k : l) {
    RndDist distLen(k.first, k.second);
    res.emplace_back(distLen(gen));
  }
  return res;
}

/*******************************************************************************
 * ParamGenerator::chooseLLCS
 * @param gen std::mt19937 random number generator
 * @param p Pair that bounds the length of the llcs
 * @param len std::vector<size_t> lengths of the strings in a problem
 * @return size_t res : res <= min_k(len[k]) and p.first <= res <= p.second
 ******************************************************************************/
size_t ParamGenerator::chooseLLCS(
    std::mt19937 &gen,
    const Pair &p,
    const std::vector<size_t> &len
) {
  size_t minLength = *std::ranges::min_element(len);
  RndDist distLLCS(std::min(p.first, minLength),
                   std::min(p.second, minLength));
  return distLLCS(gen);
}

/*******************************************************************************
 * ParamGenerator::genStrLLCS
 * @param gen std::mt19937 random number generator
 * @param llcs size_t the length of the llcs in the strings to be generated
 * @param length std::vector<size_t> contains the lengths of the strings to gen
 * @param c std::vector<Symbol> c[0] is a symbol shared across all strings
 *        for k>0 c[k] is a unique symbol for filling the kth string
 * @param one_based_positions bool if true, the positions in the strings are
 *        one-based, otherwise they are zero-based
 * @return StringVec s with s[k]=shuffles(c[0]^llcs + c[k]^length)
 ******************************************************************************/
std::pair<ParamGenerator::StringVec, ParamGenerator::Points> ParamGenerator::genStrLLCS(
    std::mt19937 &gen,
    const size_t llcs,
    const std::vector<size_t> &length,
    const std::vector<Symbol> &c,
    const bool one_based_positions
) {
  StringVec string_vec;
  const size_t n = length.size();
  string_vec.reserve(n);
  for (size_t k = 0; k < n; ++k) {
    String str(llcs, c[0]);
    str.append(length[k] - llcs, c[k + 1]);
    std::ranges::shuffle(str, gen);
    string_vec.emplace_back(str);
  }

  Points point_vec(llcs, Point(n, 0));
  for (uint k = 0; k <string_vec.size(); ++k) {
    const String &s = string_vec[k];
    uint i = 0;
    for (uint pos = 0; pos < s.length(); ++pos) {
      if (s[pos] == c[0]) {
        if (!one_based_positions) {
          point_vec[i][k] = pos;
        } else {
          point_vec[i][k] = pos + 1;
        }
        ++i;
      }
    }
    Point p(n, 0);
  }
  return {string_vec, point_vec};
}

/*******************************************************************************
 * ParamGenerator::sortByLength
 * @param vec StringVec& to sort
 ******************************************************************************/
void ParamGenerator::sortByLength(StringVec &vec) {
  std::ranges::sort(vec, [](const auto &a, const auto &b) {
    return a.size() < b.size();
  });
}

/*******************************************************************************
 * ParamGenerator::getUniqueGaps
 * @param strings StringVec with a unique llcs in each string
 * @param llcs size_t the length of the unique lcs
 * @param c Symbol used in the lcs
 * @return std::vector<Pair> containing the lower- and upper bound for the gaps
 ******************************************************************************/
std::vector<Pair> ParamGenerator::getUniqueGaps(
    const ParamGenerator::StringVec &strings,
    size_t llcs,
    Symbol c
) {
  auto res = std::vector<Pair>();
  size_t minStrLen = strings.front().size();
  res.reserve(minStrLen == 0 ? 0 : minStrLen);
  // gaps before the first (and last) match are ignored
  for (const auto &s : strings) {
    size_t prevMatchIdx = 0, matchesSeen = 0;
    for (size_t i = 0; i < s.size(); ++i) {
      if (s[i] == c) {
        matchesSeen++;
        if (matchesSeen > 1) {
          size_t gapLength = i - prevMatchIdx - 1;
          if (res.size() == llcs - 1)
            expandInterval(gapLength, res[matchesSeen - 2]);
          else
            res.emplace_back(gapLength, gapLength); // s = strings[0]
        }
        prevMatchIdx = i;
      }
    }
  }// end foreach s in string
  Pair filler = genRelaxedPair(strings);
  while (res.size() < strings.front().size() - 1) {
    res.emplace_back(filler);
  }
  return res;
}

/*******************************************************************************
 * ParamGenerator::expandInterval
 * @param value size_t that should be in the interval after the func is called
 * @param interval std::pair<size_t, size_t> that might need to be changed
 ******************************************************************************/
void ParamGenerator::expandInterval(
    size_t value,
    std::pair<size_t, size_t> &interval
) {
  if (value < interval.first)
    interval.first = value;
  if (value > interval.second)
    interval.second = value;
}

/*******************************************************************************
 * ParamGenerator::addIncProperty modifies a vector of gaps such that it is
 * increasing. Increasing means \f$gap[i] \subseteq gap[i+1]\f$ for all suitable i
 * @param gaps std::vector<Pair> to modify
 * @note In general, applying this function to change gaps for a mc constraint
 * might change the solution. But if there is only one unique lcs in each
 * string, it doesn't change the solution (because it weakens the gap tuple like
 * in a relaxation of the constraint to the normal lcs problem without affecting
 * the solution).
 ******************************************************************************/
void ParamGenerator::relaxToIncProperty(std::vector<Pair> &gaps) {
  auto subset = [](const Pair &a, const Pair &b) -> bool {
    return b.first <= a.first && a.second <= b.second; // `a` is a subset of `b`
  };

  auto curr = gaps.begin();
  for (auto next = curr + 1; next != gaps.end(); ++curr, ++next) {
    if (!subset(*curr, *next)) {
      next->first = std::min<size_t>(curr->first, next->first);
      next->second = std::max<size_t>(curr->second, next->second);
    }
  }
}

/*******************************************************************************
 * ParamGenerator::relaxToMC1CProperty
 * @param gaps std::vector<Pair> to modify
 * @note must be called on strings that have a unique shared lcs
 ******************************************************************************/
void ParamGenerator::relaxToMC1CProperty(std::vector<Pair> &gaps) {
  size_t lower = std::numeric_limits<size_t>::max();
  size_t upper = std::numeric_limits<size_t>::min();
  for (auto [a, b] : gaps) {
    lower = std::min<size_t>(lower, a);
    upper = std::max<size_t>(upper, b);
  }
  for (auto &gap : gaps)
    gap = {lower, upper};
}

/*******************************************************************************
 * ParamGenerator::relaxToSYNCProperty
 * @param gap std::vector<Pair> to modify
 * @note must be called on strings that have a unique shared lcs
 ******************************************************************************/
void ParamGenerator::relaxToSYNCProperty(std::vector<Pair> &gap) {
  bool loop = true;
  while (loop) {
    loop = false;
    for (size_t i = 0; i < gap.size(); ++i) {
      for (size_t j = i + 1; j < gap.size(); ++j) {
        if (gap[i].first == gap[j].first && gap[i].second == gap[j].second) {
          for (util::uint e = 0; e + i < gap.size() && e + j < gap.size(); e++) {
            if (!subset(gap[i + e], gap[j + e])) {
              gap[j + e] = {
                  std::min<size_t>(gap[j + e].first, gap[i + e].first),
                  std::max<size_t>(gap[j + e].second, gap[i + e].second)
              };
              i = gap.size();
              j = gap.size();
              loop = true;
            }
          } // for e
        }
      } // for j
    } // for i
  } // while
}

/*******************************************************************************
 * ParamGenerator::subset
 * @param a std::pair<size_t,size_t> closed interval
 * @param b std::pair<size_t,size_t> closed interval
 * @return true if 'a' is a subset of 'b', otherwise false
 ******************************************************************************/
bool ParamGenerator::subset(const Pair &a, const Pair &b) {
  return b.first <= a.first && a.second <= b.second; // '`a`' is a subset of `b`
}

/*******************************************************************************
 * ParamGenerator::genRndStrings
 * @param gen std::mt19937 random number generator
 * @param lengths std::vector<size_t> containing the lengths of the Strings
 * @param alphabet std::vector<Symbol> the alphabet of the Strings to generate
 * @return StringVec containing the generated Strings
 ******************************************************************************/
ParamGenerator::StringVec ParamGenerator::genRndStrings(
    std::mt19937 &gen,
    const std::vector<size_t> &lengths,
    const std::vector<Symbol> &alphabet
) {
  StringVec vec;
  vec.reserve(lengths.size());
  RndDist dist(0, alphabet.size() - 1);
  for (auto len : lengths) {
    std::basic_ostringstream<Symbol> oss;
    for (size_t i = 0; i < len; ++i) {
      oss << alphabet[dist(gen)];
    }
    vec.emplace_back(oss.str());
  }
  return vec;
}

ConstraintMap ParamGenerator::genRelaxedToStdMap(
    ConstraintType type,
    const ParamGenerator::StringVec &strings
) {
  auto vec = strings;
  sortByLength(vec);
  std::vector<Symbol> alphabet = getSymbols(vec);
  return genRelaxedToStdMap(type, vec, alphabet);
}

ConstraintMap ParamGenerator::genRelaxedToStdMap(
    ConstraintType type,
    const ParamGenerator::StringVec &vec,
    const std::vector<Symbol> &alphabet
) {
  ConstraintMap map;
  switch (type) {
    case ConstraintType::MC: {
      map[type] = std::make_shared<Constraint_MC>(genRelaxedGaps(vec));
      break;
    }
    case ConstraintType::MC_INC: {
      map[type] = std::make_shared<Constraint_MC_INC>(genRelaxedGaps(vec));
      break;
    }
    case ConstraintType::MC_1C: {
      map[type] = std::make_shared<Constraint_MC_1C>(genRelaxedGaps(vec));
      break;
    }
    case ConstraintType::MC_O1C: {
      map[type] = std::make_shared<Constraint_MC_O1C>(genRelaxedGaps(vec));
      break;
    }
    case ConstraintType::MC_O1C_SYNC: {
      map[type] = std::make_shared<Constraint_MC_O1C_SYNC>(genRelaxedGaps(vec));
      break;
    }
    case ConstraintType::SIGMA: {
      auto gaps = std::vector<Pair>(alphabet.size(), genRelaxedPair(vec));
      map[type] = std::make_shared<Constraint_Sigma>(alphabet, gaps, gaps);
      break;
    }
    case ConstraintType::SIGMA_L: {
      auto gaps = std::vector<Pair>(alphabet.size(), genRelaxedPair(vec));
      map[type] = std::make_shared<Constraint_Sigma_L>(alphabet, gaps);
      break;
    }
    case ConstraintType::SIGMA_R: {
      auto gaps = std::vector<Pair>(alphabet.size(), genRelaxedPair(vec));
      map[type] = std::make_shared<Constraint_Sigma_R>(alphabet, gaps);
      break;
    }
    case ConstraintType::Empty:break;
    case ConstraintType::BR:
    case ConstraintType::STRINGS_K:
    case ConstraintType::STRINGS_2:
    case ConstraintType::LLCS_KNOWN:
    case ConstraintType::CONST_SIG:

    default: {
      Logger::Warning() << "genRelaxedToStdMap: ConstraintType == "
                        << static_cast<unsigned int>(type)
                        << " is not implemented (nothing is added to the map)\n";
    }
  }
  return map;
}

/*******************************************************************************
 * ParamGenerator::genRelaxedGaps
 * @param s StringVec with strings sorted by length (ascending)
 * size_t len size_t length of the vector to generate
 * @return std::vector<Pair> with
 ******************************************************************************/
std::vector<Pair> ParamGenerator::genRelaxedGaps(const StringVec &s) {
  std::vector<Pair> vec = {};
  if (s.empty())
    return std::vector<Pair>{};
  const size_t minStrLen = s.front().size();
  const size_t numOfGaps = minStrLen > 0 ? minStrLen - 1 : 0;
  const Pair relaxedGap = genRelaxedPair(s);
  return std::vector(numOfGaps, relaxedGap);
}

/*******************************************************************************
 * ParamGenerator::genRelaxedPair
 * @param strings StringVec problems strings sorted by length (ascending)
 * @return Pair (0, max(|s|) )
 ******************************************************************************/
Pair ParamGenerator::genRelaxedPair(const StringVec &strings) {
  if (strings.empty())
    return {0, 0};
  size_t maxStrLen = strings.back().size();
  return {0, maxStrLen};
}

/*******************************************************************************
 * aramGenerator::getSymbols
 * @param strings StringVec contains the Strings of the problem
 * @return std::vector<Symbol> with the symbols used in the problem
 ******************************************************************************/
std::vector<Symbol> ParamGenerator::getSymbols(const StringVec &strings) {
  std::set<Symbol> set;
  for (const auto &string : strings)
    for (const auto symbol : string)
      set.insert(symbol);
  return {set.begin(), set.end()};
}

/*******************************************************************************
 * ParamGenerator::genStdFLSol
 * @param vec StringVec that contains the strings of the problem
 * @return std::unique_ptr<BaseSolution> with the solution from LLCS_STD_FL
 ******************************************************************************/
std::unique_ptr<BaseSolution> ParamGenerator::genStdFLSol(const StringVec &vec) {
  ConstraintMap emptyMap = {};
  auto spv = AlgoParam::convertToPtr(vec);
  auto solver = new algorithms::llcs::LLCS_STD_FL(spv, emptyMap);
  auto solution = solver->query();
  delete solver;
  return solution;
}

/*******************************************************************************
 * ParamGenerator::genMCSol
 * @param strings StringVec that contains the strings of the problem
 * @param gaps std::vector<Pair> the gap mc tuple to use when solving
 * @return std::unique_ptr<BaseSolution> from LLCS2_MC
 ******************************************************************************/
std::unique_ptr<BaseSolution> ParamGenerator::genMCSol(
    const ParamGenerator::StringVec &strings,
    const std::vector<Pair> &gaps
) {
  ConstraintMap map = {};
  auto spv = AlgoParam::convertToPtr(strings);
  map[ConstraintType::MC] = std::make_shared<Constraint_MC>(gaps);
  auto solver = new algorithms::llcs::LLCS2_MC(spv, map);
  auto solution = solver->query();
  delete solver;
  return solution;
}

//=== Private Methods for genWithMCSol =========================================

ConstraintMap ParamGenerator::genRndMap(
    std::mt19937 &gen,
    ConstraintType type,
    const ParamGenerator::StringVec &str,
    const Pair &bound,
    const std::vector<Symbol> &alphabet
) {
  ConstraintMap map;
  size_t numGaps = 0;
  if (!str.empty() && str[0].size() >= 2) // assuming str[0] is the smallest
    numGaps = str[0].size() - 1;
  switch (type) {
    case ConstraintType::Empty: break;
    case ConstraintType::MC: {
      std::vector<Pair> gaps = genRndGaps(gen, numGaps, bound);
      map[type] = std::make_shared<Constraint_MC>(gaps);
      break;
    }
    case ConstraintType::MC_INC: {
      std::vector<Pair> gaps = genRndGaps(gen, numGaps, bound);
      relaxToIncProperty(gaps);
      map[type] = std::make_shared<Constraint_MC_INC>(gaps);
      break;
    }
    case ConstraintType::MC_1C: {
      std::vector<Pair> gaps = genRndGaps(gen, numGaps, bound);
      relaxToMC1CProperty(gaps);
      map[type] = std::make_shared<Constraint_MC_1C>(gaps);
      break;
    }
    case ConstraintType::MC_O1C: {
      std::vector<Pair> gaps = genRndGaps(gen, numGaps, bound);
      map[type] = std::make_shared<Constraint_MC_O1C>(gaps);
      break;
    }
    case ConstraintType::MC_O1C_SYNC: {
      std::vector<Pair> gaps = genRndGaps(gen, numGaps, bound);
      relaxToSYNCProperty(gaps);
      map[type] = std::make_shared<Constraint_MC_O1C_SYNC>(gaps);
      break;
    }
    case ConstraintType::SIGMA: {
      std::vector<Pair> left = genRndGaps(gen, alphabet.size(), bound);
      std::vector<Pair> right = genRndGaps(gen, alphabet.size(), bound);
      map[type] = std::make_shared<Constraint_Sigma>(alphabet, left, right);
      break;
    }
    case ConstraintType::SIGMA_L: {
      std::vector<Pair> left = genRndGaps(gen, alphabet.size(), bound);
      map[type] = std::make_shared<Constraint_Sigma_L>(alphabet, left);
      break;
    }
    case ConstraintType::SIGMA_R: {
      std::vector<Pair> right = genRndGaps(gen, alphabet.size(), bound);
      map[type] = std::make_shared<Constraint_Sigma_R>(alphabet, right);
      break;
    }
    case ConstraintType::BR: {
      size_t maxLen = str.empty() ? 0 : str.back().size();
      RndDist dist(0, maxLen);
      map[type] = std::make_shared<Constraint_BR>(dist(gen));
    }
    case ConstraintType::STRINGS_K:
    case ConstraintType::STRINGS_2:
    case ConstraintType::LLCS_KNOWN:
    case ConstraintType::CONST_SIG:
    default: {
      Logger::Warning() << "genRndMap: ConstraintType == "
                        << static_cast<unsigned char>(type)
                        << " is not implemented (nothing is added to the map)";
    }
  }
  return map;
}

/*******************************************************************************
 * ParamGenerator::genRndGaps
 * @param gen std::mt19937 random number generator
 * @param n size_t number of gaps to generate
 * @param bounds Pair that defines the interval from with gaps are drawn
 * @return std::vector<Pair> of generated gaps
 ******************************************************************************/
std::vector<ParamGenerator::Pair> ParamGenerator::genRndGaps(
    std::mt19937 &gen,
    const size_t n,
    const ParamGenerator::Pair &bounds
) {
  std::vector<Pair> vec;
  vec.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    vec.emplace_back(genRndPair(gen, bounds));
  }
  return vec;
}

/*******************************************************************************
 * ParamGenerator::genRndPair
 * @param gen std::mt19937 random number generator
 * @param bounds interval from which subintervals are drawn
 * @return Pair {a,b} with a<=b and a<= `bounds.first` and b<= `bounds.second`
 ******************************************************************************/
Pair ParamGenerator::genRndPair(std::mt19937 &gen, const Pair &bounds) {
  RndDist dist(bounds.first, bounds.second);
  auto a = dist(gen);
  auto b = dist(gen);
  if (b < a)
    std::swap(a, b);
  return {a, b};
}

} // end of namespace