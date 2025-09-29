/*******************************************************************************
 * @file hardCoded_LLCS.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief Getter for hard coded LLCS AlgoParam Problems
 ******************************************************************************/

#include "hardCoded_LLCS.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"
#include "util/AlgoTestParam.h"
#include "util/CommonTypes.h"
#include "util/ParamGenerator.h"

#include <cassert>
#include <random>
#include <ranges>

namespace lcs_solver::testing::llcs {

/*******************************************************************************
 * Generates Identical Strings
 * @param n Number of Strings to generate
 * @param len Length of the Strings to generate
 * @return Generated Strings in a std::vector
 ******************************************************************************/
std::vector<util::String> genEqualStrings(const size_t n, const size_t len) {
  const util::Symbol c = util::to_String<util::Symbol>("a").at(0);
  const util::String s(len, c);
  return std::vector(n, s);
}

/*******************************************************************************
 * Generates Strings that are no Symbols
 * @param n Number of Strings to generate
 * @param len Length of the Strings to generate
 * @return Generated Strings in a std::vector
 ******************************************************************************/
std::vector<util::String> genDifferentStrings(const size_t n, const size_t len) {
  const auto alph = util::to_String<util::Symbol>("abc");
  assert(n <= alph.size());
  std::vector<util::String> v;
  for (size_t i = 0; i < n; i++) {
    auto s = util::String(len, alph[i]);
    v.push_back(s);
  }
  return v;
}

/*******************************************************************************
 * @brief Generates hard coded test parameters
 * @details Currently it creates 5 parameters for a problem with
 *  - zero String
 *  - one String
 *  - two Strings that do not share any Symbol
 *  - two Strings that both identical
 *  - three Strings that all identical
 * @param t ConstraintType to set in the problem
 * @param len Length of each String in the problem
 * @param strict Whether a LLCS can be greater than zero for multiple strings
 * @return std::vector<TestParam>
 ******************************************************************************/
std::vector<AlgoParam> genHardCodedParams(
    const ConstraintType t,
    const size_t len,
    const bool strict
) {
  std::vector<AlgoParam> res;
  res.reserve(4);

  auto vec1 = genEqualStrings(0, len);
  res.emplace_back(
      "ZeroStrings", // Problem name
      vec1, // Strings
      util::ParamGenerator::genRelaxedToStdMap(t, vec1),
      new algorithms::solutions::EmptySolution()
  );

  auto vec2 = genEqualStrings(1, len);
  if (strict) {
    res.emplace_back(
        "OneString",
        vec2,
        util::ParamGenerator::genRelaxedToStdMap(t, vec2),
        new algorithms::solutions::EmptySolution());
  } else {
    res.emplace_back(
        "OneString",
        vec2,
        util::ParamGenerator::genRelaxedToStdMap(t, vec2),
        new algorithms::solutions::UnsignedSolution(vec2[0].size())
    );
  }

  auto vec3 = genDifferentStrings(2, len);
  res.emplace_back(
      "TwoDifferentStrings",
      vec3,
      util::ParamGenerator::genRelaxedToStdMap(t, vec3),
      new algorithms::solutions::EmptySolution()
  );

  auto vec4 = genEqualStrings(2, len);
  res.emplace_back(
      "TwoIdenticalStrings",
      vec4,
      util::ParamGenerator::genRelaxedToStdMap(t, vec4),
      new algorithms::solutions::UnsignedSolution(vec4[0].size())
  );

  auto vec5 = genEqualStrings(3, len);
  if (strict) {
    res.emplace_back(
        "ThreeStrings",
        vec5,
        util::ParamGenerator::genRelaxedToStdMap(t, vec5),
        new algorithms::solutions::EmptySolution()
    );
  } else {
    auto sizes = vec5 | std::views::transform(
        [](const util::String &s) { return s.size(); }
    );
    size_t min = *std::ranges::min_element(sizes);
    res.emplace_back(
        "ThreeStrings",
        vec5,
        util::ParamGenerator::genRelaxedToStdMap(t, vec5),
        new algorithms::solutions::UnsignedSolution(min)
    );
  }
  return res;
}

}// namespace lcs_solver::testing::llcs