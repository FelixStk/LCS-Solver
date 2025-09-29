/*******************************************************************************
 * @file test_STD_FL.cc
 * @author Felix Steinkopp
 * @version 1.3
 * @brief GTests for the correctness of the standard llcs-folklore algorithm
 ******************************************************************************/

#include "algorithms/LLCS/LLCS_STD_FL.h"
#include "hardCoded_LLCS.h"

#include <random>

#include "util/ParamGenerator.h"
#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace lcs_solver::testing::llcs {

using ::lcs_solver::algorithms::llcs::LLCS_STD_FL;
using ::lcs_solver::constraints::ConstraintType;
using ::lcs_solver::util::AlgoParam;
using ::lcs_solver::util::ParamGenerator;

using ::testing::Eq;
using ::testing::ValuesIn;
using ::testing::PrintToString;

//==== Initialize Test Suites ==================================================
class LLCS_STD_FL_Test : public ::testing::TestWithParam<AlgoParam> {};

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters, LLCS_STD_FL_Test,
    ValuesIn(genHardCodedParams(ConstraintType::Empty, 8, false)),
    [](const ::testing::TestParamInfo<LLCS_STD_FL_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    rndTwoStrParam, LLCS_STD_FL_Test,
    ValuesIn(ParamGenerator::genWithUniqSol(
        ConstraintType::Empty,
        std::random_device{}(), // Used to initialize a random number generators
        6, // Number of AlgoParam to generate
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {0, 10}, // Bound on llcs length
        {'x', 'a', 'b'} // c[0] common symbol, symbols in the tail are unique
    )),
    [](const ::testing::TestParamInfo<LLCS_STD_FL_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    rndThreeStrParam, LLCS_STD_FL_Test,
    ValuesIn(ParamGenerator::genWithUniqSol(
        ConstraintType::Empty,
        std::random_device{}(),
        5,
        {{5, 6}, {5, 6}, {5, 10}},
        {0, 4},
        {'x', 'a', 'b', 'c'}
    )),
    [](const ::testing::TestParamInfo<LLCS_STD_FL_Test::ParamType> &info) {
      return info.param.name;
    });


//=== Definition of Value-Parameterized Tests ==================================

TEST_P(LLCS_STD_FL_Test, SolutionEq) {
  auto &[name, spv, map, sol] = GetParam();
  const auto N = spv.size();
  auto algo = LLCS_STD_FL(spv, map);
  auto solution = algo.query();
  auto expected = GetParam().sol;

  EXPECT_THAT(solution->empty(), Eq(expected->empty()));
  EXPECT_EQ(*solution, *expected)
            << (N > 0 ? "s[0]: " + util::to_string(*spv[0]) + "\n" : "")
            << (N > 1 ? "s[1]: " + util::to_string(*spv[1]) + "\n" : "")
            << (N > 2 ? "s[2]: " + util::to_string(*spv[2]) + "\n" : "")
            << (N > 3 ? "s[3]: " + util::to_string(*spv[3]) + "\n" : "")
            << "expected: " << expected->DebugString() << "\n"
            << "solution: " << solution->DebugString() << "\n";
}

} // end of namespace