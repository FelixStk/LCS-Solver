/*******************************************************************************
 * @file test_MC.cc
 * @author Felix Steinkopp
 * @version 1.3
 * @brief GTests for the correctness of the standard llcs-folklore algorithm
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_MC.h"
#include "hardCoded_LLCS.h"

#include <random>

#include "util/ParamGenerator.h"
#include "gmock/gmock-more-matchers.h"
#include "gtest/gtest.h"

namespace lcs_solver::testing::llcs {

using algorithms::llcs::LLCS2_MC;
using constraints::ConstraintType;
using util::AlgoParam;
using util::ParamGenerator;

using ::testing::Eq;
using ::testing::PrintToString;
using ::testing::ValuesIn;

//==== Initialize Test Suites ==================================================
class LLCS2_MC_Test : public ::testing::TestWithParam<AlgoParam> {};

INSTANTIATE_TEST_SUITE_P(
    hardCodedParameters,
    LLCS2_MC_Test,
    ValuesIn(genHardCodedParams(ConstraintType::MC)),
    [](const ::testing::TestParamInfo<LLCS2_MC_Test::ParamType> &info) {
      return info.param.name;
    }
);
INSTANTIATE_TEST_SUITE_P(
    UniqueSolParameters,
    LLCS2_MC_Test,
    ValuesIn(ParamGenerator::genWithUniqSol(
        ConstraintType::MC,
        std::random_device{}(), // Seed
        6, // Number of AlgoParam to generate
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}},
        {0, 10}, // Bound on llcs length
        {'x', 'a', 'b'} // c[0] common symbol, symbols in the tail are unique
    )),
    [](const ::testing::TestParamInfo<LLCS2_MC_Test::ParamType> &info) {
      return info.param.name;
    }
);
INSTANTIATE_TEST_SUITE_P(
    RelaxedParameters,
    LLCS2_MC_Test,
    ValuesIn(ParamGenerator::genWithStdSol(
        ConstraintType::MC,
        std::random_device{}(),
        6,
        {// Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}},
        {'a', 'b', 'c', 'd', 'e'})
    ),
    [](const ::testing::TestParamInfo<LLCS2_MC_Test::ParamType> &info) {
      return info.param.name;
    }
);

//=== Definition of Value-Parameterized Tests ==================================
TEST_P(LLCS2_MC_Test, SolutionEq) {
  auto &[name, spv, map, sol] = GetParam();
  const auto algo = std::make_unique<LLCS2_MC>(spv, map);
  const auto solution = algo->query();
  const auto expected = GetParam().sol;

  EXPECT_THAT(solution->empty(), Eq(expected->empty()));
  EXPECT_EQ(*solution, *expected)
            << "test name: " << name << "\n"
            << algo->DebugString()
            << "expected: " << expected->DebugString() << "\n"
            << "solution: " << solution->DebugString() << "\n";
}

}// namespace lcs_solver::testing::llcs