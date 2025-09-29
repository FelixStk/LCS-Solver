/*******************************************************************************
 * @file test_MC_1C.cc
 * @author Felix Steinkopp
 * @version 1.3
 * @brief GTests for the correctness of the LLCS2_MC_1C algorithm
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_MC_1C.h"
#include "hardCoded_LLCS.h"
#include "util/AlgoTestParam.h"
#include "util/ParamGenerator.h"

#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace lcs_solver::testing::llcs {
using algorithms::llcs::LLCS2_MC_1C;
using constraints::ConstraintType;
using util::AlgoParam;
using util::ParamGenerator;

using ::testing::Eq;
using ::testing::PrintToString;
using ::testing::ValuesIn;

//==== Initialize Test Suites ==================================================
class LLCS2_MC_1C_Test : public ::testing::TestWithParam<AlgoParam> {};

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters,
    LLCS2_MC_1C_Test,
    ValuesIn(genHardCodedParams(ConstraintType::MC_1C)),
    [](const ::testing::TestParamInfo<LLCS2_MC_1C_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    UniqueSolParameters,
    LLCS2_MC_1C_Test,
    ValuesIn(ParamGenerator::genWithUniqSol(// strings have a unique lcs in each string
        ConstraintType::MC_1C,
        std::random_device{}(), // Seed
        6, // Number of AlgoParam to generate
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {0, 10}, // Bound on llcs length
        {'x', 'a', 'b'} // c[0] common symbol, symbols in the tail are unique
    )),
    [](const ::testing::TestParamInfo<LLCS2_MC_1C_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    StdRelaxedParameters,
    LLCS2_MC_1C_Test,
    ValuesIn(ParamGenerator::genWithMCSol(// generate rnd strings, relax constraint to folklore problem
        ConstraintType::MC_1C,
        std::random_device{}(),
        6,
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {0, 5}, // Super-interval from which gaps are drawn
        {'a', 'b', 'c', 'd', 'e'}
    )),
    [](const ::testing::TestParamInfo<LLCS2_MC_1C_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    MCRelaxedParameters,
    LLCS2_MC_1C_Test,
    ValuesIn(ParamGenerator::genWithMCSol( // genWithMCSol: generate rnd strings, relax constraint to mc gaps problem
        ConstraintType::MC_1C,
        std::random_device{}(),
        6,
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {0, 5}, // Super-interval from which gaps are drawn
        {'a', 'b', 'c', 'd', 'e'}
    )),
    [](const ::testing::TestParamInfo<LLCS2_MC_1C_Test::ParamType> &info) {
      return info.param.name;
    });

//=== Definition of Value-Parameterized Tests ==================================
TEST_P(LLCS2_MC_1C_Test, SolutionEq) {
  auto &[name, spv, map, sol] = GetParam();
  auto algo = std::make_unique<LLCS2_MC_1C>(spv, map);
  auto solution = algo->query();
  auto expected = GetParam().sol;

  EXPECT_THAT(solution->empty(), Eq(expected->empty()));
  EXPECT_EQ(*solution, *expected)
            << "test name: " << name << "\n"
            << "expected: " << expected->DebugString() << "\n"
            << "solution: " << solution->DebugString() << "\n"
            << algo->DebugString() << "\n";
}

}// namespace lcs_solver::testing::llcs