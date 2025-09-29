/*******************************************************************************
 * @file test_MC_INC.cc
 * @author Felix Steinkopp
 * @version 1.3
 * @brief GTests for the correctness of the LLCS2_MC_INC algorithm
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_MC_INC.h"
#include "algorithms/LLCS/LLCS2_MC_INC_E.h"
#include "hardCoded_LLCS.h"

#include <algorithm>
#include <memory>
#include <random>
#include <string>
#include <utility>

#include "util/AlgoTestParam.h"
#include "util/ParamGenerator.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace lcs_solver::testing::llcs {

using ::lcs_solver::algorithms::llcs::LLCS2_MC_INC;
using ::lcs_solver::algorithms::llcs::LLCS2_MC_INC_E;

using ::lcs_solver::algorithms::BaseSolution;
using ::lcs_solver::constraints::ConstraintType;
using ::lcs_solver::util::AlgoParam;
using ::lcs_solver::util::String;
using ::lcs_solver::util::ParamGenerator;

using ::testing::Eq;
using ::testing::ValuesIn;
using ::testing::PrintToString;

//==== Initialize Test Suites ==================================================
class LLCS2_MC_INC_Test : public ::testing::TestWithParam<AlgoParam> {};

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters, LLCS2_MC_INC_Test,
    ValuesIn(genHardCodedParams(ConstraintType::MC_INC)),
    [](const ::testing::TestParamInfo<LLCS2_MC_INC_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    UniqueSolParameters, LLCS2_MC_INC_Test,
    ValuesIn(ParamGenerator::genWithUniqSol(
        ConstraintType::MC_INC,
        std::random_device{}(), // Used to initialize a random number generators
        6, // Number of AlgoParam to generate
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 5},
            {5, 5}
        },
        {0, 5}, // Bound on llcs length
        {'x', 'a', 'b'} // c[0] common symbol, symbols in the tail are unique
    )),
    [](const ::testing::TestParamInfo<LLCS2_MC_INC_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    StdRelaxedParameters, LLCS2_MC_INC_Test,
    ValuesIn( ParamGenerator::genWithStdSol(
        ConstraintType::MC_INC,
        std::random_device{}(),
        6,
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {'a', 'b', 'c', 'd', 'e'}
    )),
    [](const ::testing::TestParamInfo<LLCS2_MC_INC_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    MCRelaxedParameters, LLCS2_MC_INC_Test,
    ValuesIn(ParamGenerator::genWithMCSol(
        ConstraintType::MC_INC,
        std::random_device{}(),
        6,
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {0, 5}, // Super-interval from which gaps are drawn
        {'a', 'b', 'c', 'd', 'e'}
    )),
    [](const ::testing::TestParamInfo<LLCS2_MC_INC_Test::ParamType> &info) {
      return info.param.name;
    });

//=== Definition of Value-Parameterized Tests ==================================

TEST_P(LLCS2_MC_INC_Test, SolutionEq_SegTree) {
  auto &[name, spv, map, sol] = GetParam();
  auto algo = std::make_unique<LLCS2_MC_INC>(spv, map);
  auto solution = algo->query();
  auto *expected = GetParam().sol;

  EXPECT_THAT(solution->empty(), Eq(expected->empty()));
  EXPECT_EQ(*solution, *expected)
            << "test name: " << name << "\n"
            << "expected: " << expected->DebugString() << "\n"
            << "solution: " << solution->DebugString() << "\n"
            << algo->DebugString() << "\n";
}

TEST_P(LLCS2_MC_INC_Test, SolutionEq_CRBuffer) {
  auto &[name, spv, map, sol] = GetParam();
  auto algo = std::make_unique<LLCS2_MC_INC_E>(spv, map);
  auto solution = algo->query();
  auto *expected = GetParam().sol;
  auto gaps = algo->getGaps();

  EXPECT_THAT(solution->empty(), Eq(expected->empty()));
  EXPECT_EQ(*solution, *expected)
            << "test name: " << name << "\n"
            << "expected: " << expected->DebugString() << "\n"
            << "solution: " << solution->DebugString() << "\n"
            << algo->DebugString() << "\n"
            << PrintToString(gaps) << "\n";
}

} // end of namespace