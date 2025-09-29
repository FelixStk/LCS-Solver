/*******************************************************************************
 * @file test_MC_INC_Match.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief Compares LLCS_MC_INC1 with LLCS_MC_INC2
 ******************************************************************************/

#include <vector>

#include "algorithms/LLCS/LLCS2_STD_FL.h"
#include "algorithms/LLCS/LLCS2_STD_H.h"
#include "algorithms/solutions/EmptySolution.h"
#include "gmock/gmock-more-matchers.h"
#include "gtest/gtest.h"
#include "util/AlgoTestParam.h"
#include "util/ParamGenerator.h"

namespace lcs_solver::testing::llcs {
using algorithms::llcs::LLCS2_STD_FL;
using algorithms::llcs::LLCS2_STD_H;

using constraints::ConstraintType;
using util::AlgoParam;
using util::ParamGenerator;

using ::testing::Eq;
using ::testing::ValuesIn;
using ::testing::PrintToString;

//==== Initialize Test Suites ==================================================
class LLCS2_STD_Match_Test : public ::testing::TestWithParam<AlgoParam> {};

INSTANTIATE_TEST_SUITE_P(
    LLCS2_STD_Match_Test, LLCS2_STD_Match_Test,
    ValuesIn(ParamGenerator::genWithoutSol(
        ConstraintType::Empty,
        std::random_device{}(), //3827720579, //std::random_device{}(),
        5,
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {0, 3}, // bounds for gap lengths
        {'a', 'b', 'c', 'd', 'e'})),
    [](const ::testing::TestParamInfo<LLCS2_STD_Match_Test::ParamType> &info) {
      return info.param.name + "_NoSol";
    });

//=== Definition of Value-Parameterized Tests ==================================
TEST_P(LLCS2_STD_Match_Test, SolutionEqual) {
  auto &[name, spv, map, sol] = GetParam();
  const auto n = spv.size();

  const auto algo1 = std::make_unique<LLCS2_STD_FL>(spv, map);
  const auto algo2 = std::make_unique<LLCS2_STD_H>(spv, map);
  const auto solution1 = algo1->query();
  const auto solution2 = algo2->query();

  EXPECT_THAT(solution1->empty(), Eq(solution2->empty()));
  EXPECT_EQ(*solution1, *solution2)
            << "Name: " << GetParam().name << "\n"
            << "Lengths: " << spv[0]->size() << " " << spv[1]->size() << "\n"
            << (n > 0 ? "s[0]: " + util::to_string(*spv[0]) + "\n" : "")
            << (n > 1 ? "s[1]: " + util::to_string(*spv[1]) + "\n" : "")
            << "solution_FL: " << solution1->DebugString() << "\n"
            << "solution_H : " << solution2->DebugString() << "\n"
            << "LLCS2_STD_LF:\n" << algo1->DebugString() << "\n"
            << "============================================================\n"
            << "LLCS2_STD_H:\n " << algo2->DebugString() << "\n";
}

} // end of namespace