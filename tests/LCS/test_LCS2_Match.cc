/*******************************************************************************
 * @file test_LCS2_Match.cc
 * @author Felix Steinkopp
 * @version 1.2
 * @brief GTests to check LCS2_RT against LCS2_STD_S (using random inputs)
 ******************************************************************************/

#include "algorithms/LCS/LCS2_RT.h"
#include "algorithms/LCS/LCS2_STD_H.h"
#include "algorithms/LCS/LCS2_STD_S.h"
#include "gmock/gmock-more-matchers.h"
#include "gtest/gtest.h"
#include "util/AlgoTestParam.h"
#include "util/ParamGenerator.h"

namespace {
using ::lcs_solver::constraints::ConstraintType;
using ::lcs_solver::util::AlgoParam;
using ::lcs_solver::util::ParamGenerator;

using ::lcs_solver::algorithms::AlgoType;
using ::lcs_solver::algorithms::lcs::LCS2_RT;
using ::lcs_solver::algorithms::lcs::LCS2_STD_H;
using ::lcs_solver::algorithms::lcs::LCS2_STD_S;

using ::testing::Eq;
using ::testing::PrintToString;
using ::testing::ValuesIn;

//==== Initialize Test Suites ==================================================
class LCS2_Match_Test : public testing::TestWithParam<AlgoParam> {};

INSTANTIATE_TEST_SUITE_P(
    LCS2_Reconstructor_Match, LCS2_Match_Test,
    ValuesIn(ParamGenerator::genWithoutSol(
        ConstraintType::Empty,
        std::random_device{}(), // 1542007748
        5,
        {// Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}},
        {0, 3},// bounds for gap lengths
        {'a', 'b', 'c', 'd', 'e'})),
    [](const testing::TestParamInfo<LCS2_Match_Test::ParamType> &info) {
      return info.param.name + "_NoSol";
    });

//=== Definition of Value-Parameterized Tests ==================================
TEST_P(LCS2_Match_Test, STD_Stack_Tree_Eq) {
  //Read out test parameter
  const auto name = GetParam().name;
  const auto &spv = GetParam().pointers;
  // const auto N = spv.size();
  auto map = GetParam().map;

  // Run algorithm
  const auto algo_with_stack = std::make_unique<LCS2_STD_S>(spv,map);
  const auto algo_with_tree = std::make_unique<LCS2_RT>(
    spv,
    map,
    AlgoType::LLCS2_STD_FL
  );
  const auto result_with_stack = algo_with_stack->query();
  const auto result_with_tree = algo_with_tree->query();

  // Compare results
  EXPECT_EQ(*result_with_stack == *result_with_tree, true)
            << "\nDebug Info (" << name << " )\n"
            << "spv[0]:"  << lcs_solver::util::to_string(*spv[0]) << "\n"
            << "spv[1]:"  << lcs_solver::util::to_string(*spv[1]) << "\n"
            << "algoWithStack: \n" << algo_with_stack->DebugString() << "\n"
            << "resultWithStack: \n" << result_with_stack->DebugString() << "\n"
            << "============================================================\n"
            << "algoWithTree: \n" << algo_with_tree->DebugString() << "\n"
            << "resultWithTree: \n" << result_with_tree->DebugString() << "\n"
            << "test: \n" << (*result_with_stack == *result_with_tree) << "\n";
}

TEST_P(LCS2_Match_Test, STD_Stack_Hirsch_Eq) {
  //Read out test parameter
  const auto name = GetParam().name;
  const auto &spv = GetParam().pointers;
  // const auto N = spv.size();
  auto map = GetParam().map;

  // Run algorithm
  const auto algo_with_stack = std::make_unique<LCS2_STD_S>(spv,map);
  const auto algo_hirschberg = std::make_unique<LCS2_STD_H>(spv,map);
  const auto result_with_stack = algo_with_stack->query();
  const auto result_with_hirsch = algo_hirschberg->query();

  // Compare results
  EXPECT_EQ(*result_with_stack == *result_with_hirsch, true)
            << "\nDebug Info (" << name << " )\n"
            << "spv[0]:"  << lcs_solver::util::to_string(*spv[0]) << "\n"
            << "spv[1]:"  << lcs_solver::util::to_string(*spv[1]) << "\n"
            << "algo_stack: \n" << algo_with_stack->DebugString() << "\n"
            << "resul_stack: \n" << result_with_stack->DebugString() << "\n"
            << "============================================================\n"
            // << "algo_hirsch: \n" << algo_hirschberg->DebugString() << "\n"
            << "result_hirsch: \n" << result_with_hirsch->DebugString() << "\n";
}

}// namespace