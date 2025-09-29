/*******************************************************************************
 * @file test_MC_INC_Match.cc
 * @author Felix Steinkopp
 * @version 1.3
 * @brief Compares LLCS_MC_INC1 with LLCS_MC_INC2
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_MC_INC.h"
#include "algorithms/LLCS/LLCS2_MC_INC_E.h"
#include "algorithms/solutions/EmptySolution.h"
#include "algorithms/solutions/UnsignedSolution.h"
#include "util/AlgoTestParam.h"
#include "util/ParamGenerator.h"

#include<vector>

#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace lcs_solver::testing::llcs {
using ::lcs_solver::algorithms::llcs::LLCS2_MC_INC;
using ::lcs_solver::algorithms::llcs::LLCS2_MC_INC_E;

using ::lcs_solver::constraints::ConstraintType;
using ::lcs_solver::util::AlgoParam;
using ::lcs_solver::util::ParamGenerator;

using ::testing::Eq;
using ::testing::ValuesIn;
using ::testing::PrintToString;

//==== Initialize Test Suites ==================================================
class LLCS2_MC_INC_Match_Test : public ::testing::TestWithParam<AlgoParam> {};

INSTANTIATE_TEST_SUITE_P(
    LLCS2_MC_INC_Match, LLCS2_MC_INC_Match_Test,
    ValuesIn(ParamGenerator::genWithoutSol(
        ConstraintType::MC_INC,
        std::random_device{}(),
        1,
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {0, 3}, // bounds for gap lengths
        {'a', 'b', 'c', 'd', 'e'})),
    [](const ::testing::TestParamInfo<LLCS2_MC_INC_Match_Test::ParamType> &info) {
      return info.param.name + "_NoSol";
    });

//=== Definition of Value-Parameterized Tests ==================================
TEST_P(LLCS2_MC_INC_Match_Test, SolutionEqual) {
  auto &[name, spv, map, sol] = GetParam();
  const auto N = spv.size();

  auto algo1 = std::make_unique<LLCS2_MC_INC>(spv, map);
  auto algo2 = std::make_unique<LLCS2_MC_INC_E>(spv, map);
  auto solution1 = algo1->query();
  auto solution2 = algo2->query();

  EXPECT_THAT(solution1->empty(), Eq(solution2->empty()));
  EXPECT_EQ(*solution1, *solution2)
            << "Name: " << GetParam().name << "\n"
            << "Gaps: " << PrintToString(algo2->getGaps()) << "\n"
            << "Lengths: " << spv[0]->size() << " " << spv[1]->size() << "\n"
            << (N > 0 ? "s[0]: " + util::to_string(*spv[0]) + "\n" : "")
            << (N > 1 ? "s[1]: " + util::to_string(*spv[1]) + "\n" : "")
            << (N > 2 ? "s[2]: " + util::to_string(*spv[2]) + "\n" : "")
            << (N > 3 ? "s[3]: " + util::to_string(*spv[3]) + "\n" : "")
            << "solutionSegT: " << solution1->DebugString() << "\n"
            << "solutionEff : " << solution2->DebugString() << "\n"
            << "LLCS2_MC_INC: " << algo1->DebugString() << "\n"
            << "============================================================\n"
            << "LLCS2_MC_INC_E: " << algo2->DebugString() << "\n";
}

} // end of namespace