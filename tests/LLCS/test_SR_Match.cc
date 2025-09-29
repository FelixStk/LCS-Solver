/*******************************************************************************
 * @file test_SR_Match.cc
 * @author Felix Steinkopp
 * @version 1.3
 * @brief GTests if the llcs of algorithms (llcs for different local sigma
 * constrained problems) are consistent with each other:
 * - Compares LLCS_SR_MQ with LLCS_SR_RMQ
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_SR_RMQ.h"
#include "algorithms/LLCS/LLCS2_SR_MQ.h"
#include "hardCoded_LLCS.h"

#include <random>

#include "util/ParamGenerator.h"
#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace lcs_solver::testing::llcs {
using ::lcs_solver::algorithms::llcs::LLCS2_SR_MQ;
using ::lcs_solver::algorithms::llcs::LLCS2_SR_RMQ;
using ::lcs_solver::constraints::ConstraintType;
using ::lcs_solver::util::AlgoParam;
using ::lcs_solver::util::ParamGenerator;

using ::testing::Eq;
using ::testing::ValuesIn;
using ::testing::PrintToString;

//==== Initialize Test Suites ==================================================
class LLCS2_SR_Match_Test : public ::testing::TestWithParam<AlgoParam> {};

INSTANTIATE_TEST_SUITE_P(
    LLCS2_SR_MQ_Match, LLCS2_SR_Match_Test,
    ValuesIn(ParamGenerator::genWithoutSol(
        ConstraintType::SIGMA_R,
        std::random_device{}(),
        6,
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {0, 3}, // bounds for gap lengths
        {'a', 'b', 'c', 'd', 'e'}
    )),
    [](const ::testing::TestParamInfo<LLCS2_SR_Match_Test::ParamType> &info) {
      return info.param.name + "_NoSol";
    });

//=== Definition of Value-Parameterized Tests ==================================
TEST_P(LLCS2_SR_Match_Test, RMQSolEqMQSol) {
  auto &[name, spv, map, sol] = GetParam();
  const auto N = spv.size();
  auto algo1 = std::make_unique<LLCS2_SR_RMQ>(spv, map);
  auto algo2 = std::make_unique<LLCS2_SR_MQ>(spv, map);
  auto solutionRMQ = algo1->query();
  auto solutionMQ = algo2->query();

  EXPECT_THAT(solutionRMQ->empty(), Eq(solutionMQ->empty()));
  EXPECT_EQ(*solutionRMQ, *solutionMQ)
            << "Name: " << name << "\n"
            << "Right: " << PrintToString(algo2->right) << "\n"
            << "Lengths: " << spv[0]->size() << " " << spv[1]->size() << "\n"
            << (N > 0 ? "s[0]: " + util::to_string(*spv[0]) + "\n" : "")
            << (N > 1 ? "s[1]: " + util::to_string(*spv[1]) + "\n" : "")
            << (N > 2 ? "s[2]: " + util::to_string(*spv[2]) + "\n" : "")
            << (N > 3 ? "s[3]: " + util::to_string(*spv[3]) + "\n" : "")
            << "solutionRMQ: " << solutionRMQ->DebugString() << "\n"
            << "solutionMQ : " << solutionMQ->DebugString() << "\n"
            << "LLCS2_SR_RMQ: " << algo1->DebugString() << "\n"
            << "============================================================\n"
            << "LLCS2_SR_MQ: " << algo2->DebugString() << "\n";
}

} // end of namespace