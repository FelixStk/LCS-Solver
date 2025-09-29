/*******************************************************************************
 * @file test_SR_RMQ.cc
 * @author Felix Steinkopp
 * @version 1.3
 * @brief GTests for the correctness of the LLCS2_SR_MQ algorithm
 ******************************************************************************/

#include "algorithms/LLCS/LLCS2_SR_RMQ.h"
#include "algorithms/LLCS/LLCS2_SR_MQ.h"
#include "hardCoded_LLCS.h"

#include <random>

#include "util/ParamGenerator.h"
#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace lcs_solver::testing::llcs {
using ::lcs_solver::algorithms::llcs::LLCS2_SR_RMQ;
using ::lcs_solver::constraints::ConstraintType;
using ::lcs_solver::util::AlgoParam;
using ::lcs_solver::util::ParamGenerator;

using ::testing::Eq;
using ::testing::ValuesIn;
using ::testing::PrintToString;

//==== Initialize Test Suites ==================================================
class LLCS2_SR_RMQ_Test : public ::testing::TestWithParam<AlgoParam> {};

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters, LLCS2_SR_RMQ_Test,
    ValuesIn(genHardCodedParams(ConstraintType::SIGMA_R)),
    [](const ::testing::TestParamInfo<LLCS2_SR_RMQ_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    UniqueSolParameters, LLCS2_SR_RMQ_Test,
    ValuesIn(ParamGenerator::genWithUniqSol(
        ConstraintType::SIGMA_R,
        std::random_device{}(), // Used to initialize a random number generators
        6, // Number of AlgoParam to generate
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {0, 10}, // Bound on llcs length
        {'x', 'a', 'b'} // c[0] common symbol, symbols in the tail are unique
    )),
    [](const ::testing::TestParamInfo<LLCS2_SR_RMQ_Test::ParamType> &info) {
      return info.param.name;
    });

INSTANTIATE_TEST_SUITE_P(
    StdRelaxedParameters, LLCS2_SR_RMQ_Test,
    ValuesIn(ParamGenerator::genWithStdSol(
        ConstraintType::SIGMA_R,
        std::random_device{}(),
        6,
        { // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
            {5, 10},
            {5, 10}
        },
        {'a', 'b', 'c', 'd', 'e'}
    )),
    [](const ::testing::TestParamInfo<LLCS2_SR_RMQ_Test::ParamType> &info) {
      return info.param.name;
    });

//=== Definition of Value-Parameterized Tests ==================================

TEST_P(LLCS2_SR_RMQ_Test, SolutionEq) {
  auto &[name, spv, map, sol] = GetParam();
  const auto N = spv.size();
  auto algo = std::make_unique<LLCS2_SR_RMQ>(spv, map);
  const auto & s = algo->getStringViewVec();
  auto solution = algo->query();
  auto expected = GetParam().sol;

  EXPECT_THAT(solution->empty(), Eq(expected->empty()));
  EXPECT_EQ(*solution, *expected)
            << "Name: " << GetParam().name << "\n"
            << "Right: " << PrintToString(algo->right) << "\n"
            << "Lengths: " << s[0].size() << " " << s[1].size() << "\n"
                << (N > 0 ? "s[0]: " + util::to_string(*spv[0]) + "\n" : "")
                << (N > 1 ? "s[1]: " + util::to_string(*spv[1]) + "\n" : "")
                << (N > 2 ? "s[2]: " + util::to_string(*spv[2]) + "\n" : "")
                << (N > 3 ? "s[3]: " + util::to_string(*spv[3]) + "\n" : "")
            << "expected: " << expected->DebugString() << "\n"
            << "solution: " << solution->DebugString() << "\n"
            << algo->DebugString() << "\n";
}

} // end of namespace