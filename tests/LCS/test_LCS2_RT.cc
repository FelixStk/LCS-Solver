/*******************************************************************************
 * @file test_LCS_RT2.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief GTests to check LCS2_RT against Vector3DSolutions generated with
 *        ParamGenerator::genWithUniqSol
 * @see ParamGenerator::genWithUniqSol
 ******************************************************************************/

#include "algorithms/LCS/LCS2_RT.h"
#include "algorithms/solutions/Vector3DSolution.h"
#include "gmock/gmock-more-matchers.h"
#include "gtest/gtest.h"
#include "util/AlgoTestParam.h"
#include "util/ParamGenerator.h"

namespace {
using lcs_solver::algorithms::AlgoType;
using lcs_solver::algorithms::BaseAlgorithm;
using lcs_solver::algorithms::BaseSolution;
using lcs_solver::algorithms::lcs::LCS2_RT;
using lcs_solver::algorithms::solutions::Vector3DSolution;
using lcs_solver::constraints::ConstraintMap;
using lcs_solver::constraints::ConstraintType;
using lcs_solver::util::AlgoParam;
using lcs_solver::util::ParamGenerator;
using testing::Eq;
using testing::PrintToString;
using testing::ValuesIn;

//==== Generate Test Parameters=================================================

using TestParam = std::pair<AlgoType, AlgoParam>;
std::vector<TestParam> GetParams(AlgoType t) {
  static std::map<AlgoType, std::vector<AlgoParam>> map;
  static constexpr int n = 5;
  // static constexpr unsigned int seed = 2630895230;
  static const auto seed = std::random_device{}();
  static constexpr std::array types =
      std::to_array<std::pair<AlgoType, ConstraintType>>({
          std::make_pair(AlgoType::LLCS2_STD_FL, ConstraintType::Empty),
          std::make_pair(AlgoType::LLCS2_MC, ConstraintType::MC),
          std::make_pair(AlgoType::LLCS2_MC_1C, ConstraintType::MC_1C),
          std::make_pair(AlgoType::LLCS2_MC_INC, ConstraintType::MC_INC),
          std::make_pair(AlgoType::LLCS2_MC_INC_E, ConstraintType::MC_INC),
          std::make_pair(AlgoType::LLCS2_MC_O1_SYNC,
                         ConstraintType::MC_O1C_SYNC),
          std::make_pair(AlgoType::LLCS2_SR_MQ, ConstraintType::SIGMA_R),
          std::make_pair(AlgoType::LLCS2_SR_RMQ, ConstraintType::SIGMA_R),
          std::make_pair(AlgoType::LLCS2_SL_MQ, ConstraintType::SIGMA_L),
          std::make_pair(AlgoType::LLCS2_SL_RMQ, ConstraintType::SIGMA_L),
          std::make_pair(AlgoType::LLCS2_SA_MQ, ConstraintType::SIGMA),
          std::make_pair(AlgoType::LLCS2_SA_RMQ, ConstraintType::SIGMA),
      });
  if (map.empty()) {
    for (const auto& [algo_type, constraint_type] : types) {
      map.emplace(
          algo_type,
          ParamGenerator::genWithUniqSol(
              constraint_type,
              seed,
              n,  // Number of parameters to generate
              {
                  // Bounds for string length: l[k]={a,b} <=> a <=s [k].size <=
                  // b
                  {5, 10},
                  {5, 10},
              },
              {0, 3},                     // bounds for gap lengths
              {'a', 'b', 'c', 'd', 'e'},  // 'a' == match between strings
              false  // switch solution type to Vector3DSolution
              ));
    }
  }
  std::vector<TestParam> result;
  for (auto& param : map.at(t)) {
    result.emplace_back(t, param);
  }
  return result;
}

//==== Initialize Test Suites ==================================================
class LCS2_RT2_Test : public testing::TestWithParam<TestParam> {};

INSTANTIATE_TEST_SUITE_P(
    Uniq_STD_FL,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_STD_FL)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_MC,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_MC)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_MC_1C,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_MC_1C)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_MC_INC,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_MC_INC)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_MC_INC_E,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_MC_INC_E)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_MC_O1_SYNC,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_MC_O1_SYNC)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_SR_MQ,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_SR_MQ)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_SR_RMQ,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_SR_RMQ)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_SL_MQ,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_SL_MQ)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_SL_RMQ,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_SL_RMQ)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_SA_MQ,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_SA_MQ)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    Uniq_SA_RMQ,
    LCS2_RT2_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_SA_RMQ)),
    [](const testing::TestParamInfo<LCS2_RT2_Test::ParamType>& info) {
      return info.param.second.name;
    });

//=== Definition of Value-Parameterized Tests ==================================
TEST_P(LCS2_RT2_Test, SolutionEq) {
  // Read out test parameter
  const AlgoType& algo_type = GetParam().first;
  const std::string& name = GetParam().second.name;
  const BaseAlgorithm::StrPtrVector& spv = GetParam().second.pointers;
  const BaseSolution* expected_ptr = GetParam().second.sol;
  const auto* test_ptr = dynamic_cast<const Vector3DSolution*>(expected_ptr);
  ConstraintMap map = GetParam().second.map;

  // Run algorithm
  const auto algo = std::make_unique<LCS2_RT>(spv, map, algo_type);
  const auto sol_ptr = algo->query();

  // Compare results
  EXPECT_EQ(*sol_ptr == *test_ptr, true)
      << "\nDebug Info (" << name << ")\n"
      << "spv[0]: " << lcs_solver::util::to_string(*spv[0]) << "\n"
      << "spv[1]: " << lcs_solver::util::to_string(*spv[1]) << "\n"
      << "algo->DebugString(): \n"
      << algo->DebugString() << "\n"
      << "============================================================\n"
      << "sol_ptr->DebugString(): \n"
      << sol_ptr->DebugString() << "\n\n"
      << "expected_ptr->DebugString(): \n"
      << expected_ptr->DebugString() << "\n";
}

}  // namespace