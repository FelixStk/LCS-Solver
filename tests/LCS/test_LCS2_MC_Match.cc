/*******************************************************************************
 * @file test_LCS_MC_Match.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief GTests to check LCS2_RT reconstructions against each other
 ******************************************************************************/

#include "algorithms/LCS/LCS2_RT.h"
#include "algorithms/solutions/BaseCollector.h"
#include "algorithms/solutions/Vector3DSolution.h"
#include "constraints/local/Constraint_MC_1C.h"
#include "gmock/gmock-more-matchers.h"
#include "gtest/gtest.h"
#include "util/AlgoTestParam.h"
#include "util/ParamGenerator.h"

namespace {
using lcs_solver::algorithms::AlgoType;
using lcs_solver::algorithms::BaseAlgorithm;
using lcs_solver::algorithms::BaseSolution;
using lcs_solver::algorithms::solutions::BaseCollector;
using lcs_solver::algorithms::lcs::LCS2_RT;
using lcs_solver::algorithms::solutions::Vector3DSolution;
using lcs_solver::constraints::BaseConstraint;
using lcs_solver::constraints::ConstraintMap;
using lcs_solver::constraints::ConstraintType;
using lcs_solver::constraints::local::Constraint_MC;
using lcs_solver::util::AlgoParam;
using lcs_solver::util::uint;
using lcs_solver::util::ParamGenerator;
using testing::Eq;
using testing::PrintToString;
using testing::ValuesIn;
using Pair  = std::pair<uint, uint>;

//==== Generate Test Parameters=================================================

using TestParam = std::pair<AlgoType, AlgoParam>;
std::vector<TestParam> GetParams(AlgoType t, unsigned int seed = 0, unsigned int n = 0) {
  if (!seed)
    seed = std::random_device{}();
  if (!n)
    n = 5;
  ConstraintType ct = ConstraintType::Empty;
  if (t == AlgoType::LLCS2_MC_INC || t == AlgoType::LLCS2_MC_INC_E) {
    ct = ConstraintType::MC_INC;
  }
  if (t == AlgoType::LLCS2_MC_1C) {
    ct = ConstraintType::MC_1C;
  }
  if (t == AlgoType::LLCS2_MC_O1_SYNC) {
    ct = ConstraintType::MC_O1C_SYNC;
  }
  std::vector<AlgoParam> algo_params = ParamGenerator::genWithoutSol(
      ct,
      seed,
      n,
      {// Bounds for string length: l[k]={a,b} <=> a <=s [k].size <= b
       {5, 10},
       {5, 10}},
      {0, 3},  // bounds for gap lengths
      {'a', 'b', 'c', 'd', 'e'});

  // Add MC Constraint
  for (auto& param : algo_params) {
    BaseConstraint* ptr = nullptr;
    if (param.map.contains(ConstraintType::MC_INC)) {
      ptr = param.map[ConstraintType::MC_INC].get();
    }
    else if (param.map.contains(ConstraintType::MC_1C)) {
      ptr = param.map[ConstraintType::MC_1C].get();
    }
    else if (param.map.contains(ConstraintType::MC_O1C_SYNC)) {
      ptr = param.map[ConstraintType::MC_O1C_SYNC].get();
    }
    auto constraint_ptr = dynamic_cast<Constraint_MC*>(ptr);
    std::vector<Pair> gaps = constraint_ptr->GetGapVector();
    param.map[ConstraintType::MC] = std::make_shared<Constraint_MC>(gaps);
  }

  std::vector<TestParam> result;
  for (auto& param : algo_params) {
    result.emplace_back(t, param);
  }
  return result;
}

//==== Initialize Test Suites ==================================================
class LCS2_MC_Equal_Test : public testing::TestWithParam<TestParam> {};
class LCS2_MC_Subset_Test : public testing::TestWithParam<TestParam> {};

INSTANTIATE_TEST_SUITE_P(
    MC_INC,
    LCS2_MC_Equal_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_MC_INC)),
    [](const testing::TestParamInfo<LCS2_MC_Equal_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    MC_INC_E,
    LCS2_MC_Equal_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_MC_INC_E)),
    [](const testing::TestParamInfo<LCS2_MC_Equal_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    MC_C1,
    LCS2_MC_Equal_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_MC_1C)),
    [](const testing::TestParamInfo<LCS2_MC_Equal_Test::ParamType>& info) {
      return info.param.second.name;
    });

INSTANTIATE_TEST_SUITE_P(
    MC_O1C_SYNC,
    LCS2_MC_Subset_Test,
    ValuesIn(GetParams(AlgoType::LLCS2_MC_O1_SYNC)),
    [](const testing::TestParamInfo<LCS2_MC_Subset_Test::ParamType>& info) {
      return info.param.second.name;
    });

//=== Definition of Value-Parameterized Tests ==================================
TEST_P(LCS2_MC_Equal_Test, SolutionEq) {
  // Read out test parameter
  const AlgoType& algo_type = GetParam().first;
  const std::string& name = GetParam().second.name;
  const BaseAlgorithm::StrPtrVector& spv = GetParam().second.pointers;
  ConstraintMap map = GetParam().second.map;

  // Run algorithm
  const auto algo1 = std::make_unique<LCS2_RT>(spv, map, AlgoType::LLCS2_MC);
  const auto algo2 = std::make_unique<LCS2_RT>(spv, map, algo_type);
  auto solution1 = algo1->query();
  auto solution2 = algo2->query();
  const auto sol_mc_ptr = algo1->query();
  const auto sol_t_ptr = algo2->query();

  // Compare results
  EXPECT_EQ(*sol_mc_ptr == *sol_t_ptr, true)
      << "\nDebug Info (" << name << ")\n"
      << "spv[0]: " << lcs_solver::util::to_string(*spv[0]) << "\n"
      << "spv[1]: " << lcs_solver::util::to_string(*spv[1]) << "\n"
      << "algo2->DebugString(): \n"
      << algo2->DebugString() << "\n"
      << "============================================================\n"
      << "sol_mc_ptr->DebugString(): \n"
      << sol_mc_ptr->DebugString() << "\n\n"
      << "sol_t_ptr->DebugString(): \n"
      << sol_t_ptr->DebugString() << "\n"
      << "gaps->DebugString(): \n "
      << map.at(ConstraintType::MC)->DebugString() ;
}

TEST_P(LCS2_MC_Subset_Test, SubsetEq) {
  // Read out test parameter
  const AlgoType& algo_type = GetParam().first;
  const std::string& name = GetParam().second.name;
  const BaseAlgorithm::StrPtrVector& spv = GetParam().second.pointers;
  ConstraintMap map = GetParam().second.map;

  // Run algorithm
  const auto algo1 = std::make_unique<LCS2_RT>(spv, map, AlgoType::LLCS2_MC);
  const auto algo2 = std::make_unique<LCS2_RT>(spv, map, algo_type);
  const auto sol_mc_ptr = algo1->query();
  const auto sol_t_ptr = algo2->query();

  // Compare results
  const auto raw_mc_ptr = dynamic_cast<BaseCollector*>(sol_mc_ptr.get());
  const auto raw_t_ptr  = dynamic_cast<BaseCollector*>(sol_t_ptr.get());
  EXPECT_EQ(raw_t_ptr->isSubset(*raw_mc_ptr), true)
      << "\nDebug Info (" << name << ")\n"
      << "spv[0]: " << lcs_solver::util::to_string(*spv[0]) << "\n"
      << "spv[1]: " << lcs_solver::util::to_string(*spv[1]) << "\n"
      << "algo2->DebugString(): \n"
      << algo2->DebugString() << "\n"
      << "============================================================\n"
      << "sol_mc_ptr->DebugString(): \n"
      << sol_mc_ptr->DebugString() << "\n\n"
      << "sol_t_ptr->DebugString(): \n"
      << sol_t_ptr->DebugString() << "\n"
      << "gaps->DebugString(): \n "
      << map.at(ConstraintType::MC)->DebugString() ;
}

}  // namespace