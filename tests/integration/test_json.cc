#include "constraints/ConstraintMap.h"
#include "constraints/local/Constraint_MC.h"
#include "constraints/local/Constraint_MC_INC.h"
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

namespace {

using lcs_solver::constraints::ConstraintType;
using lcs_solver::constraints::ConstraintMap;
using lcs_solver::constraints::local::Constraint_MC;
using lcs_solver::constraints::local::Constraint_MC_INC;
using GapVector = std::vector<std::pair<size_t, size_t>>;

TEST(JsonIntegrationTest, RoundTrip) {
  GapVector gap_vector;
  gap_vector.emplace_back(0, 0);
  gap_vector.emplace_back(1, 1);

  ConstraintMap map;
  map.emplace(ConstraintType::MC, std::make_shared<Constraint_MC>(gap_vector));
  map.emplace(ConstraintType::MC_INC, std::make_shared<Constraint_MC_INC>(gap_vector));

  Constraint_MC original(gap_vector);
  //
  //   nlohmann::json j;
  //   j["constraint"] = original;
  //   const auto deserialized = j["constraint"].get<Constraint_MC>();


  // for (auto & [type, constraint] : map) {
  //   nlohmann::json j;
  //   j["constraint"] = constraint.get();
  //   const auto deserialized = j["constraint"].get<Constraint_MC>();
  //   EXPECT_EQ(*constraint, deserialized);
  // }
}

}// namespace
