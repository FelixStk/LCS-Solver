#ifndef LCS_SOLVER_CONSTRAINTS_CONSTRAINTMAP_H_
#define LCS_SOLVER_CONSTRAINTS_CONSTRAINTMAP_H_

#include "constraints/BaseConstraint.h"
#include "constraints/ConstraintType.h"
#include "algorithms/AlgoType.h"
#include "algorithms/AlgoCategory.h"
#include <concepts>
#include <map>
#include <nlohmann/json.hpp>

namespace lcs_solver::constraints {

using ConstraintMap = std::map<ConstraintType, std::shared_ptr<BaseConstraint>>;

template<typename T>
concept ConstraintConcept = std::is_base_of_v<BaseConstraint, T> && requires { { T::StaticType() } -> std::same_as<ConstraintType>; };

template<ConstraintConcept T>
T *Get(const ConstraintMap& map) {
  const ConstraintType t = T::StaticType();
  if (!map.contains(t)) return nullptr;
  return dynamic_cast<T *>(map.at(t).get());
}

std::pair<algorithms::AlgoType, algorithms::AlgoType> GetLcsAlgoTypes(
  const ConstraintMap& map,
  bool prefer_max_queues = true,
  bool avoid_segment_tree = true,
  algorithms::AlgoCategory category = algorithms::AlgoCategory::LCS2
);

void to_json(nlohmann::json &j, const BaseConstraint* p);
void to_json(nlohmann::json &j, const std::shared_ptr<BaseConstraint>& p);
void from_json(nlohmann::json &j, std::shared_ptr<BaseConstraint>& p);

void from_json(nlohmann::json &j, ConstraintMap& map);

}// namespace lcs_solver::constraints
#endif//LCS_SOLVER_CONSTRAINTS_CONSTRAINTMAP_H_
