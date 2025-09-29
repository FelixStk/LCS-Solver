/******************************************************************************
 * @file ConstraintMap.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Implementation of to_json and from_json for ConstraintMap
 *****************************************************************************/

#include "constraints/ConstraintMap.h"

#include "constraints/ConstraintType.h"
#include "constraints/global/Constraint_BR.h"
#include "constraints/global/Constraint_LLCS_Known.h"
#include "constraints/input/Constraint_Const_Sigma.h"
#include "constraints/input/Constraint_LCS_2.h"
#include "constraints/input/Constraint_LCS_K.h"
#include "constraints/local/Constraint_MC.h"
#include "constraints/local/Constraint_MC_1C.h"
#include "constraints/local/Constraint_MC_INC.h"
#include "constraints/local/Constraint_MC_O1C.h"
#include "constraints/local/Constraint_MC_O1C_SYNC.h"
#include "constraints/local/Constraint_Sigma.h"
#include "constraints/local/Constraint_Sigma_L.h"
#include "constraints/local/Constraint_Sigma_R.h"

namespace lcs_solver::constraints {

/*******************************************************************************
 * @brief Determines the appropriate algorithm types for solving a constrained
 * LCS problem without user interaction.
 *
 * @details Given a set of constraints, this function selects a pair of
 * algorithm types that can be used with `AlgoFactory::Create()` to construct a
 * suitable `BaseAlgorithm` for solving problems with constraints.
 *
 * The returned pair represents:
 * - `main`: the primary algorithm type used to solve the problem.
 * - `child`: the secondary algorithm type used by the main algorithm, if
 *   applicable. For example, `LCS2_RT` may wrap a specific LLCS2 Algorith to
 *   reconstruct an LCS from
 *
 * The choice depends on the constraint types present in the map and the
 * specified preferences. The result is only well-defined if there is at most
 * one local constraint in the map.
 *
 * @param map A map of constraints (`ConstraintType` to `BaseConstraint`).
 * @param prefer_max_queues Whether to use algorithms based on the `MaxQueue2D`
 *        class (default: true)
 * @param avoid_segment_tree Whether to avoid algorithms that require segment
 *        trees (default: true).
 * @param category The category of algorithms to consider (e.g., `LLCS2` or
 *        `LCS2`).
 *
 * @return A pair of `AlgoType`s: `{main, child}`. If the selected algorithm
 * category is a wrapper (like `LCS2_RT`), the `child` type indicates the
 * internal LLCS algorithm to use.
 *
 * @throws std::runtime_error if the algorithm category is unsupported, or if
 * unsupported constraints like `Constraint_LCS_K` with `k != 2` are present.
 ******************************************************************************/
std::pair<algorithms::AlgoType, algorithms::AlgoType> GetLcsAlgoTypes(
    const ConstraintMap &map,
    const bool prefer_max_queues,
    const bool avoid_segment_tree,
    const algorithms::AlgoCategory category) {
  using algorithms::AlgoCategory;
  using algorithms::AlgoType;

  switch (category) {
    case AlgoCategory::Empty:
    case AlgoCategory::LLCS:
    case AlgoCategory::LCS:
    case AlgoCategory::Oracle:
    case AlgoCategory::Online:
    case AlgoCategory::Other:
      throw std::runtime_error("GetLcsAlgoTypes called with unsupported category");
    case AlgoCategory::LLCS2:
    case AlgoCategory::LCS2:
      // Enabled AlgoCategories
      break;
  }

  // Check if there is a need to use an algorithm for more than two strings
  if (map.contains(ConstraintType::STRINGS_K)) {
    if (const auto p = Get<input::Constraint_LCS_K>(map); p->GetK() != 2) {
      // Algorithms for reconstruction with more than two strings haven't been
      // implemented yet.
      throw std::runtime_error(
          "GetLcsAlgoTypes does not support Constraint_LCS_K");
    }
  }

  // Mapping from the constraint map to the right llcs algorithm
  auto llcs = AlgoType::LLCS2_STD_FL;
  if (map.contains(ConstraintType::MC)) {
    llcs = AlgoType::LLCS2_MC;
  } else if (map.contains(ConstraintType::MC_INC)) {
    if (!avoid_segment_tree) {
      // LLCS2_MC_INC identifies a algorithm that uses segment trees
      llcs = AlgoType::LLCS2_MC_INC;
    } else {
      llcs = AlgoType::LLCS2_MC_INC_E;
    }
  } else if (map.contains(ConstraintType::MC_1C)) {
    llcs = AlgoType::LLCS2_MC_1C;
  } else if (map.contains(ConstraintType::MC_O1C)) {
    throw std::runtime_error("No Default Set for ConstraintType::MC_O1C");
  } else if (map.contains(ConstraintType::MC_O1C_SYNC)) {
    llcs = AlgoType::LLCS2_MC_O1_SYNC;
  }
  if (map.contains(ConstraintType::SIGMA)) {
    if (prefer_max_queues) {
      llcs = AlgoType::LLCS2_SA_MQ;
    } else {
      llcs = AlgoType::LLCS2_SA_RMQ;
    }
  }
  if (map.contains(ConstraintType::SIGMA_R)) {
    if (prefer_max_queues) {
      llcs = AlgoType::LLCS2_SR_MQ;
    } else {
      llcs = AlgoType::LLCS2_SR_RMQ;
    }
  }
  if (map.contains(ConstraintType::SIGMA_L)) {
    if (prefer_max_queues) {
      llcs = AlgoType::LLCS2_SL_MQ;
    } else {
      llcs = AlgoType::LLCS2_SL_RMQ;
    }
  }

  AlgoType main, child;
  switch (category) {
    case AlgoCategory::LLCS2:
      main = llcs;
      child = AlgoType::Unknown;
      break;
    case AlgoCategory::LCS2:
      main = AlgoType::LCS2_RT;
      child = llcs;
      break;
    default:;
      main = AlgoType::Unknown;
      child = AlgoType::Unknown;
  }

  return std::make_pair(main, child);
}

/*******************************************************************************
 * to_json (for BaseConstraint *)
 * @param j nlohmann::json object
 * @param p const std::shared_ptr<BaseConstraint>
 ******************************************************************************/
void to_json(nlohmann::json &j, const BaseConstraint *p) {
  if (p == nullptr) {
    j.clear();
    return;
  }
  j["identifier"] = GetMapTypeToName().at(p->GetType());
  const std::string key = "data";
  switch (p->GetType()) {
    case ConstraintType::MC:
      j[key] = *dynamic_cast<const local::Constraint_MC *>(p);
      break;
    case ConstraintType::MC_INC:
      j[key] = *dynamic_cast<const local::Constraint_MC_INC *>(p);
      break;
    case ConstraintType::Empty:
      j[key] = nullptr;
      break;
    case ConstraintType::MC_1C:
      j[key] = *dynamic_cast<const local::Constraint_MC_1C *>(p);
      break;
    case ConstraintType::MC_O1C:
      j[key] = *dynamic_cast<const local::Constraint_MC_O1C *>(p);
      break;
    case ConstraintType::MC_O1C_SYNC:
      j[key] = *dynamic_cast<const local::Constraint_MC_O1C_SYNC *>(p);
      break;
    case ConstraintType::SIGMA:
      j[key] = *dynamic_cast<const local::Constraint_Sigma *>(p);
      break;
    case ConstraintType::SIGMA_R:
      j[key] = *dynamic_cast<const local::Constraint_Sigma_R *>(p);
      break;
    case ConstraintType::SIGMA_L:
      j[key] = *dynamic_cast<const local::Constraint_Sigma_L *>(p);
      break;
    case ConstraintType::BR:
      j[key] = *dynamic_cast<const global::Constraint_BR *>(p);
      break;
    case ConstraintType::LLCS_KNOWN:
      j[key] = *dynamic_cast<const global::Constraint_LLCS_Known *>(p);
      break;
    case ConstraintType::STRINGS_K:
      j[key] = *dynamic_cast<const input::Constraint_LCS_K *>(p);
      break;
    case ConstraintType::STRINGS_2:
      j[key] = *dynamic_cast<const input::Constraint_LCS_2 *>(p);
      break;
    case ConstraintType::CONST_SIG:
      j[key] = *dynamic_cast<const input::Constraint_Const_Sigma *>(p);
      break;
    default:
      throw std::runtime_error("Unknown ConstraintType");
  }
}
/*******************************************************************************
 * to_json (for std::shared_ptr<BaseConstraint>)
 * @param j nlohmann::json object
 * @param p const std::shared_ptr<BaseConstraint>
 ******************************************************************************/
void to_json(nlohmann::json &j, const std::shared_ptr<BaseConstraint> &p) {
  const BaseConstraint *c = p.get();
  to_json(j, c);
}

/*******************************************************************************
 * from_json (for std::shared_ptr<BaseConstraint>)
 * @param j nlohmann::json object
 * @param p const std::shared_ptr<BaseConstraint>
 ******************************************************************************/
void from_json(nlohmann::json &j, std::shared_ptr<BaseConstraint> &p) {
  if (j.empty()) {
    p = nullptr;
    return;
  }
  std::string name = j["identifier"];
  ConstraintType t = GetMapNameToType().at(name);
  std::string key = "data";
  switch (t) {
    case ConstraintType::Empty: {
      p = nullptr;
      break;
    }

    case ConstraintType::MC: {
      auto constraint = j[key].get<local::Constraint_MC>();
      p = std::make_shared<local::Constraint_MC>(constraint);
      break;
    }
    case ConstraintType::MC_INC: {
      auto constraint = j[key].get<local::Constraint_MC_INC>();
      p = std::make_shared<local::Constraint_MC_INC>(constraint);
      break;
    }
    case ConstraintType::MC_1C: {
      auto constraint = j[key].get<local::Constraint_MC_1C>();
      p = std::make_shared<local::Constraint_MC_1C>(constraint);
      break;
    }
    case ConstraintType::MC_O1C: {
      auto constraint = j[key].get<local::Constraint_MC_O1C>();
      p = std::make_shared<local::Constraint_MC_O1C>(constraint);
      break;
    }
    case ConstraintType::MC_O1C_SYNC: {
      auto constraint = j[key].get<local::Constraint_MC_O1C_SYNC>();
      p = std::make_shared<local::Constraint_MC_O1C_SYNC>(constraint);
      break;
    }
    case ConstraintType::SIGMA: {
      auto constraint = j[key].get<local::Constraint_Sigma>();
      p = std::make_shared<local::Constraint_Sigma>(constraint);
      break;
    }
    case ConstraintType::SIGMA_R: {
      auto constraint = j[key].get<local::Constraint_Sigma_R>();
      p = std::make_shared<local::Constraint_Sigma_R>(constraint);
      break;
    }
    case ConstraintType::SIGMA_L: {
      auto constraint = j[key].get<local::Constraint_Sigma_L>();
      p = std::make_shared<local::Constraint_Sigma_L>(constraint);
      break;
    }
    case ConstraintType::BR: {
      auto constraint = j[key].get<global::Constraint_BR>();
      p = std::make_shared<global::Constraint_BR>(constraint);
      break;
    }
    case ConstraintType::LLCS_KNOWN: {
      auto constraint = j[key].get<global::Constraint_LLCS_Known>();
      p = std::make_shared<global::Constraint_LLCS_Known>(constraint);
      break;
    }
    case ConstraintType::STRINGS_K: {
      auto constraint = j[key].get<input::Constraint_LCS_K>();
      p = std::make_shared<input::Constraint_LCS_K>(constraint);
      break;
    }
    case ConstraintType::STRINGS_2: {
      auto constraint = j[key].get<input::Constraint_LCS_2>();
      p = std::make_shared<input::Constraint_LCS_2>(constraint);
      break;
    }
    case ConstraintType::CONST_SIG: {
      auto constraint = j[key].get<input::Constraint_Const_Sigma>();
      p = std::make_shared<input::Constraint_Const_Sigma>(constraint);
      break;
    }
  }
}

/*******************************************************************************
 * from_json (for ConstraintMap)
 * @param j nlohmann::json object
 * @param map const std::shared_ptr<BaseConstraint>
 * @note nlohmann::json supports the container of the standard library. But it
 * has some problems with the correct constructions of the BaseConstraint due to
 * the share pointers. This function can be called directly to avoid these
 * issues
 ******************************************************************************/
void from_json(nlohmann::json &j, ConstraintMap &map) {
  map.clear();
  for (auto &[key, value] : j.items()) {
    auto t = value[0].get<ConstraintType>();

    std::shared_ptr<BaseConstraint> ptr;
    constraints::from_json(value[1], ptr);
    map.emplace(t, ptr);
  }
}

// /*******************************************************************************
//  * to_json (for ConstraintMap)
//  * @param j nlohmann::json object
//  * @param map const ConstraintMap
//  ******************************************************************************/
// void to_json(nlohmann::json &j, const ConstraintMap &map) {
//   for (const auto &[key, ptr] : map) {
//     j["type_"] = key;
//     j["data"] = ptr;
//   }
// }
//
// /*******************************************************************************
//  * from_json (for ConstraintMap)
//  * @param j nlohmann::json object
//  * @param map const ConstraintMap
//  ******************************************************************************/
// void from_json(const nlohmann::json &j, ConstraintMap *map) {
//
// }

}  // namespace lcs_solver::constraints