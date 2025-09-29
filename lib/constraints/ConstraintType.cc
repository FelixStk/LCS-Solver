/******************************************************************************
 * @file ConstraintType.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Implementation of a type to name, description map
 * @note If you add a constraint, remember to also add a switch case in
 * ConstraintMap.cc for JSON conversion of std::shared_ptr<BaseConstraint>
 *****************************************************************************/

#include "constraints/ConstraintType.h"

#include "constraints/ConstraintMap.h"
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
 * Definition of Mapping: ConstraintType => Names
 * @return Reference std::map<ConstraintType, std::string_view>
 ******************************************************************************/
const std::map<ConstraintType, std::string_view>& GetMapTypeToName() {
  static std::map<ConstraintType, std::string_view> map_name = {
    {ConstraintType::MC, local::Constraint_MC::kName},
    {ConstraintType::MC_INC, local::Constraint_MC_INC::kName},
    {ConstraintType::MC_1C, local::Constraint_MC_1C::kName},
    {ConstraintType::MC_O1C, local::Constraint_MC_O1C::kName},
    {ConstraintType::MC_O1C_SYNC, local::Constraint_MC_O1C_SYNC::kName},
    {ConstraintType::SIGMA, local::Constraint_Sigma::kName},
    {ConstraintType::SIGMA_R, local::Constraint_Sigma_R::kName},
    {ConstraintType::SIGMA_L, local::Constraint_Sigma_L::kName},
    {ConstraintType::BR, global::Constraint_BR::kName},
    {ConstraintType::LLCS_KNOWN, global::Constraint_LLCS_Known::kName},
    {ConstraintType::STRINGS_K, input::Constraint_LCS_K::kName},
    {ConstraintType::STRINGS_2, input::Constraint_LCS_2::kName},
    {ConstraintType::CONST_SIG, input::Constraint_Const_Sigma::kName}};
  return map_name;
}

/*******************************************************************************
 * Definition of Mapping: ConstraintType => Description
 * @return Reference std::map<ConstraintType, std::string_view>
 ******************************************************************************/
const std::map<ConstraintType, std::string_view>& GetMapTypeToDesc() {
  static std::map<ConstraintType, std::string_view> map_description = {
    {ConstraintType::MC, local::Constraint_MC::kDescription},
    {ConstraintType::MC_INC, local::Constraint_MC_INC::kDescription},
    {ConstraintType::MC_1C, local::Constraint_MC_1C::kDescription},
    {ConstraintType::MC_O1C, local::Constraint_MC_O1C::kDescription},
    {ConstraintType::MC_O1C_SYNC, local::Constraint_MC_O1C_SYNC::kDescription},
    {ConstraintType::SIGMA, local::Constraint_Sigma::kDescription},
    {ConstraintType::SIGMA_R, local::Constraint_Sigma_R::kDescription},
    {ConstraintType::SIGMA_L, local::Constraint_Sigma_L::kDescription},
    {ConstraintType::BR, global::Constraint_BR::kDescription},
    {ConstraintType::LLCS_KNOWN, global::Constraint_LLCS_Known::kDescription},
    {ConstraintType::STRINGS_K, input::Constraint_LCS_K::kDescription},
    {ConstraintType::STRINGS_2, input::Constraint_LCS_2::kDescription},
    {ConstraintType::CONST_SIG, input::Constraint_Const_Sigma::kDescription}};
  return map_description;
}

/*******************************************************************************
 * Getter for a Mapping: Names => ConstraintType
 * @return Reference std::map<std::string_view, ConstraintType>
 ******************************************************************************/
const std::map<std::string_view, ConstraintType>& GetMapNameToType() {
  static std::map<std::string_view, ConstraintType> map_r_name;
  if (map_r_name.empty()) {
    for (const auto &[key, value] : GetMapTypeToName()) {
      map_r_name[value] = key;
    }
  }
  return map_r_name;
}

void to_json(nlohmann::json &j, const ConstraintType &c) {
  j["ConstraintType"] = GetMapTypeToName().at(c);
}

void from_json(const nlohmann::json &j, ConstraintType &c) {
  std::string key = j["ConstraintType"];
  c = GetMapNameToType().at(key);
}

}// namespace lcs_solver::constraints