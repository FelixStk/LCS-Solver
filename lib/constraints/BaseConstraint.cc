/******************************************************************************
 * @file BaseConstraint.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Impl. BaseConstraint: Getter for category, type & constructor
 *****************************************************************************/

#include "constraints/BaseConstraint.h"

#include <fstream>
#include <nlohmann/json.hpp>

#include "constraints/ConstraintMap.h"
#include "constraints/ConstraintType.h"

namespace lcs_solver::constraints {

BaseConstraint::BaseConstraint(ConstraintCategory cat, ConstraintType t)
    : category_(cat), type_(t) {}

std::shared_ptr<BaseConstraint> BaseConstraint::FromJson(
    const std::filesystem::path& path) {
  nlohmann::json j(path);
  std::shared_ptr<BaseConstraint> ptr;
  constraints::from_json(j, ptr);
  return ptr;
}

void BaseConstraint::ToJson(const std::filesystem::path& path, const BaseConstraint* constraint) {
  std::ofstream fs;
  fs.open(path, std::ofstream::out | std::ofstream::trunc);
  if (!fs) {
    throw std::runtime_error("Could not open file " + path.string() + " for writing");
  }
  nlohmann::json j;
  constraints::to_json(j, constraint);
  fs << j.dump(4);
  fs.close();
}

void BaseConstraint::ToJson(const std::filesystem::path& path) const {
  ToJson(path, this);
}

ConstraintCategory BaseConstraint::GetCategory() const {
  return category_;
}

ConstraintType BaseConstraint::GetType() const {
  return type_;
}

}// namespace lcs_solver::constraints

//== Old Code for subsequence checking =========================================
// bool LCS_Gap::isSubsequence(std::string sub, std::string str)
// {
//     int subLength = sub.length();
//     int strLength = str.length();
//
//     int i = 0; // Index for the subsequence
//     int j = 0; // Index for the string
//     int g = 0; // Count of characters between aligned characters
//     bool isFirstMatch = true; // Indicator for the first match
//
//     while (i < subLength && j < strLength) {
//         if (sub[i] == str[j]) {
//             if(!isFirstMatch && (g < gap[i-1].first || g > gap[i-1].second)){
//                 return false;
//             } else {
//                 isFirstMatch = false;
//             }
//             ++i; // Move to the next character in the subsequence
//             g = 0; // Reset the gap
//         }
//         else {
//             if (!isFirstMatch) { // Count gaps only after the first match
//                 ++g; // Increase the gap
//             }
//         }
//         ++j; // Always move to the next character in the string
//     }
//     return i == subLength;  // Does not consider the gap after the last matching
// }