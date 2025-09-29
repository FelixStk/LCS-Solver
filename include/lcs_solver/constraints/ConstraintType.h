#ifndef LCS_SOLVER_CONSTRAINTS_CONSTRAINTTYPE_H_
#define LCS_SOLVER_CONSTRAINTS_CONSTRAINTTYPE_H_

#include <nlohmann/json.hpp>

namespace lcs_solver::constraints {
enum class ConstraintType : unsigned char {
  Empty = 0x00,
  MC = 0x01,
  MC_INC = 0x02,
  MC_1C = 0x03,
  MC_O1C = 0x04,
  MC_O1C_SYNC = 0x05,
  SIGMA = 0x65,
  SIGMA_R = 0x07,
  SIGMA_L = 0x08,
  BR = 0x09,
  LLCS_KNOWN = 0x0A,
  STRINGS_K = 0x0B,
  STRINGS_2 = 0x0C,
  CONST_SIG = 0x0D,
};

const std::map<std::string_view, ConstraintType> &GetMapNameToType();
const std::map<ConstraintType, std::string_view> &GetMapTypeToName();
const std::map<ConstraintType, std::string_view> &GetMapTypeToDesc();

void to_json(nlohmann::json &j, const ConstraintType &c);
void from_json(const nlohmann::json &j, ConstraintType &c);

}// namespace lcs_solver::constraints
#endif//LCS_SOLVER_INCLUDE_CONSTRAINTS_CONSTRAINTTYPE_H_
