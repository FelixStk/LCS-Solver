#ifndef LCS_SOLVER_CONSTRAINTS_CONSTRAINTCATEGORY_H_
#define LCS_SOLVER_CONSTRAINTS_CONSTRAINTCATEGORY_H_

#include "nlohmann/json.hpp"

namespace lcs_solver::constraints{

enum class ConstraintCategory : unsigned char {
  LOCAL = 1,    ///< Longest Common Subsequence (LCS)
  GLOBAL,       ///< LCS with Multiple Gap Constraints given by m-1 tuples (LCS_MC)
  INPUT,        ///< LCS with Multiple Gap Constraints given by m-1 increasing tuples (LCS_MC_INC)
  Invalid,      ///< Constraint type is undefined (this must be the last element in the enum for typesafe reads)
};

NLOHMANN_JSON_SERIALIZE_ENUM( ConstraintCategory, {
    {ConstraintCategory::Invalid, nullptr},
    {ConstraintCategory::LOCAL, "local"},
    {ConstraintCategory::GLOBAL, "global"},
    {ConstraintCategory::INPUT, "input"},
})

}
#endif //LCS_SOLVER_CONSTRAINTS_CONSTRAINTCATEGORY_H_
