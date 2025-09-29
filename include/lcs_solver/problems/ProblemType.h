#ifndef LCS_SOLVER_PROBLEMS_PROBLEMTYPE_H_
#define LCS_SOLVER_PROBLEMS_PROBLEMTYPE_H_

#include <nlohmann/json.hpp>

namespace lcs_solver::problems {

/*******************************************************************************
 * ProblemType
 * @details Defines different types of problems implemented in classes inherited
 *  from BaseProblem. Has two constants First and Last for looping through it.
 ******************************************************************************/
enum class ProblemType : unsigned char {
  LCS_Base = 0x00,
  LCS_Classic = 0x01,
  LCS_MC = 0x02,
  LCS_MC_INC = 0x03,
  LCS_MC_1C = 0x04,
  LCS_MC_O1C_SYNC = 0x05,
  LCS_Sigma_L = 0x06,
  LCS_Sigma_R = 0x07,
  LCS_Sigma = 0x08
};

NLOHMANN_JSON_SERIALIZE_ENUM( ProblemType, {
    {ProblemType::LCS_Base, "LCS_Base"},
    {ProblemType::LCS_Classic, "LCS_Classic"},
    {ProblemType::LCS_MC, "LCS_MC"},
    {ProblemType::LCS_MC_INC, "LCS_MC_INC"},
    {ProblemType::LCS_MC_1C, "LCS_MC_1C"},
    {ProblemType::LCS_MC_O1C_SYNC, "LCS_MC_O1C_SYNC"},
    {ProblemType::LCS_Sigma_L, "LCS_Sigma_L"},
    {ProblemType::LCS_Sigma_R, "LCS_Sigma_R"},
    {ProblemType::LCS_Sigma, "LCS_Sigma"}
})

} //end of namespace lcs_solver::problems
#endif //LCS_SOLVER_PROBLEMS_PROBLEMTYPE_H_
