#ifndef LCS_SOLVER_INCLUDE_ALGORITHMS_ALGOTYPE_H_
#define LCS_SOLVER_INCLUDE_ALGORITHMS_ALGOTYPE_H_

#include <nlohmann/json.hpp>

namespace lcs_solver::algorithms {

enum class AlgoType : unsigned char {
  Unknown = 0x00,
  LLCS2_STD_FL = 0x01,
  LLCS2_MC = 0x02,
  LLCS2_MC_INC = 0x03,
  LLCS2_MC_INC_E = 0x04,
  LLCS2_MC_1C = 0x05,
  LLCS2_MC_O1_SYNC = 0x06,
  LLCS2_SR_MQ = 0x07,
  LLCS2_SR_RMQ = 0x08,
  LLCS2_SL_MQ = 0x09,
  LLCS2_SL_RMQ = 0x0A,
  LLCS2_SA_MQ = 0x0B,
  LLCS2_SA_RMQ = 0x0C,
  LLCS_STD_FL = 0x0D,
  LCS2_STD_S = 0x0E,
  LCS2_RT = 0x0F,
  LLCS2_STD_H,
  LCS2_STD_H,
  First = LLCS2_STD_FL,
  Last = LLCS2_SA_RMQ
};

NLOHMANN_JSON_SERIALIZE_ENUM(AlgoType, {
    {AlgoType::Unknown, "Unknown"},
    {AlgoType::LLCS2_STD_FL, "LLCS2_STD_FL"},
    {AlgoType::LLCS2_MC, "LLCS2_MC"},
    {AlgoType::LLCS2_MC_INC, "LLCS2_MC_INC"},
    {AlgoType::LLCS2_MC_1C, "LLCS2_MC_1C"},
    {AlgoType::LLCS2_MC_O1_SYNC, "LLCS2_MC_O1_SYNC"},
    {AlgoType::LLCS2_SR_MQ, "LLCS2_SR_MQ"},
    {AlgoType::LLCS2_SR_RMQ, "LLCS2_SR_RMQ"},
    {AlgoType::LLCS2_SL_MQ, "LLCS2_SL_R_LLCS2_SR_MQ"},
    {AlgoType::LLCS2_SL_RMQ, "LLCS2_SL_R_LLCS2_SR_RMQ"},
    {AlgoType::LLCS2_SA_MQ, "LLCS2_SA_MQ"},
    {AlgoType::LLCS2_SA_RMQ, "LLCS2_SA_RMQ"},
    {AlgoType::LLCS_STD_FL, "LLCS_STD_FL"},
    {AlgoType::LCS2_STD_S, "LCS2_STD_S"},
    {AlgoType::LCS2_RT, "LCS2_RT"}
})

}// namespace lcs_solver::algorithms
#endif//LCS_SOLVER_INCLUDE_ALGORITHMS_ALGOTYPE_H_
