#ifndef LCS_SOLVER_INCLUDE_ALGORITHMS_ALGOCATEGORY_H_
#define LCS_SOLVER_INCLUDE_ALGORITHMS_ALGOCATEGORY_H_

#include <nlohmann/json.hpp>

namespace lcs_solver::algorithms {

enum class AlgoCategory : unsigned char {
  Empty = 0x00,
  LLCS = 0x01,
  LLCS2 = 0x02,
  LCS = 0x03,
  LCS2 = 0x04,
  Oracle = 0x05,
  Online = 0x06,
  Other = 0x07
};

NLOHMANN_JSON_SERIALIZE_ENUM( AlgoCategory, {
    {AlgoCategory::Empty, nullptr},
    {AlgoCategory::LLCS, "LLCS"},
    {AlgoCategory::LLCS2, "LLCS2"},
    {AlgoCategory::LCS, "LCS"},
    {AlgoCategory::Oracle, "Oracle"},
    {AlgoCategory::Online, "Online"},
    {AlgoCategory::Other, "Other"}
})

}
#endif //LCS_SOLVER_INCLUDE_ALGORITHMS_ALGOCATEGORY_H_
