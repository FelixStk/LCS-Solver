#ifndef LCS_SOLVER_TEST_LLCS_TEST_H
#define LCS_SOLVER_TEST_LLCS_TEST_H

#include "util/AlgoTestParam.h"

namespace lcs_solver::testing::llcs {

using util::AlgoParam;
using constraints::ConstraintType;

std::vector<AlgoParam> genHardCodedParams(ConstraintType t, size_t len=8, bool strict = true);

}
#endif //TEST_H
