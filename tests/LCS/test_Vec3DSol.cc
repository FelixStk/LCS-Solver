/*******************************************************************************
 * @file test_Vec3DSol.cc
 * @author Felix Steinkopp
 * @version 1.1
 * @brief GTests to check construction and use of the Vector3DSolution class
 ******************************************************************************/

#include "algorithms/solutions/Vector3DSolution.h"

#include <random>

#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace lcs_solver::testing::lcs  {

using lcs_solver::algorithms::solutions::Vector3DSolution;
using lcs_solver::util::to_String;
using lcs_solver::util::Symbol;
using uint = Vector3DSolution::uint;
using Matrix = std::vector<std::vector<uint>>;
using Points = Vector3DSolution::Points;

TEST(Vector3DSolution, IteratorUsage) {
  // Define Vector3DSolution Solution
  Vector3DSolution::StrPtrVec spv;
  spv.emplace_back(std::make_shared<const util::String >(to_String <Symbol>("ABCD")));
  spv.emplace_back(std::make_shared<const util::String >(to_String <Symbol>("ACBD")));
  Matrix m1 = {{0, 0}, {1, 2}, {3, 3}};
  Matrix m2 = {{0, 0}, {2, 1}, {3, 3}};
  auto p1 = Points(spv, m1);
  auto p2 = Points(spv, m1);
  auto pointVec = std::vector{p1, p2};

  // Use Iterator
  auto expected = Vector3DSolution(pointVec);
  auto it = expected.begin();
  auto end = expected.end();

  std::vector<Points> vec = {};
  for (; it != end; ++it) {
    Points sol = *it;
    vec.push_back(sol);
  }
  Vector3DSolution testVec = Vector3DSolution(vec);

  EXPECT_EQ(testVec, expected)
      << "test Vector3DSolution:\n" << testVec.DebugString() << "\n"
      << "expected Vector3DSolution:\n" << expected.DebugString() << std::endl;
}

} // namespace
