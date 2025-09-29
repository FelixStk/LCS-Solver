#include "structures/Matrix.h"// Include your Matrix and MatrixConstIter headers
#include "structures/MatrixConstIter.h"
#include "gmock/gmock-more-matchers.h"
#include "gtest/gtest.h"
#include <vector>

namespace lcs_solver::structures {

using ::testing::ContainerEq;

// Test fixture for MatrixConstIter
class MatrixIterTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create a 2x3x4 matrix for testing
    matrix_ = Matrix<int>({2, 3, 4});
    for (size_t i = 0; i < matrix_.size(); ++i) {
      matrix_[i] = static_cast<int>(i);// Fill with sequential values
    }
  }

  Matrix<int> matrix_;
};

// Test parameterized constructor and traversal
TEST_F(MatrixIterTest, ForwardTraversal) {
  std::vector<std::vector<size_t>> axis_orders;
  axis_orders.push_back({0, 1, 2});
  axis_orders.push_back({1, 0, 2});
  axis_orders.push_back({0, 2, 1});
  axis_orders.push_back({2, 0, 1});
  axis_orders.push_back({1, 2, 0});
  axis_orders.push_back({2, 1, 0});
  for (auto axis_order : axis_orders) {
    std::vector<std::vector<size_t>> expected_idx;
    std::vector<size_t> pos(3, 0);
    for (; pos[0] < matrix_.dim.at(axis_order[0]); ++pos[0], pos[1] = 0, pos[2] = 0) {
      for (; pos[1] < matrix_.dim.at(axis_order[1]); ++pos[1], pos[2] = 0) {
        for (; pos[2] < matrix_.dim.at(axis_order[2]); ++pos[2]) {
          std::vector<size_t> idx(3, 0);
          for (size_t i = 0; i < pos.size(); ++i) {
            idx[axis_order[i]] = pos[i];
          }
          expected_idx.push_back(idx);
        }
      }
    }
    std::vector<int> expected_values;
    for (const auto &idx : expected_idx) {
      expected_values.push_back(matrix_.linearIndex(idx));
    }

    // Verify that the iterator produces the expected values
    std::vector<int> values;
    auto expectedIter = expected_values.begin();
    auto iter = MatrixConstIter<int>(matrix_, axis_order);
    auto end = MatrixConstIter<int>();
    while (iter != end) {
      values.push_back(*iter);
      // EXPECT_EQ(*iter, *expectedIter) << "it:" << iter.DebugString() << std::endl;
      ++expectedIter;
      ++iter;
    }
    EXPECT_THAT(values, ContainerEq(expected_values))
        << "axis_order" << testing::PrintToString(axis_order) << std::endl
        << "matrix: " << matrix_.DebugString() << std::endl
        << "expected_idx" << testing::PrintToString(expected_idx) << std::endl;

    EXPECT_EQ(iter, MatrixConstIter<int>());// Verify that the iterator reaches the end
  }
}

TEST_F(MatrixIterTest, BackWardTraversal) {
  std::vector<std::vector<size_t>> axis_orders;
  axis_orders.push_back({0, 1, 2});
  axis_orders.push_back({1, 0, 2});
  axis_orders.push_back({0, 2, 1});
  axis_orders.push_back({2, 0, 1});
  axis_orders.push_back({1, 2, 0});
  axis_orders.push_back({2, 1, 0});
  for (auto axis_order : axis_orders) {
    std::vector<std::vector<size_t>> expected_idx;
    std::vector<size_t> pos(3, 0);
    for (; pos[0] < matrix_.dim.at(axis_order[0]); ++pos[0], pos[1] = 0, pos[2] = 0) {
      for (; pos[1] < matrix_.dim.at(axis_order[1]); ++pos[1], pos[2] = 0) {
        for (; pos[2] < matrix_.dim.at(axis_order[2]); ++pos[2]) {
          std::vector<size_t> idx(3, 0);
          for (size_t i = 0; i < pos.size(); ++i) {
            idx[axis_order[i]] = pos[i];
          }
          expected_idx.push_back(idx);
        }
      }
    }
    std::vector<int> expected_values;
    std::ranges::reverse(expected_idx);
    for (const auto &idx : expected_idx) {
      expected_values.push_back(matrix_.linearIndex(idx));
    }

    // Verify that the iterator produces the expected values
    std::vector<int> values;
    auto expectedIter = expected_values.crbegin();
    auto iter = MatrixConstIter<int>(matrix_, axis_order, true);
    auto crend = MatrixConstIter<int>();
    while (iter != crend) {
      values.push_back(*iter);
      // EXPECT_EQ(*iter, *expectedIter) << "it:" << iter.DebugString() << std::endl;
      ++expectedIter;
      ++iter;
    }
    EXPECT_THAT(values, ContainerEq(expected_values))
        << "axis_order" << testing::PrintToString(axis_order) << std::endl
        << "matrix: " << matrix_.DebugString() << std::endl
        << "expected_idx" << testing::PrintToString(expected_idx) << std::endl;

    EXPECT_EQ(iter, MatrixConstIter<int>());// Verify that the iterator reaches the end
  }
}

TEST_F(MatrixIterTest, EqualityAndInequality) {
  std::vector<size_t> axis_order = {1, 0, 2};
  MatrixConstIter iter1(matrix_, axis_order);
  MatrixConstIter iter2(matrix_, axis_order);

  EXPECT_EQ(iter1, iter2);  // Iterators should be equal initially
  EXPECT_NE(++iter1, iter2);// Increment one iterator and check inequality
  EXPECT_EQ(iter1, ++iter2);// Increment the second iterator and check equality
}

// Test dereferencing
TEST_F(MatrixIterTest, Dereferencing) {
  std::vector<size_t> axis_order = {1, 0, 2};
  MatrixConstIter<int> iter(matrix_, axis_order);

  // Check the first value
  EXPECT_EQ(*iter, 0);

  // Check the second value
  ++iter;
  EXPECT_EQ(*iter, 1);
}

TEST_F(MatrixIterTest, PostIncrement) {
  std::vector<size_t> axis_order = {1, 0, 2};
  MatrixConstIter<int> iter(matrix_, axis_order);
  MatrixConstIter<int> old_iter = iter++;
  EXPECT_EQ(*old_iter, 0);
  EXPECT_EQ(*iter, 1);
}

TEST_F(MatrixIterTest, EndIterator) {
  const std::vector<size_t> axis_order = {1, 0, 2};
  MatrixConstIter<int> iter(matrix_, axis_order);
  for (size_t i = 0; i < matrix_.size(); ++i) {
    ++iter;
  }

  // Verify that the iterator is equal to a default-constructed iterator
  auto end = MatrixConstIter<int>();
  EXPECT_EQ(iter == end, true)
      << "matrix: " << matrix_.DebugString() << "\n"
      << "matrix.size: " << matrix_.size() << "\n"
      << "iter: " << iter.DebugString() << "\n"
      << "end: " << end.DebugString() << std::endl;
}

TEST_F(MatrixIterTest, InvalidAxisOrder) {
  std::vector<size_t> invalid_axis_order = {1, 0};// Missing axis 2
  EXPECT_THROW(MatrixConstIter<int>(matrix_, invalid_axis_order), std::invalid_argument);
}

}// namespace lcs_solver::structures