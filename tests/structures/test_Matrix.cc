/*******************************************************************************
 * @file test_Matrix.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief GTests for the correctness of the Matrix struct
 ******************************************************************************/

#include "util/CommonTypes.h"
#include "structures/Matrix.h"
#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"
#include <string>
#include <algorithm>
#include <random>

using ::lcs_solver::structures::Matrix;

using ::testing::ContainerEq;
using ::testing::Eq;
using ::testing::PrintToString;
using ::testing::StrEq;

namespace {
using ::lcs_solver::util::uint;
TEST(MatrixTest, Constructor) {
  EXPECT_NO_THROW(Matrix<uint>());
  EXPECT_THROW(Matrix<uint>({2, 2, 0}), std::invalid_argument);
}

TEST(MatrixTest, empty) {
  auto M0 = Matrix<uint>();
  auto M1 = Matrix<uint>({2, 2});
  auto M2 = Matrix<uint>(std::vector<size_t>{});
  EXPECT_TRUE(M0.empty());
  EXPECT_FALSE(M1.empty());
  EXPECT_TRUE(M2.empty());
}

TEST(MatrixTest, size) {
  auto M0 = Matrix<uint>();
  auto M1 = Matrix<uint>({2, 2});
  auto M2 = Matrix<uint>({2, 4, 3});
  EXPECT_EQ(M0.size(), 0);
  EXPECT_EQ(M1.size(), 4);
  EXPECT_EQ(M2.size(), 24);
}

TEST(MatrixTest, clear) {
  auto M = Matrix<uint>({2, 2});
  for (size_t i = 0; i < M.size(); i++)
    M[i] = 1;
  M.clear();
  for (size_t i = 0; i < M.size(); i++)
    EXPECT_EQ(M[i], 0);
}

TEST(MatrixTest, linearIndex) {
  auto d = std::vector<size_t>{2, 4, 3};
  auto M = Matrix<uint>(d);
  for (size_t i = 0; i < M.size(); i++)
    M[i] = i; // filled like a normal vector

  for (size_t k = 0, i = 0; k < d[0]; ++k) {
    for (size_t l = 0; l < d[1]; ++l) {
      for (size_t m = 0; m < d[2]; ++m, ++i) {
        auto idx = std::vector<size_t>({k, l, m});
        size_t pos1 = M.linearIndex({k, l, m});
        size_t pos2 = M.linearIndex(idx);
        EXPECT_EQ(pos1, i);
        EXPECT_EQ(pos2, i);
        EXPECT_EQ(M[pos1],
                  i); // meaning if idx1 < idx2 lexicographically then pos1 < pos2, with pos_i = M.linearIndex(idx_i)
      }
    }
  }
}

TEST(MatrixTest, vectorizeIndex) {
  auto d = std::vector<size_t>{2, 4, 4};
  auto M = Matrix<uint>(d);
  for (size_t i = 0; i < M.size(); i++)
    M[i] = i;

  size_t pos = 0;
  for (size_t k = 0; k < d[0]; ++k) {
    for (size_t l = 0; l < d[1]; ++l) {
      for (size_t m = 0; m < d[2]; ++m, ++pos) {
        auto indices = M.vectorizeIndex(pos);
        auto expected = std::vector<size_t>({k, l, m});
        EXPECT_THAT(indices, ContainerEq(expected))
                  << "Call: vectorizeIndex(" << pos << ")";
        EXPECT_THAT(M.linearIndex(indices),
                    Eq(pos)); // M.linearIndex(M.vectorizeIndex(pos) == pos
      }
    }
  }
}

TEST(MatrixTest, stride) {
  auto d = std::vector<size_t>{2, 4, 3};
  auto M = Matrix<uint>(d);
  auto expected =
      std::vector<size_t>({1, d[2], d[2] * d[1], d[2] * d[1] * d[0]});
  EXPECT_THAT(M.stride, ContainerEq(expected));
  EXPECT_THAT(M.stride[M.stride.size() - 1], Eq(M.size()));
}

TEST(MatrixTest, strideDiag) {
  auto d = std::vector<size_t>{2, 4, 3};
  auto D = d.size();
  auto M = Matrix<uint>(d);
  auto product = std::vector<size_t>({1});

  // Another stride test
  for (auto it = d.rbegin(); it != d.rend(); ++it) {
    size_t t = 1;
    for (auto it2 = d.rbegin(); it2 != it + 1; ++it2) {
      t *= *it2;
    }
    product.push_back(t);
  }
  EXPECT_THAT(M.stride, ContainerEq(product));

  // Expect that strideDiag = 1 + d[D-1] + d[D-1]*d[D-2] + ... + (d[D-1]*d[D-2]*...*d[1]) with D = d.size()
  size_t expected = 0;
  for (size_t i = 0; i < D; ++i)
    expected += product[i];
  //EXPECT_THAT(expected, Eq(1 + d[D-1] + d[D-1]*d[D-2]) );
  EXPECT_THAT(M.strideDiag, expected);

}

TEST(MatrixTest, stepDown) {
  auto d = std::vector<size_t>{2, 4, 3};
  auto D = d.size();
  auto M = Matrix<uint>(d);
  for (size_t i = 0; i < M.size(); i++)
    M[i] = i;

  size_t pos = 0;
  for (size_t k = 0; k < d[0]; ++k) {
    for (size_t l = 0; l < d[1]; ++l) {
      for (size_t m = 0; m < d[2]; ++m, ++pos) {
        auto idx = std::vector<size_t>({k, l, m});
        for (size_t i = 0; i < D; ++i) {
          if (idx[D - i - 1] == 0) {
            /**
             * If idx = [a,b,c,d] and the a,b,c or d is that we with i for a stepDown is zero
             * , then M.stepDown(idx, i) = idx
             *
             * idx[D-i-1] is the index for the ith Dimension of the matrix. For example if M[{x,y,z}]
             * and idx = {x,y,z}, then different z define the dimension along i==0 and are at idx[3-1]
             */

            EXPECT_THAT(M.vectorizeIndex(M.stepDown(pos, i)), ContainerEq(idx))
                      << "pos ~ " << PrintToString(idx)
                      << PrintToString(M.vectorizeIndex(pos)) << "\n"
                      << "i == " << i << " should keep position";
          } else {
            /**
             * if idx = [a,b,c,d] and i = 2 then M.stepDown(idx, i) = {a,b,c-1,d}
             */
            --idx[D - i - 1];
            EXPECT_THAT(M.vectorizeIndex(M.stepDown(pos, i)), ContainerEq(idx))
                      << "pos ~ " << PrintToString(idx)
                      << PrintToString(M.vectorizeIndex(pos)) << "\n"
                      << "i == " << i << " should not keep position \n";
            ++idx[D - i - 1];
          }
        }
      }
    }
  }
}

TEST(MatrixTest, stepDownDiag) {
  auto d = std::vector<size_t>{2, 4, 3};
  auto M = Matrix<uint>(d);
  for (size_t i = 0; i < M.size(); i++)
    M[i] = i;

  size_t pos = 0;
  auto zeros = std::vector<size_t>(d.size(), 0);
  for (size_t k = 0; k < d[0]; ++k) {
    for (size_t l = 0; l < d[1]; ++l) {
      for (size_t m = 0; m < d[2]; ++m, ++pos) {
        auto idx = std::vector<size_t>({k, l, m});
        bool isBorder = false;
        for (const auto &i : idx)
          if (i == 0)
            isBorder = true;
        if (isBorder) {
          EXPECT_THAT(M.vectorizeIndex(M.stepDownDiag(pos)), ContainerEq(idx))
                    << "pos = " << pos << " ~ "
                    << PrintToString(M.vectorizeIndex(pos))
                    << " border case" "\n";
        } else {
          auto expected = std::vector<size_t>({k - 1, l - 1, m - 1});
          EXPECT_THAT(M.vectorizeIndex(M.stepDownDiag(pos)),
                      ContainerEq(expected))
                    << "pos ~" << pos << " ~ "
                    << PrintToString(M.vectorizeIndex(pos)) << "\n";
        }
      }
    }
  }
}

TEST(MatrixTest, DebugString) {
  auto M1 = Matrix<uint>();
  auto M2 = Matrix<uint>({2});
  auto M3 = Matrix<uint>({2, 2, 2});
  for (size_t i = 0; i < M3.size(); i++)
    M3[i] = i;
  std::string s1 = "Empty Matrix";
  std::string s2 = "{0, 0}";
  std::string
      s3 = "M[0,0:1,0:1]: \n0 1 \n2 3 \n\nM[1,0:1,0:1]: \n4 5 \n6 7 \n\n";
  EXPECT_THAT(M1.DebugString(), StrEq(s1))
            << "The raw DebugString is:\n" << M1.DebugString();
  EXPECT_THAT(M2.DebugString(), StrEq(s2))
            << "The raw DebugString is:\n" << M2.DebugString();
  EXPECT_THAT(M3.DebugString(), StrEq(s3))
            << "The raw DebugString is:\n" << M3.DebugString();
}

TEST(MatrixTest, lastIndex) {
  auto M = Matrix<uint>({2, 4, 3});
  auto expected = std::vector<size_t>({1, 3, 2});
  EXPECT_THAT(M.lastIndex, ContainerEq(expected));
}

TEST(MatrixTest, incr) {
  auto M1 = Matrix<uint>({2, 3, 4});
  auto indices1 = std::vector<size_t>({0, 0, 0});
  auto expected1 = std::vector<size_t>({0, 0, 1});
  auto expected2 = std::vector<size_t>({0, 1, 0});
  EXPECT_TRUE(M1.incr(indices1));
  EXPECT_THAT(indices1, ContainerEq(expected1));

  M1.incr(indices1);
  M1.incr(indices1);
  M1.incr(indices1);
  EXPECT_THAT(indices1, ContainerEq(expected2));

  auto M2 = Matrix<uint>({2, 2});
  auto indices2 = std::vector<size_t>({1, 1});
  EXPECT_FALSE(M2.incr(indices2));

  auto M3 = Matrix<uint>();
  auto indices3 = std::vector<size_t>({1, 1});
  EXPECT_FALSE(M3.incr(indices3));

  auto M4 = Matrix<uint>({1, 1, 1});
  auto indices4 = std::vector<size_t>({0, 0, 0});
  EXPECT_FALSE(M4.incr(indices4));

}

TEST(MatrixTest, getIdxAbs) {
  auto d = std::vector<size_t>{2, 4, 3};
  auto M = Matrix<uint>(d);
  for (size_t i = 0; i < M.size(); i++)
    M[i] = i;

  for (size_t k = 0; k < d[0]; ++k) {
    for (size_t l = 0; l < d[1]; ++l) {
      for (size_t m = 0; m < d[2]; ++m) {
        // for each cell with the index {k,l,m}:
        for (size_t dim = 0; dim < d.size(); ++dim) {
          for (size_t i = 0; i < d[dim]; ++i) {
            auto idx = std::vector<size_t>({k, l, m});
            size_t start = M.linearIndex(idx);
            idx[dim] = i;
            size_t expected = M.linearIndex(idx);
            size_t pos = M.getIdxAbs(start, dim, i);
            EXPECT_EQ(pos, expected)
                      << "start: " << start << std::endl
                      << "dim: " << dim << std::endl
                      << "i: " << i << std::endl
                      << M.DebugString();
          }
        }
      }
    }
  }

}

TEST(MatrixTest, getIdxRel) {
  auto d = std::vector<size_t>{2, 4, 3};
  auto M = Matrix<uint>(d);
  for (size_t i = 0; i < M.size(); i++)
    M[i] = i;

  size_t pos = 0;
  for (size_t k = 0; k < d[0]; ++k) {
    for (size_t l = 0; l < d[1]; ++l) {
      for (size_t m = 0; m < d[2]; ++m, ++pos) {
        // for each cell with the index {k,l,m}:
        for (size_t dim = 0; dim < d.size(); ++dim) {
          for (size_t i = 0; i < d[dim]; ++i) {
            auto idx = std::vector<size_t>({k, l, m});
            size_t start = M.linearIndex(idx);
            int move = static_cast<int>(i) - idx[dim];
            idx[dim] = i;
            size_t expected = M.linearIndex(idx);
            EXPECT_EQ(M.getIdxRel(start, dim, move), expected);
          }
        }
      }
    }
  }

}

TEST(MatrixTest, getIdxFst) {
  auto d = std::vector<size_t>{2, 4, 3};
  auto M = Matrix<uint>(d);
  for (size_t i = 0; i < M.size(); i++)
    M[i] = i;

  auto expected = std::vector<size_t>(d.size(), 0);
  for (size_t k = 0; k < d[0]; ++k) {
    expected[0] = M.linearIndex({k, 0, 0});
    for (size_t l = 0; l < d[1]; ++l) {
      expected[1] = M.linearIndex({k, l, 0});
      for (size_t m = 0; m < d[2]; ++m) {
        expected[2] = M.linearIndex({k, l, m});
        // for each cell with the index {k,l,m}:
        std::vector<size_t> idx = {k, l, m};
        size_t pos = M.linearIndex(idx);
        for (size_t i = 0; i < expected.size(); ++i)
          EXPECT_EQ(M.getIdxFst(pos, i), expected[i])
                    << "pos: " << pos << "\n"
                    << "i: " << i << "\n"
                    << "idx(pos): " << PrintToString(idx) << "\n"
                    << M.DebugString();
      }
    }
  }
}

TEST(MatrixTest, getIdxVectorized) {
  auto d = std::vector<size_t>{2, 4, 3};
  auto M = Matrix<uint>(d);

  for (size_t k = 0; k < d[0]; ++k) {
    for (size_t l = 0; l < d[1]; ++l) {
      for (size_t m = 0; m < d[2]; ++m) { // for each cell {k,l,m}:
        std::vector<size_t> idx = {k, l, m};
        size_t pos = M.linearIndex(idx);
        for (size_t i = 0; i < d.size(); ++i)
          EXPECT_EQ(M.vectorizeIndex(pos, i), idx[i])
                    << "pos: " << pos << "\n"
                    << "i: " << i << "\n"
                    << "idx(pos): " << PrintToString(idx) << "\n";
      }
    }
  }
}

}  // end of anonym namespace