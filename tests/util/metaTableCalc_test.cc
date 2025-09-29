/*******************************************************************************
 * @file metaTableCalc_test.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief GTests make_log2_table and make_pow2_table
 ******************************************************************************/

#include "util/MetaTableCalc.h"

#include <cmath>

#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace {

using ::lcs_solver::util::Meta_Log2_table;
using ::lcs_solver::util::Meta_Pow2_table;
using ::lcs_solver::util::Log2_table;
using ::lcs_solver::util::Pow2_table;

TEST(MetaPow2Log2, MetaLog2Table_LEN0){
  constexpr auto table = Meta_Log2_table<0>();
  EXPECT_EQ(table.arr.empty(), true);
}

TEST(MetaPow2Log2, MetaLog2Table_LEN1){
  constexpr size_t N = 1;
  constexpr auto LOG2 = Meta_Log2_table<N>();
  EXPECT_EQ(LOG2.arr.size(), N);
  EXPECT_EQ(LOG2[1], 0);
}

TEST(MetaPow2Log2, MetaLog2Table_LEN2){
  constexpr size_t N = 2;
  constexpr auto LOG2 = Meta_Log2_table<N>();
  EXPECT_EQ(LOG2.arr.size(), N);
  EXPECT_EQ(LOG2[1], 0);
  EXPECT_EQ(LOG2[2], 1);
}

TEST(MetaPow2Log2, MetaLog2Table_LEN256){
  constexpr size_t N = 256;
  constexpr auto LOG2 = Meta_Log2_table<N>();
  EXPECT_EQ(LOG2.arr.size(), N);
  for(size_t i = 1; i<=N; ++i){
    EXPECT_EQ(LOG2[i], static_cast<size_t>(log2(i)))
    << "i: " << i << std::endl;
  }
}

TEST(MetaPow2Log2, MetaPow2Table_LEN1){
  constexpr size_t N = 1;
  constexpr auto POW2 = Meta_Pow2_table<N>();
  EXPECT_EQ(POW2.arr.size(), N);
  EXPECT_EQ(POW2[0], 1);
}

TEST(MetaPow2Log2, MetaPow2Table_LEN24){
  constexpr size_t N = 24;
  constexpr auto POW2 = Meta_Pow2_table<N>();
  EXPECT_EQ(POW2.arr.size(), N);
  for(size_t i = 0; i<N; ++i){
    EXPECT_EQ(POW2[i], static_cast<size_t>(pow(2,i)))
              << "i: " << i << std::endl;
  }
}


TEST(MetaPow2Log2, Log2Table_LEN0){
  auto table = Log2_table<size_t>(0);
  EXPECT_EQ(table.arr.empty(), true);
}

TEST(MetaPow2Log2, Log2Table_LEN1){
  constexpr size_t N = 1;
  auto LOG2 = Log2_table<size_t >(N);
  EXPECT_EQ(LOG2.arr.size(), N);
  EXPECT_EQ(LOG2[1], 0);
}

TEST(MetaPow2Log2, Log2Table_LEN2){
  constexpr size_t N = 2;
  auto LOG2 = Log2_table<size_t >(N);
  EXPECT_EQ(LOG2.arr.size(), N);
  EXPECT_EQ(LOG2[1], 0);
  EXPECT_EQ(LOG2[2], 1);
}

TEST(MetaPow2Log2, Log2Table_LEN256){
  constexpr size_t N = 256;
  auto LOG2 = Log2_table<size_t >(N);
  EXPECT_EQ(LOG2.arr.size(), N);
  for(size_t i = 1; i<=N; ++i){
    EXPECT_EQ(LOG2[i], static_cast<size_t>(log2(i)))
              << "i: " << i << std::endl;
  }
}

TEST(MetaPow2Log2, Pow2Table_LEN1){
  constexpr size_t N = 1;
  auto POW2 = Pow2_table<size_t>(N);
  EXPECT_EQ(POW2.arr.size(), N);
  EXPECT_EQ(POW2[0], 1);
}

TEST(MetaPow2Log2, Pow2Table_LEN24){
  constexpr size_t N = 24;
  auto POW2 = Pow2_table<size_t>(N);
  EXPECT_EQ(POW2.arr.size(), N);
  for(size_t i = 0; i<N; ++i){
    EXPECT_EQ(POW2[i], static_cast<size_t>(pow(2,i)))
              << "i: " << i << std::endl;
  }
}

}