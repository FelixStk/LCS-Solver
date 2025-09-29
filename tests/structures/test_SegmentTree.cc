/******************************************************************************
 * @file test_SegmentTree.cc
 * @author Felix Steinkopp
 * @version 1.1
 * @brief Test the correctness of the SegmentTree class
 * @note SegmentTree uses function pointers: q:TxT->T is the query function
 * and u:TxT->T in U are update functions. Updates are done elementwise
 * node[i] = u(data[i], x). (T,q) should be a commutative semi group this means
 * q(a,b) = q(b,a). Additionally its assumed that, functions in U are be closed
 * under composition and distribute over q: u(q(a,b)) = q(u(a), u(b))
 * This notion is based on: https://doi.org/10.4230/LIPIcs.ITCS.2024.35
 ******************************************************************************/

#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

#include "structures/SegmentTree.h"

namespace {

using ::lcs_solver::structures::MemoryOrder;
using ::lcs_solver::structures::SegmentTree;

using ::testing::Lt;
using ::testing::Le;
using ::testing::Eq;
using ::testing::ContainerEq;
using ::testing::PrintToString;
using ::testing::ValuesIn;

//== Parameter type for the tests =========================================
using T = int;
using Pair = std::pair<T, T>;
using SegmentTreeTestParam = std::tuple<
    std::string,                                  // 0: name
    std::vector<T>,                               // 1: initialValues
    std::vector<std::tuple<size_t, T>>,           // 2: updates
    std::vector<std::tuple<size_t, size_t, T>>,   // 3: rangeUpdates
    std::vector<size_t>,                          // 4: queries
    std::vector<std::tuple<size_t, size_t>>       // 5: rangeQueries
>;

//== Definition of Testparameter ===============================================
std::vector<SegmentTreeTestParam> getHardCodedParams() {
  std::vector<SegmentTreeTestParam> params = {
      {
          "EmptyInitialValues", // name of parameter struct
          {},  // Vector for Initialization
          {{2, 10}},  // Simple Updates (pos, val}
          {{1, 3, 5}},  // Ranged Updates (frontIndex, backIndex, val)
          {0},  // Simple Queries  (index)
          {{0, 0}}  // Ranged Queries (frontIndex, backIndex)
      },
      {
          "HardCodedValue1",        // Name
          std::vector<int>(10, 0),  // Initialization Values
          {{2, 10}},                // Simple Updates
          {{5, 8, 48}, {3, 9, 22}}, // Ranged Updates
          {0},                      // Simple Queries
          {{0, 0}}                  // Ranged Queries
      },
      {
          "HardCodedValue2",   // Name
          {80, 61, 5, 9, 91, 85, 33, 16, 36, 30}, // Initialization Values
          {{2, 10}},              // Simple Updates
          {{5, 8, 48}, {3, 9, 22}}, // Ranged Updates
          {0},                    // Simple Queries
          {{1, 4}}                // Ranged Queries
      }
  };
  return params;
}

/*******************************************************************************
 * genRndParameters
 * @param initSeed unsigned int used for random number generation
 * @param nParams size_t for the number of parameters to generate
 * @param nInit size_t for the number of values to initialize in the data vector
 * @param nSimpleUpdates size_t for the number of point updates
 * @param nRangeUpdates size_t for the number of range updates
 * @param nSimpleQueries size_t for the number of point queries
 * @param nRangeQueries size_t for the number of range queries
 * @param r std::pair<T,T> range of possible values used data initialization
 * @return std::vector<SegmentTreeTestParam>
 ******************************************************************************/
std::vector<SegmentTreeTestParam> genRndParameters(unsigned int initSeed,
                                                   size_t nParams = 5,
                                                   size_t nInit = 10,
                                                   size_t nSimpleUpdates = 3,
                                                   size_t nRangeUpdates = 2,
                                                   size_t nSimpleQueries = 3,
                                                   size_t nRangeQueries = 3,
                                                   Pair r = {1, 100}) {
  std::uniform_int_distribution<T> dist(r.first, r.second);

  std::vector<SegmentTreeTestParam> params;
  params.reserve(nParams);
  for (size_t i = 0; i < nParams; ++i) {
    unsigned int seed = initSeed + i;
    std::mt19937 gen(seed);
    SegmentTreeTestParam x;
    auto &[name, vec, us, ur, sq, rq] = x;
    // Name
    name = "Random_" + std::to_string(seed);
    // Values for initialization
    std::generate_n(std::back_inserter(vec), nInit,
                    [&]() { return dist(gen); });
    // Simple Update (changes one element)
    std::generate_n(std::back_inserter(us), nSimpleUpdates,
                    [&]() -> std::tuple<size_t, int> {
                      return {dist(gen) % nInit, dist(gen)};
                    });
    // Range Updates
    std::generate_n(std::back_inserter(ur), nRangeUpdates,
                    [&]() -> std::tuple<size_t, size_t, int> {
                      size_t start = dist(gen) % nInit;
                      size_t end = start + (dist(gen) % (nInit - start));
                      auto value = dist(gen);
                      return {start, end, value};
                    });
    // Simple Query
    std::generate_n(std::back_inserter(sq), nSimpleQueries,
                    [&]() -> size_t {
                      return dist(gen) % nInit;
                    });
    // Range Query
    std::generate_n(std::back_inserter(rq), nRangeQueries,
                    [&]() -> std::tuple<size_t, size_t> {
                      size_t start = dist(gen) % nInit;
                      size_t end = start + (dist(gen) % (nInit - start));
                      return {start, end};
                    });
    params.push_back(std::move(x));
  }
  return params;
}


//==== Initialize Test Suites ==================================================

class SegmentTreeTest : public testing::TestWithParam<SegmentTreeTestParam> {};

const auto rndSeed = std::random_device{}();
//const unsigned int seed = 2275096221;

auto hardCodedParameters = getHardCodedParams();
auto rndParameters = genRndParameters(
    rndSeed,  // seed
    5,        // nParam
    10,       // nInit
    3,        // nSimpleUpdates
    2,        // nRangeUpdates
    3,        // nSimpleQueries
    3,        // nRangeQueries
    {1, 100}  // value range
);

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters, SegmentTreeTest,
    ValuesIn(hardCodedParameters),
    [](const testing::TestParamInfo<SegmentTreeTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    RandomParameters, SegmentTreeTest,
    ValuesIn(rndParameters),
    [](const testing::TestParamInfo<SegmentTreeTest::ParamType> &info) {
      return std::get<0>(info.param);
    });


//==== Parametrized Test =======================================================

TEST_P(SegmentTreeTest, InitializationTest) {
  auto &data = std::get<1>(GetParam());
  auto max = [](int a, int b) { return a < b ? a : b; };
  auto order = MemoryOrder::Euler;
  SegmentTree<int> st(max, data, order);
  size_t expected = data.empty() ? 0 : data.size();
  size_t out = st.size();
  EXPECT_THAT(out, Eq(expected));
}

TEST_P(SegmentTreeTest, SimpleUpdateTest) {
  auto len = std::get<1>(GetParam()).size();
  auto updates = std::get<2>(GetParam());
  auto max = [](int a, int b) { return a > b ? a : b; };
  auto set = [](int a, int b) { return b; };
  SegmentTree<int> st(max, len, MemoryOrder::Euler, 0);

  auto expected = std::vector<int>(len, 0);
  auto out = std::vector<int>(len, 0);
  for (const auto &update : updates) {
    auto [index, val] = update;
    st.update(set, index, val);
    if (len > 0)
      expected[index] = val;
    for (size_t i = 0; i < len; ++i)
      out[i] = st.query(i);
    EXPECT_THAT(out, ContainerEq(expected))
              << "Update: " << PrintToString(update) << "\n"
              << "Initial Data: vector of zeros";
  }
}

TEST_P(SegmentTreeTest, RangeUpdateTest) {
  auto len = std::get<1>(GetParam()).size();
  auto updates = std::get<3>(GetParam());
  auto q = [](int a, int b) { return a > b ? a : b; };
  auto u = [](int a, int b) { return b; };
  SegmentTree<int> st(q, len, MemoryOrder::Euler, 0);

  auto out = std::vector<int>(len, 0);
  auto expected = std::vector<int>(len, 0);
  for (const auto &update : updates) {
    // do update
    auto [a, b, val] = update;
    st.update(u, a, b, val);

    // get expected data and data in segment tree
    for (size_t i = a; i <= b; ++i)
      if (len > 0) expected[i] = u(expected[i], val);
    for (size_t i = 0; i < len; ++i)
      out[i] = st.query(i);

    EXPECT_THAT(out, ContainerEq(expected))
              << "rangeUpdate was: " << PrintToString(update);
  }
}

TEST_P(SegmentTreeTest, BuildQueryTest) {
  auto data = std::get<1>(GetParam());
  auto q = [](int a, int b) { return a > b ? a : b; };
  SegmentTree<int> st(q, data, MemoryOrder::Euler);
  auto leaf = std::vector<int>(data.size(), 0);
  if (!data.empty())
    for (size_t i = 0; i < data.size(); ++i)
      leaf[i] = st.query(i);

  EXPECT_THAT(leaf, ContainerEq(data));
}

TEST_P(SegmentTreeTest, RangeQueryTest) {
  auto data = std::get<1>(GetParam());
  auto queries = std::get<5>(GetParam());
  auto q = [](int a, int b) { return a > b ? a : b; };
  SegmentTree<int> st(q, data, MemoryOrder::Euler);

  for (const auto &query : queries) {
    auto [a, b] = query;
    auto result = st.query(a, b);
    int expected = 0;
    for (size_t i = a; i <= b; ++i)
      if (!data.empty()) expected = q(expected, data[i]);
    EXPECT_THAT(result, Eq(expected))
              << "Query: " << PrintToString(query) << "\n"
              << "Initial Data: " << PrintToString(data);
  }
}

}  // end of anonym namespace