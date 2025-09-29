/*******************************************************************************
 * @file rmq_test.cc
 * @author Felix Steinkopp
 * @version 1.2
 * @brief GTests for the RMQ based on LCA
 ******************************************************************************/

#include "util/RMQ.h"

#include <random>
#include <utility>
#include <set>

#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace {

using ::lcs_solver::util::RMQ_ONN;
using ::lcs_solver::util::RMQ_nlogn;
using ::lcs_solver::util::RMQ_pm1;
using ::lcs_solver::util::RMQ_ON;
using ::lcs_solver::util::RMQ_TYPE;
using ::lcs_solver::util::CartesianTree;

using ::testing::Eq;
using ::testing::Gt;
using ::testing::ValuesIn;
using ::testing::PrintToString;
using Pair = std::pair<size_t, size_t>; // for specifying data ranges in vectors

//== Parameter structure for the tests =========================================
using TestParam = std::tuple<
    std::string,       ///< Name of test Parameter Combination
    std::vector<int>,  ///< data vector
    std::vector<Pair>  ///< multiple ranges for range max/min queries
>;

//== Getters for Testparameter Vector ==========================================
std::vector<TestParam> getHardCodedParams1() { // for RMQ without PM1 property
  std::vector<TestParam> params = {
      {"EmptyVector", {}, {}},
      {"Len1Vector", {0}, {{0, 0}}},
      {"Len2Vector", {0, 1}, {{0, 1}}},
      {"Len3Vector", {0, 1, 2}, {{0, 2}}},
      {"Len4Vector", {0, 1, 2, 3}, {{0, 3}}},
      {"BadQueryOrder", {0, 1, 2, 3, 4, 5}, {{5, 1}}},
      {"BadQueryIndex", {0, 1, 2, 3, 4, 5}, {{6, 10}}}
  };
  return params;
}

std::vector<TestParam> getHardCodedParams2() { // for RMQ with PM1 property
  std::vector<TestParam> params = {
      {"EmptyVector", {}, {}},
      {"Len1Vector", {0}, {{0, 0}}},
      {"Len2Vector", {0, 1}, {{0, 1}}},
      {"Len3Vector", {0, 1, 2}, {{0, 2}}},
      {"Len4Vector", {0, 1, 2, 3}, {{0, 3}}},
      {"BadQueryOrder", {0, 1, 2, 3, 4, 5}, {{5, 1}}},
      {"BadQueryIndex", {0, 1, 2, 3, 4, 5}, {{6, 10}}}
  };
  return params;
}

/*******************************************************************************
 * getRandomVector generates random data vectors used for testing
 * @param initSeed Seed for random number generator
 * @param nParams Number of TestParam structs to generate
 * @param len bounds for the length of the data vector
 * @param r bounds for the values of the data vector
 * @param nQueries number of queries to generate
 * @param pm1Property whether the data vector satisfies the pm1 property
 * @param pOneStep if pm1Property, probability(data[i+1]-data[i]==1)
 * @return vector of TestParam (param name, data vector and vector of queries)
 ******************************************************************************/
std::vector<TestParam> getRandomVector(const unsigned int initSeed,
                                       const size_t nParams = 5,
                                       const std::pair<int, int> len = {8, 10},
                                       const std::pair<int, int> r = {0, 100},
                                       const size_t nQueries = 5,
                                       const bool pm1Property = false,
                                       const double pOneStep = 0.5) {
  std::vector<TestParam> params;
  std::uniform_int_distribution<int> d0(len.first, len.second);
  std::uniform_int_distribution<int> d1(r.first, r.second);
  std::bernoulli_distribution d2(pOneStep); // pOneStep = p(true)

  for (size_t i = 0; i < nParams; ++i) {
    const size_t seed = initSeed + i;
    std::mt19937 gen(seed);
    const size_t length = d0(gen);
    std::uniform_int_distribution<size_t> d3(0, length - 1);
    TestParam x;

    // Generate Name
    std::get<0>(x) = "Rnd_" + std::to_string(seed);

    // Generate Vector
    if (!pm1Property) {
      std::generate_n(std::back_inserter(std::get<1>(x)), length, [&]() -> int {
        return d1(gen);
      });
    } else {
      std::get<1>(x).reserve(length);
      int n = static_cast<int>(length);
      for (int j = 0, val = d1(gen); j < n; ++j, val += d2(gen) ? 1 : -1) {
        std::get<1>(x).emplace_back(val);
      }
    }
    std::generate_n(std::back_inserter(std::get<2>(x)), nQueries, [&]() -> Pair {
      auto a = d3(gen);
      auto b = d3(gen);
      return a < b ? Pair(a, b) : Pair(b, a);
    });
    params.push_back(std::move(x));
  }
  return params;
}

//==== Initialize Test Suites ==================================================
class RMQTest : public testing::TestWithParam<TestParam> {};
class RMQPM1Test : public testing::TestWithParam<TestParam> {};

const auto seed = std::random_device{}();
// const unsigned int seed = 2275096221;
auto hardCodedParam1 = getHardCodedParams1();
auto rndVector = getRandomVector(
    seed,     // initSeed
    10,        // nParams
    {1, 4},   // min and max of elements in the data vector to sort
    {0, 10}   // range for random values in the vector
);
auto rndVectorBig = getRandomVector(
    seed,          // initSeed
    5,             // nParams
    {1023, 1025},  // min and max of elements in the data vector to sort
    {0, 90}        // range for random values in the vector
);

auto hardCodedParam2 = getHardCodedParams2();
auto rndPM1Vector = getRandomVector(
    seed,     // initSeed
    5,       // nParams
    {1, 4},   // min and max of elements in the data vector to sort
    {0, 10},  // range for the head value of the data vector
    4,        // Number of random queries
    true,     // vector should have pm1-property: |vec[i+1]-vec[i]| == 1
    0.5       // probability for vec[i+1] - vec[i] = +1
);
auto rndPM1VectorBig = getRandomVector(
    seed,          // initSeed
    5,             // nParams
    {1023, 1025},  // min and max of elements in the data vector to sort
    {0, 80},       // range for the head value of the data vector
    4,             // Number of random queries
    true,          // vector should have pm1-property: |vec[i+1]-vec[i]| == 1
    0.50001        // probability for vec[i+1] - vec[i] = +1
);

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters, RMQTest,
    ValuesIn(hardCodedParam1),
    [](const testing::TestParamInfo<RMQTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    RandomParameters, RMQTest,
    ValuesIn(rndVector),
    [](const testing::TestParamInfo<RMQTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    RandomParametersBig, RMQTest,
    ValuesIn(rndVectorBig),
    [](const testing::TestParamInfo<RMQTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters, RMQPM1Test,
    ValuesIn(hardCodedParam2),
    [](const testing::TestParamInfo<RMQTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    RandomParameters, RMQPM1Test,
    ValuesIn(rndPM1Vector),
    [](const testing::TestParamInfo<RMQTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    RandomParametersBig, RMQPM1Test,
    ValuesIn(rndPM1VectorBig),
    [](const testing::TestParamInfo<RMQTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

//==== Initialize Test Suites ==================================================

TEST_P(RMQTest, ONN_Max) {
  const auto &v = std::get<1>(GetParam());
  const auto &queries = std::get<2>(GetParam());
  auto rmq = new RMQ_ONN<int>(RMQ_TYPE::MAX, v);

  for (auto [a, b] : queries) {
    auto out = rmq->query(a, b);
    bool badQuery = !(a < v.size() && b < v.size());
    if (!badQuery) {
      if (b < a) std::swap(a, b);
      auto start = std::next(v.begin(), static_cast<int>(a));
      auto end = std::next(v.begin(), static_cast<int>(b + 1));
      auto expected = std::max_element(start, end);
      auto pos = static_cast<size_t>(std::distance(v.begin(), expected));
      EXPECT_THAT(v[out], Eq(*expected));
      EXPECT_THAT(out, Eq(pos));
    } else
      EXPECT_THAT(out, Gt(v.size()));
  }
  delete rmq;
}

TEST_P(RMQTest, ONN_Min) {
  const auto &v = std::get<1>(GetParam());
  const auto &queries = std::get<2>(GetParam());
  auto rmq = new RMQ_ONN<int>(RMQ_TYPE::MIN, v);

  for (auto [a, b] : queries) {
    auto out = rmq->query(a, b);
    bool badQuery = !(a < v.size() && b < v.size());
    if (!badQuery) {
      if (b < a) std::swap(a, b);
      auto start = std::next(v.begin(), static_cast<int>(a));
      auto end = std::next(v.begin(), static_cast<int>(b + 1));
      auto expected = std::min_element(start, end);
      auto pos = static_cast<size_t>(std::distance(v.begin(), expected));
      EXPECT_THAT(v[out], Eq(*expected));
      EXPECT_THAT(out, Eq(pos));
    } else
      EXPECT_THAT(out, Gt(v.size()));
  }
  delete rmq;
}
//
TEST_P(RMQTest, NLOGN_Max) {
  const auto &v = std::get<1>(GetParam());
  const auto &queries = std::get<2>(GetParam());
  auto rmq = new RMQ_nlogn<int>(RMQ_TYPE::MAX, v);

  for (auto [a, b] : queries) {
    auto out = rmq->query(a, b);
    bool badQuery = !(a < v.size() && b < v.size());
    if (!badQuery) {
      if (b < a) std::swap(a, b);
      auto start = std::next(v.begin(), static_cast<int>(a));
      auto end = std::next(v.begin(), static_cast<int>(b + 1));
      auto expected = std::max_element(start, end);
      auto pos = static_cast<size_t>(std::distance(v.begin(), expected));
      EXPECT_THAT(v[out], Eq(*expected));
      EXPECT_THAT(out, Eq(pos));
    } else
      EXPECT_THAT(out, Gt(v.size()));
  }
  delete rmq;
}

TEST_P(RMQTest, NLOGN_Min) {
  const auto &v = std::get<1>(GetParam());
  const auto &queries = std::get<2>(GetParam());
  auto rmq = new RMQ_nlogn<int>(RMQ_TYPE::MIN, v);

  for (auto [a, b] : queries) {
    auto out = rmq->query(a, b);
    bool badQuery = !(a < v.size() && b < v.size());
    if (!badQuery) {
      if (b < a) std::swap(a, b);
      auto start = std::next(v.begin(), static_cast<int>(a));
      auto end = std::next(v.begin(), static_cast<int>(b + 1));
      auto expected = std::min_element(start, end);
      auto pos = static_cast<size_t>(std::distance(v.begin(), expected));
      EXPECT_THAT(v[out], Eq(*expected));
      EXPECT_THAT(out, Eq(pos));
    } else
      EXPECT_THAT(out, Gt(v.size()));
  }
  delete rmq;
}

TEST_P(RMQPM1Test, RMQPM1_Max) {
  const auto &v = std::get<1>(GetParam());
  const auto &queries = std::get<2>(GetParam());
  auto rmq = new RMQ_pm1<int>(RMQ_TYPE::MAX, v);

  for (auto [a, b] : queries) {
    auto out = rmq->query(a, b);
    bool badQuery = !(a < v.size() && b < v.size());
    if (!badQuery) {
      if (b < a) std::swap(a, b);
      auto start = std::next(v.begin(), static_cast<int>(a));
      auto end = std::next(v.begin(), static_cast<int>(b + 1));
      auto expected = std::max_element(start, end);
      auto pos = static_cast<size_t>(std::distance(v.begin(), expected));
      EXPECT_THAT(v[out], Eq(*expected));
      EXPECT_THAT(out, Eq(pos)) << "query: (" << a << " " << b << ")\n";
    } else
      EXPECT_THAT(out, Gt(v.size()));
  }
  delete rmq;
}

TEST_P(RMQPM1Test, RMQPM1_Min) {
  const auto &v = std::get<1>(GetParam());
  const auto &queries = std::get<2>(GetParam());
  auto rmq = new RMQ_pm1<int>(RMQ_TYPE::MIN, v);

  for (auto [a, b] : queries) {
    auto out = rmq->query(a, b);
    bool badQuery = !(a < v.size() && b < v.size());
    if (!badQuery) {
      if (b < a) std::swap(a, b);
      auto start = std::next(v.begin(), static_cast<int>(a));
      auto end = std::next(v.begin(), static_cast<int>(b + 1));
      auto expected = std::min_element(start, end);
      auto pos = static_cast<size_t>(std::distance(v.begin(), expected));
      EXPECT_THAT(v[out], Eq(*expected));
      EXPECT_THAT(out, Eq(pos));
    } else
      EXPECT_THAT(out, Gt(v.size()));
  }
  delete rmq;
}

TEST_P(RMQTest, ON_Max) {
  const auto &v = std::get<1>(GetParam());
  const auto &queries = std::get<2>(GetParam());
  auto rmq = new RMQ_ON<int>(RMQ_TYPE::MAX, v);
//  auto rmq = std::make_unique<CartesianTree<int>>(RMQ_TYPE::MAX, v);  // Use std::make_unique

  for (auto [a, b] : queries) {
    auto out = rmq->query(a, b);
    bool badQuery = !(a < v.size() && b < v.size());
    if (!badQuery) {
      if (b < a) std::swap(a, b);
      auto start = std::next(v.begin(), static_cast<int>(a));
      auto end = std::next(v.begin(), static_cast<int>(b + 1));
      auto expected = std::max_element(start, end);
      auto pos = static_cast<size_t>(std::distance(v.begin(), expected));
      EXPECT_THAT(v[out], Eq(*expected));
      EXPECT_THAT(out, Eq(pos));
    } else
      EXPECT_THAT(out, Gt(v.size()));
  }
  delete rmq;
}

TEST_P(RMQTest, ON_Min) {
  const auto &v = std::get<1>(GetParam());
  const auto &queries = std::get<2>(GetParam());
  auto rmq = new RMQ_ON<int>(RMQ_TYPE::MIN, v);

  for (auto [a, b] : queries) {
    auto out = rmq->query(a, b);
    bool badQuery = !(a < v.size() && b < v.size());
    if (!badQuery) {
      if (b < a) std::swap(a, b);
      auto start = std::next(v.begin(), static_cast<int>(a));
      auto end = std::next(v.begin(), static_cast<int>(b + 1));
      auto expected = std::min_element(start, end);
      auto pos = static_cast<size_t>(std::distance(v.begin(), expected));
      EXPECT_THAT(v[out], Eq(*expected)) << "query(" << a << ", " << b << ")";
      EXPECT_THAT(out, Eq(pos)) << "query(" << a << ", " << b << ")";
    } else
      EXPECT_THAT(out, Gt(v.size()));
  }
  delete rmq;
}

} // end of namespace