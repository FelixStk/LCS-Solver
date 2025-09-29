/*******************************************************************************
 * @file cantor_test.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief GTests for cantor function (as introduction to GTest. Test is trivial.
 ******************************************************************************/

#include "util/Cantor.h"

#include <random>
#include <utility>

#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace {

using ::lcs_solver::util::Cantor;

using ::testing::Eq;
using ::testing::ValuesIn;
using Pair = std::pair<size_t, size_t>;
using TestParam = std::pair<std::string, Pair>;


//== Getters for Testparameter Vector ==========================================
std::vector<TestParam> getHardCodedParams() {
  std::vector<TestParam> params = {
      {"ZeroZero", {0, 0}},
      {"ZeroOne", {0, 1}},
      {"OneZero", {1, 0}}
  };
  return params;
}

std::vector<TestParam> getRandomParams(const unsigned int initSeed,
                                       const size_t nParams,
                                       const Pair numRange) {
  std::vector<TestParam> params;
  std::uniform_int_distribution<size_t> dist(numRange.first, numRange.second);
  for (size_t i = 0; i < nParams; ++i) {
    const size_t seed = initSeed + i;
    std::mt19937 gen(seed);
    TestParam x;
    std::string name = "Rnd_" + std::to_string(seed);
    Pair value = {dist(gen), dist(gen)};
    params.emplace_back(std::move(name), std::move(value));
  }
  return params;
}


//==== Initialize Test Suites ==================================================
class CantorTest : public testing::TestWithParam<TestParam> {};

const auto seed = std::random_device{}();
auto hardCodedParam = getHardCodedParams();
auto rndParam = getRandomParams(
    seed,     // initSeed
    6,        // nParams
    {0, 100}  // range for val in pair
);

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters, CantorTest,
    ValuesIn(hardCodedParam),
    [](const testing::TestParamInfo<CantorTest::ParamType> &info) {
      return info.param.first;
    });

INSTANTIATE_TEST_SUITE_P(
    RandomParameters, CantorTest,
    ValuesIn(rndParam),
    [](const testing::TestParamInfo<CantorTest::ParamType> &info) {
      return info.param.first;
    });

//==== Initialize Test Suites ==================================================

TEST_P(CantorTest, ReturnValue) {
  auto in = GetParam().second;
  auto a = in.first;
  auto b = in.second;
  auto expected = (a+b)*(a+b+1)/2 +b;
  auto out = Cantor::calc(in);
  EXPECT_THAT(out, Eq(expected));
  EXPECT_FALSE(false);
}

} // end of namespace