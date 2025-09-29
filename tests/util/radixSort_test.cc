/*******************************************************************************
 * @file radixSort_test.cc
 * @author Felix Steinkopp
 * @version 1.2
 * @brief GTests for radix sort
 ******************************************************************************/

#include "util/RadixSort.h"

#include <random>
#include <utility>
#include <set>

#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

namespace {

using ::lcs_solver::util::RadixSort;

using ::testing::Eq;
using ::testing::ContainerEq;
using ::testing::ValuesIn;
using ::testing::PrintToString;

using Pair = std::pair<size_t, size_t>; // (Key, Value) pairs in radix sort
using Vec = std::vector<Pair>;
using TestParam = std::pair<std::string, Vec>;

//== Getters for Testparameter Vector ==========================================
std::vector<TestParam> getHardCodedParams() {
  std::vector<TestParam> params = {
      {"EmptyVector", {}}
  };
  return params;
}

std::vector<TestParam> getRandomParams(const unsigned int initSeed,
                                       const size_t nParams = 5,
                                       const size_t len = 10,
                                       const Pair numRange = {0, 100}) {
  std::vector<TestParam> params;
  params.reserve(nParams);
  std::uniform_int_distribution<size_t> dist(numRange.first, numRange.second);
  for (size_t i = 0; i < nParams; ++i) {
    const size_t seed = initSeed + i;
    std::mt19937 gen(seed);
    std::string name = "Rnd_" + std::to_string(seed);
    Vec v;
    std::generate_n(std::back_inserter(v), len, [&]() -> Pair {
      return {dist(gen), dist(gen)};
    });
    params.emplace_back(std::move(name),std::move(v));
  }
  return params;
}

//==== Initialize Test Suites ==================================================
class RadixSortTest : public testing::TestWithParam<TestParam> {};

const auto seed = std::random_device{}();
//const unsigned int seed = 3734371328;
auto hardCodedParam = getHardCodedParams();
auto rndParam = getRandomParams(
    seed,     // initSeed
    6,        // nParams
    10,       // size of vectors to sort
    {0, 100}  // range for random values in (key,value) pairs
);

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters, RadixSortTest,
    ValuesIn(hardCodedParam),
    [](const testing::TestParamInfo<RadixSortTest::ParamType> &info) {
      return info.param.first;
    });

INSTANTIATE_TEST_SUITE_P(
    RandomParameters, RadixSortTest,
    ValuesIn(rndParam),
    [](const testing::TestParamInfo<RadixSortTest::ParamType> &info) {
      return info.param.first;
    });

//==== Initialize Test Suites ==================================================

TEST_P(RadixSortTest, NormalSort) {
  auto &in = GetParam().second;
//  auto span = std::span<const Pair>(in.begin(), in.end());
  auto [sorted, label] = RadixSort<size_t, size_t>::sort(in);
  auto expected = in;
  std::sort(expected.begin(), expected.end(), [](Pair a, Pair b) {
    return a.first < b.first;
  });
  EXPECT_THAT(sorted, ContainerEq(expected));
}

TEST_P(RadixSortTest, LabelCorrectness) {
  auto &in = GetParam().second;
  auto [sorted, label] = RadixSort<size_t, size_t>::sort(in);
  auto expected = in;
  std::sort(expected.begin(), expected.end(), [](Pair a, Pair b) {
    return a.first < b.first;
  });

  for (size_t i = 0; i < in.size(); ++i) {
    EXPECT_THAT(sorted[label[i]], Eq(in[i]));
  }
}

TEST_P(RadixSortTest, HandlesDuplicatesCorrectly) {
  auto &in = GetParam().second;
  auto [sorted, label] = RadixSort<size_t, size_t>::sort(in, true);

  // Erase duplicates with a set (in radixsort erase and unique is used)
  auto comp = [](const Pair &a, const Pair &b) { return a.first < b.first; };
  std::set<Pair, decltype(comp)> s(comp);
  size_t size = in.size();
  for (size_t i = 0; i < size; ++i)
    s.insert(in[i]);
  std::vector<Pair> expected(s.begin(), s.end());

  EXPECT_THAT(sorted, ContainerEq(expected));
  for (size_t i = 0; i < in.size(); ++i) {
    /* note: only compare the keys, because uniqueness of keys is not given */
    EXPECT_THAT(sorted[label[i]].first, Eq(in[i].first))
              << "labels: " << PrintToString(label) << std::endl
              << "sorted: " << PrintToString(sorted) << std::endl
              << "i: " << i;
  }
}

TEST_P(RadixSortTest, Monostate) {
  std::vector<size_t> keys;
  for (const auto &[key, val] : GetParam().second) {
    keys.push_back(key);
  }

  auto [sorted, label] = lcs_solver::util::RadixSort<size_t>::sort(keys);

  auto expected = keys;
  std::sort(expected.begin(), expected.end());
  std::vector<std::pair<size_t, std::monostate>> expected_pairs;
  for (auto e : expected) {
    expected_pairs.push_back({e, {}});
  }

  EXPECT_THAT(sorted, ContainerEq(expected_pairs));
}

} // end of namespace