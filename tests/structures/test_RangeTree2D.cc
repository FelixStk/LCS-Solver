/*******************************************************************************
 * @file test_RangeTree2D.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief GTests for checking the correctness of the RangeTree2D struct
 ******************************************************************************/

#include <algorithm>
#include <random>
#include <string>
#include <vector>


#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

#include "structures/RangeTree2D.h"

namespace {

using ::lcs_solver::structures::RangeTree2D;

using ::testing::Lt;
using ::testing::Le;
using ::testing::Eq;
using ::testing::ContainerEq;
using ::testing::PrintToString;
using ::testing::ValuesIn;

using Pair = std::pair<int, int>;
using Points = std::vector<Pair>;
using Ranges = std::vector<Pair>;

//== Parameter structure for the tests =========================================
using TestParam = std::tuple<
    std::string,                            // 0: name
    std::vector<Pair>,                    // 1: points
    std::vector<Ranges>,                    // 2: queries
    std::vector<Points>                     // 3: solutions for queries
>;

/*******************************************************************************
 * getHardCodedParams
 * @return std::vector<TestParam> with some important edge cases
 ******************************************************************************/
std::vector<TestParam> getHardCodedParams() {
  std::vector<TestParam> params;
  params.emplace_back(
      "QueryEmptyTree",
      std::vector<Pair>{},
      std::vector<Ranges>{
          {{-3,3},{3,3}},     // Query 1
          {{0,0},{0,0}},      // Query 2
      },
      std::vector<Points>{
          {},                 // Solution 1
          {}                  // Solution 2
      }
  );
  return params;
}

/*******************************************************************************
 * getHardCodedParams
 * @param initSeed unsigned int used for random number generation
 * @param nParams number of queries and length of the resulting TestParam-Vector
 * @return std::vector<TestParam> with some important edge cases
 ******************************************************************************/
std::vector<TestParam> genTextBookExamples(
  const unsigned int initSeed,
  const size_t nParams
) {
  std::vector<TestParam> params;
  for (size_t numParam = 0; numParam < nParams; ++numParam) {
    const size_t seed = initSeed + numParam;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> dist(0, 99);
    TestParam x;
    auto &[name, points, queries, solutions] = x;

    name = "Rnd_" + std::to_string(seed);
    points = std::vector<Pair>{
      {2,19}, {7,10}, {12, 3}, {17,62}, {21,49}, {41,95}, {58,59}, {93,70},
      {5,80}, {8,37}, {15,99}, {33,30}, {52,23}, {67,89}
    };

    Ranges r = {{dist(gen),dist(gen)},{dist(gen),dist(gen)}};
    if(r[0].first > r[0].second) { std::swap(r[0].first, r[0].second); }
    if(r[1].first > r[1].second) { std::swap(r[1].first, r[1].second); }
    queries.push_back(r);

    // std::vector<Points> answer;
    solutions.emplace_back();
    for(auto &p : points) {
      if(r[0].first <= p.first && p.first <= r[0].second) {
        if(r[1].first <= p.second && p.second <= r[1].second) {
          solutions[0].push_back(p);
        }
      }
    }
    params.push_back(x);
  }
  return params;
}

/*******************************************************************************
 * genTestParams
 * @param initSeed unsigned int used for random number generation
 * @param nPointClouds size_t for the number of different point clouds to create
 * @param nPointsRange pair specifying the min/max of points in a point clouds
 * @param nQueries size_t specifying the number of test queries for each cloud
 * @param numRange Pair specifying a closed interval to draw random numbers
 * @return std::vector<TestParam>
 ******************************************************************************/
std::vector<TestParam> genTestParams(
    const unsigned int initSeed,
    const size_t nPointClouds,
    const std::pair<size_t,size_t> nPointsRange,
    const size_t nQueries,
    const Pair   numRange
 ) {
  std::vector<TestParam> params;
  for (size_t numParam = 0; numParam < nPointClouds; ++numParam) {
    // Setup of generators, constants and empty TestParam x
    const size_t seed = initSeed + numParam;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> dist(numRange.first, numRange.second);
    TestParam x;
    auto &[name, points, queries, solutions] = x;

    // Generate Parameter: Name
    name = "Rnd_" + std::to_string(seed);

    // Generate Parameter: Points
    std::uniform_int_distribution<size_t> dist2(nPointsRange.first, nPointsRange.second);
    const size_t nPoints = dist2(gen);
    std::generate_n(std::back_inserter(points), nPoints,
                    [&]() -> std::pair<int, int> {
                      return {dist(gen), dist(gen)};
                    });

    // Generate Parameter: Queries
    std::generate_n(std::back_inserter(queries), nQueries,
                    [&]() -> Ranges {
                      Ranges r = {{dist(gen),dist(gen)},{dist(gen),dist(gen)}};
                      if(r[0].first > r[0].second) { std::swap(r[0].first, r[0].second); }
                      if(r[1].first > r[1].second) { std::swap(r[1].first, r[1].second); }
                      return r;
                    });

    // Generate Parameter: Solutions
    solutions.resize(nQueries);
    for(size_t i = 0; i < nQueries; ++i) {
      const auto & r = queries[i];
      const auto & [lowerX, upperX] = r[0];
      const auto & [lowerY, upperY] = r[1];
      for(auto &p : points) {
        const int x = p.first;
        const int y = p.second;
        if(lowerX <= x && x <= upperX && lowerY <= y && y <= upperY) {
            solutions[i].push_back(p);
        }
      }
    }

    params.push_back(x);
  } // for each PointCloud
  return params;
}

//==== Initialize Test Suites ==================================================

class RangeTree2DTest : public testing::TestWithParam<TestParam> {};

const auto seed = std::random_device{}();
// const auto seed = 838653516;
// const auto seed = 3688269970;
// const auto seed = 3548924442;

std::vector<TestParam> hardCodedParam = getHardCodedParams();
std::vector<TestParam> textBookExamples = genTextBookExamples(seed, 5);
std::vector<TestParam> randomParams = genTestParams(
    seed,     // seed
    5,        // nPointClouds,
    {5,10},       // nPoints,
    5,        // nQueries,
    {1, 100}  // numRange
 );

// INSTANTIATE_TEST_SUITE_P(
//     HardCodedParameters, RangeTree2DTest,
//     ValuesIn(hardCodedParam),
//     [](const testing::TestParamInfo<RangeTree2DTest::ParamType> &info) {
//       return std::get<0>(info.param);
//     });

// INSTANTIATE_TEST_SUITE_P(
//     TexBookExample, RangeTree2DTest,
//     ValuesIn(textBookExamples),
//     [](const testing::TestParamInfo<RangeTree2DTest::ParamType> &info) {
//       return std::get<0>(info.param);
//     });

INSTANTIATE_TEST_SUITE_P(
    RandomParameters, RangeTree2DTest,
    ValuesIn(randomParams),
    [](const testing::TestParamInfo<RangeTree2DTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

//==== Parametrized Test =======================================================
TEST_P(RangeTree2DTest, Query) {
  auto name = std::get<0>(GetParam());
  auto points = std::get<1>(GetParam());
  auto queries = std::get<2>(GetParam());
  auto expected = std::get<3>(GetParam());

  auto tree = RangeTree2D<int,int>(points);
  int i = 0;
  for(const auto &r : queries) {
    auto result = tree.query(r[0], r[1]);
    Points solution;
    for (const auto &point : result) {
      solution.emplace_back(point);
    }
    std::ranges::sort(solution);
    std::ranges::sort(expected[i]);
    EXPECT_THAT(solution, ContainerEq(expected[i]))
      << "Debug Info (" << name << ", i == "<< i << " )\n"
      << "Points: " << PrintToString(points) << "\n"
      << "Ranges: " << PrintToString(r) << "\n"
      << "Tree: \n" << tree.DebugString() << "\n";
    i++;
  }
}


}  // end of anonymous namespace