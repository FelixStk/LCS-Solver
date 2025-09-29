/*******************************************************************************
 * @file benchmark/bm_rmq.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief Benchmarking of the RMQ Datastructures (Setup and Query)
 ******************************************************************************/

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <iterator>
#include <random>
#include <utility>
#include <vector>

#include "benchmark/benchmark.h"
#include "util/RMQ.h"

namespace {
using ::lcs_solver::util::RMQ_ONN;
using ::lcs_solver::util::RMQ_nlogn;
using ::lcs_solver::util::RMQ_pm1;
using ::lcs_solver::util::RMQ_ON;
using ::lcs_solver::util::RMQ_TYPE;
using uint = unsigned int;

/*******************************************************************************
 * @brief Generates a vector of random numbers
 * @param seed uint used to initialize the random number generator
 * @param length unit size of the vector to generate
 * @param r std::pair<int,int> determines the range of possible numbers
 * @param pm1Property bool if the generate vector must have the pm1 property
 * @param pOneStep probability for a step up in case `pm1Property==true`
 * @return std::vector<int> of random numbers in `r`
 ******************************************************************************/
std::vector<int> genRandomVector(
    const uint seed,
    const uint length,
    const std::pair<int, int> r = {0, 100},
    const bool pm1Property = false,
    const double pOneStep = 0.5
) {
  std::uniform_int_distribution<int> d1(r.first, r.second);
  std::bernoulli_distribution d2(pOneStep);
  std::uniform_int_distribution<size_t> const d3(0, length - 1);
  std::mt19937 gen(seed);

  std::vector<int> vec;
  vec.reserve(length);
  if (pm1Property) {
    int val = d1(gen);
    int const n = static_cast<int>(length);
    for (int j = 0; j < n; ++j) {
      vec.emplace_back(val);
      const bool oneStep = d2(gen);
      if (oneStep) {
        val = val + (val < r.second ? 1 : 0);
      } else {
        val = val - (val > r.first ? 1 : 0);
      }
    } // end for
  } else {
    // PM1-Property is not active
    std::generate_n(std::back_inserter(vec), length, [&]() -> int {
      return d1(gen);
    });
  }
  return vec;
}

/*******************************************************************************
 * @brief Generates a random Range for a RMQ Query
 * @param seed to init random number generators
 * @param length the size of the data vector
 * @param random weather to generate a random pair in [0, length -1]**2 or
 * otherwise return {0,length -1}
 * @return std::pair<size_t, size_t> in [0, length -1]**2
 ******************************************************************************/
std::pair<size_t, size_t> genQueryPair(
    const uint seed,
    const size_t length,
    const bool random = false
) {
  std::pair<size_t, size_t> res;
  if (random) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> dist(0, length - 1);
    size_t a = dist(gen);
    size_t b = dist(gen);
    if (a > b) std::swap(a, b);
    res = {a, b};
  } else {
    res = {0, length - 1};
  }
  return res;
}

//int64_t LOG2(size_t n) {
//  auto d = static_cast<double>(n);
//  return static_cast<int64_t>(log2(d));
//}

void BM_SETUP_RMQ_ONN(::benchmark::State &state) {
  size_t seed = std::random_device{}();
  const auto n = state.range(0);
  for (auto _ : state) {
    state.PauseTiming();
    seed++;
    auto vec = genRandomVector(seed, n);
    state.ResumeTiming();
    auto *rmq = new RMQ_ONN<int>(RMQ_TYPE::MAX, vec);
    state.PauseTiming();
    delete rmq;
    state.ResumeTiming();
  }
}

void BM_SETUP_RMQ_nlogn(::benchmark::State &state) {
  size_t seed = std::random_device{}();
  const auto n = state.range(0);
  for (auto _ : state) {
    state.PauseTiming();
    seed++;
    auto vec = genRandomVector(seed, n);
    state.ResumeTiming();
    auto rmq = new RMQ_nlogn<int>(RMQ_TYPE::MAX, vec);
    state.PauseTiming();
    delete rmq;
    state.ResumeTiming();
  }
}

void BM_SETUP_RMQ_pm1(::benchmark::State &state) {

  size_t seed = std::random_device{}();
  const auto n = state.range(0);
  for (auto _ : state) {
    state.PauseTiming();
    seed++;
    auto vec = genRandomVector(seed, n, {0, 100}, true, 0.50);
    state.ResumeTiming();
    auto rmq = new RMQ_pm1<int>(RMQ_TYPE::MAX, vec);
    state.PauseTiming();
    delete rmq;
    state.ResumeTiming();
  }
}

void BM_SETUP_RMQ_ON(::benchmark::State &state) {
  size_t seed = std::random_device{}();
  const auto n = state.range(0);
  for (auto _ : state) {
    state.PauseTiming();
    seed++;
    auto vec = genRandomVector(seed, n);
    state.ResumeTiming();
    auto *rmq = new RMQ_ONN<int>(RMQ_TYPE::MAX, vec);
    state.PauseTiming();
    delete rmq;
    state.ResumeTiming();
  }
}

void BM_QUERY_RMQ_ONN(::benchmark::State &state) {
  const auto n = state.range(0);
  size_t seed = std::random_device{}();
  auto vec = genRandomVector(seed, n);
  auto *rmq = new RMQ_ONN<int>(RMQ_TYPE::MAX, vec);
  for (auto _ : state) {
    state.PauseTiming();
    const auto r = genQueryPair(seed++, n);
    state.ResumeTiming();
    rmq->query(r.first, r.second);
  }
  delete rmq;
}

void BM_QUERY_RMQ_nlogn(::benchmark::State &state) {
  const auto n = state.range(0);
  size_t seed = std::random_device{}();
  auto vec = genRandomVector(seed, n);
  auto rmq = new RMQ_nlogn<int>(RMQ_TYPE::MAX, vec);
  for (auto _ : state) {
    state.PauseTiming();
    const auto r = genQueryPair(seed++, n);
    state.ResumeTiming();
    rmq->query(r.first, r.second);
  }
  delete rmq;
}

void BM_QUERY_RMQ_pm1(::benchmark::State &state) {
  const auto n = state.range(0);
  size_t seed = std::random_device{}();
  auto vec = genRandomVector(seed, n, {0, 100}, true, 0.50);
  auto rmq = new RMQ_pm1<int>(RMQ_TYPE::MAX, vec);
  for (auto _ : state) {
    state.PauseTiming();
    const auto r = genQueryPair(seed++, n);
    state.ResumeTiming();
    rmq->query(r.first, r.second);
  }
  delete rmq;
}

void BM_QUERY_RMQ_ON(::benchmark::State &state) {
  const auto n = state.range(0);
  size_t seed = std::random_device{}();
  auto vec = genRandomVector(seed, n);
  auto rmq = new RMQ_ON<int>(RMQ_TYPE::MAX, vec);
  for (auto _ : state) {
    state.PauseTiming();
    const auto r = genQueryPair(seed++, n);
    state.ResumeTiming();
    rmq->query(r.first, r.second);
  }
  delete rmq;
}

void BM_QUERY_NO_RMQ(::benchmark::State &state) {
  const auto n = state.range(0);
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    // Setup of Benchmark
    state.PauseTiming();
    seed++;
    auto vec = genRandomVector(seed, n);
    //const auto r = genQueryPair(seed, n);
    const auto r = std::pair<int, int>(0, n - 1);
    auto start = std::next(vec.begin(), static_cast<int>(r.first));
    auto end = std::next(vec.begin(), static_cast<int>(r.second + 1));
    benchmark::ClobberMemory(); // Don't use caching
    state.ResumeTiming();
    // Function to Benchmark
    int res = 0;
    benchmark::DoNotOptimize(res = *std::max_element(start, end));
  }
}

constexpr uint START = 100;
constexpr uint END = 10000;
constexpr uint STEP = 500;

} // end namespace

BENCHMARK(BM_SETUP_RMQ_ONN)->DenseRange(START,END,STEP);
BENCHMARK(BM_SETUP_RMQ_nlogn)->DenseRange(START,END,STEP);
BENCHMARK(BM_SETUP_RMQ_pm1)->DenseRange(START,END,STEP);
BENCHMARK(BM_SETUP_RMQ_ON)->DenseRange(START,END,STEP);

BENCHMARK(BM_QUERY_RMQ_ONN)->DenseRange(START,END,STEP);
BENCHMARK(BM_QUERY_RMQ_nlogn)->DenseRange(START,END,STEP);
BENCHMARK(BM_QUERY_RMQ_pm1)->DenseRange(START,END,STEP);
BENCHMARK(BM_QUERY_RMQ_ON)->DenseRange(START,END,STEP);

BENCHMARK(BM_QUERY_NO_RMQ)->DenseRange(START,END,STEP);

/*******************************************************************************
 * @note Use the `--benchmark_format=<console|json|csv>` flag to set the format
 * type. `console` is the default format
 ******************************************************************************/
BENCHMARK_MAIN();
