/*******************************************************************************
 * @file benchmark/bm_segtree.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief Benchmarking of Segment Trees
 ******************************************************************************/
#include <functional>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include "benchmark/benchmark.h"
#include "structures/Matrix.h"
#include "structures/MaxSegTreeND.h"

namespace {
using ::lcs_solver::structures::MaxSegTreeND;
using uint = MaxSegTreeND::uint;
using Matrix = lcs_solver::structures::Matrix<uint>;

/*******************************************************************************
 * Helper function to generate random matrices
 * @param dimensions std::vector<uint> dimensions of the matrix to generate
 * @return lcs_solver::structures::Matrix<uint> a matrix with random values
 ******************************************************************************/
Matrix generateMatrix(const std::vector<uint> &dimensions) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  auto multi = std::accumulate(begin(dimensions), end(dimensions), 1, std::multiplies<int>());

  std::uniform_int_distribution<> dist(0, multi);
  Matrix mat = Matrix(dimensions);
  for (auto &val : mat.data) {
    val = dist(gen);
  }
  return mat;
}

/*******************************************************************************
 * Helper function to generate random indices
 * @param dimensions vector<uint> the area that the tree manages. For example if
 * dimensions = {d0,d1} the tree manages [0:d0-1]x[0:d1-1]
 * @return vector<uint> to access a random leave in the segment tree
 ******************************************************************************/
std::vector<uint> generateRandomIndices(
    const std::vector<uint> &dimensions
) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::vector<uint> indices;
  for (uint const dim : dimensions) {
    std::uniform_int_distribution<> dis(0, dim - 1);
    indices.push_back(dis(gen));
  }
  return indices;
}

/*******************************************************************************
 * Helper function to generate random range
 * @param dimensions vector<uint> the area that the tree manages
 * @return a random subarea of the area that the tree manages
 ******************************************************************************/
std::vector<MaxSegTreeND::Pair> generateRandomRange(
    const std::vector<uint> &dimensions
) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::vector<MaxSegTreeND::Pair> range;
  for (uint const dim : dimensions) {
    std::uniform_int_distribution<> dis(0, dim - 1);
    uint start = dis(gen);
    uint end = dis(gen);
    if (start > end) std::swap(start, end);
    range.emplace_back(start, end);
  }
  return range;
}

/*******************************************************************************
 * BM_PointQuery
 * @param state benchmark::State reference
 ******************************************************************************/
void BM_PointQuery(benchmark::State &state) {
  const auto n = static_cast<uint>(state.range(0));
  const std::vector<uint> dimensions = {n, n}; // 2D tree of size 100x100
  auto mat = generateMatrix(dimensions);
  MaxSegTreeND tree(mat);

  for (auto _ : state) {
    state.PauseTiming();
    auto indices = generateRandomIndices(dimensions);
    state.ResumeTiming();
    benchmark::DoNotOptimize(tree.query(indices));
  }
}

/*******************************************************************************
 * BM_RangeUpdate
 * @param state benchmark::State reference
 ******************************************************************************/
void BM_RangeUpdate(benchmark::State &state) {
  const auto n = static_cast<uint>(state.range(0));
  const std::vector<uint> dimensions = {n, n}; // 2D tree of size nxn
  MaxSegTreeND tree(dimensions);
  uint value = 1;
  for (auto _ : state) {
    state.PauseTiming();
    auto range = generateRandomRange(dimensions);
    state.ResumeTiming();
    tree.update(range, value++);
    benchmark::ClobberMemory();
  }
}

// Definition of Test Range
constexpr uint BM_ST_START = 4;
constexpr uint BM_ST_END = 130;
constexpr uint BM_ST_STEP = 2;

} // end namespace


BENCHMARK(BM_PointQuery)->DenseRange(BM_ST_START, BM_ST_END, BM_ST_STEP);
BENCHMARK(BM_RangeUpdate)->DenseRange(BM_ST_START, BM_ST_END, BM_ST_STEP);

//BENCHMARK(BM_PointQuery)->RangeMultiplier(2)->Range(1 << 2, 1 << 20);
//BENCHMARK(BM_RangeUpdate)->RangeMultiplier(2)->Range(1 << 2, 1 << 20);

/*******************************************************************************
 * @note Use the `--benchmark_format=<console|json|csv>` flag to set the format
 * type. `console` is the default format
 ******************************************************************************/
BENCHMARK_MAIN();
