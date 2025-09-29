/*******************************************************************************
 * @file benchmark/bm_llcs.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief Benchmarking of the LLCS Algorithms
 ******************************************************************************/

#include <random>

#include "algorithms/BaseAlgorithm.h"
#include "algorithms/LLCS/LLCS_STD_FL.h"
#include "algorithms/LLCS/LLCS2_MC.h"
#include "algorithms/LLCS/LLCS2_MC_INC.h"
#include "algorithms/LLCS/LLCS2_MC_INC_E.h"
#include "algorithms/LLCS/LLCS2_MC_1C.h"
#include "algorithms/LLCS/LLCS2_MC_O1_SYNC.h"
#include "algorithms/LLCS/LLCS2_SR_MQ.h"
#include "algorithms/LLCS/LLCS2_SR_RMQ.h"
#include "algorithms/LLCS/LLCS2_SA_MQ.h"
#include "algorithms/LLCS/LLCS2_SA_RMQ.h"

#include "benchmark/benchmark.h"

#include "util/AlgoTestParam.h"
#include "util/ParamGenerator.h"

using ::lcs_solver::algorithms::llcs::LLCS_STD_FL;
using ::lcs_solver::algorithms::llcs::LLCS2_MC;
using ::lcs_solver::algorithms::llcs::LLCS2_MC_INC;
using ::lcs_solver::algorithms::llcs::LLCS2_MC_INC_E;
using ::lcs_solver::algorithms::llcs::LLCS2_MC_1C;
using ::lcs_solver::algorithms::llcs::LLCS2_MC_O1_SYNC;
using ::lcs_solver::algorithms::llcs::LLCS2_SR_MQ;
using ::lcs_solver::algorithms::llcs::LLCS2_SR_RMQ;
using ::lcs_solver::algorithms::llcs::LLCS2_SA_MQ;
using ::lcs_solver::algorithms::llcs::LLCS2_SA_RMQ;

using ::lcs_solver::util::ParamGenerator;
using ::lcs_solver::algorithms::BaseSolution;
using ConstraintType = ::lcs_solver::constraints::ConstraintType;

/*******************************************************************************
 * Benchmarks LLCS2_STD_FL
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_STD_FL(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::Empty,
        seed, // Seed used to init random number generators
        {{m, m}, {n, n}}, // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'} // Alphabet for strings
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS_STD_FL(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}

/*******************************************************************************
 * Benchmarks LLCS2_MC
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_MC(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::MC,
        seed, // Seed used to init random number generators
        {{m, m}, {n, n}}, // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'} // Alphabet for strings
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS2_MC(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}

/*******************************************************************************
 * Benchmarks LLCS2_MC_INC
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_MC_INC(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::MC_INC,
        seed,               // Seed used to init random number generators
        {{m, m}, {n, n}},   // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'}     // c[0] common symbol, elem in the tail are unique
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS2_MC_INC(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}

/*******************************************************************************
 * Benchmarks LLCS2_MC_INC_E
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_MC_INC_E(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::MC_INC,
        seed,               // Seed used to init random number generators
        {{m, m}, {n, n}},   // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'}     // c[0] common symbol, elem in the tail are unique
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS2_MC_INC_E(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}

/*******************************************************************************
 * Benchmarks LLCS2_MC_1C
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_MC_1C(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::MC_1C,
        seed,               // Seed used to init random number generators
        {{m, m}, {n, n}},   // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'}     // c[0] common symbol, elem in the tail are unique
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS2_MC_1C(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}

/*******************************************************************************
 * Benchmarks LLCS2_MC_1C
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_MC_O1_SYNC(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::MC_O1C_SYNC,
        seed,               // Seed used to init random number generators
        {{m, m}, {n, n}},   // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'}     // c[0] common symbol, elem in the tail are unique
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS2_MC_O1_SYNC(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}

/*******************************************************************************
 * Benchmarks LLCS2_SR_MQ
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_SR_MQ(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::SIGMA_R,
        seed,               // Seed used to init random number generators
        {{m, m}, {n, n}},   // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'}     // c[0] common symbol, elem in the tail are unique
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS2_SR_MQ(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}

/*******************************************************************************
 * Benchmarks LLCS2_SR_RMQ
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_SR_RMQ(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::SIGMA_R,
        seed,               // Seed used to init random number generators
        {{m, m}, {n, n}},   // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'}     // c[0] common symbol, elem in the tail are unique
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS2_SR_RMQ(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}


/*******************************************************************************
 * Benchmarks LLCS2_SR_MQ
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_SA_MQ(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::SIGMA,
        seed,               // Seed used to init random number generators
        {{m, m}, {n, n}},   // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'}     // c[0] common symbol, elem in the tail are unique
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS2_SA_MQ(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}

/*******************************************************************************
 * Benchmarks LLCS2_SR_RMQ
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_SA_RMQ(::benchmark::State &state) {
  const auto m = static_cast<size_t>(state.range(0));
  const auto n = static_cast<size_t>(state.range(1));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithoutSol(
        ConstraintType::SIGMA,
        seed,               // Seed used to init random number generators
        {{m, m}, {n, n}},   // Bounds on string lengths
        {0, std::min(m, n)}, // Bound for gaps
        {'x', 'a', 'b'}     // c[0] common symbol, elem in the tail are unique
    );

    const auto &spv = param.pointers;
    auto &map = param.map;
    benchmark::ClobberMemory(); // Don't use caching

    //=== Start Clock ==========================================================
    state.ResumeTiming();
    auto algo = new LLCS2_SA_RMQ(spv, map);
    std::unique_ptr<BaseSolution> sol;
    benchmark::DoNotOptimize(sol = algo->query());

    //=== Tear Down (Stop Clock) ===============================================
    state.PauseTiming();
    delete algo;
    state.ResumeTiming();
  }
}

/*******************************************************************************
 * CustomArguments sets argument ranges for the benchmark.
 * @param b Pointer to the Benchmark instance where the arguments are to be set.
 ******************************************************************************/
static void CustomArguments(benchmark::internal::Benchmark *b) {
  for (int i = 10; i <= 10; i += 10)
    for (int j = 8; j <= 128; j *= 2)
      b->Args({i, j});
}

static void CustomArguments2(benchmark::internal::Benchmark *b) {
  for (int i = 8; i <= 256; i *= 2)
    for (int j = 8; j <= 256; j *= 2)
      b->Args({i, j});
}

BENCHMARK(BM_LLCS2_STD_FL)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_MC)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_MC_INC)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_MC_INC_E)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_MC_1C)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_MC_O1_SYNC)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_SR_MQ)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_SR_RMQ)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_SA_MQ)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_SA_RMQ)->Apply(CustomArguments2);

BENCHMARK_MAIN();
