/*******************************************************************************
 * @file bm_llcs_mc.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief Benchmarking for showing k dependency in the runtimes of LLCS2_MC and
 *        LLCS2_MC_O1_SYNC
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
#include "constraints/local/Constraint_MC_O1C_SYNC.h"
#include "constraints/local/Constraint_MC_1C.h"

#include "benchmark/benchmark.h"

#include "util/AlgoTestParam.h"
#include "util/ParamGenerator.h"

using ::lcs_solver::algorithms::llcs::LLCS2_MC;
using ::lcs_solver::algorithms::llcs::LLCS2_MC_O1_SYNC;
using ::lcs_solver::algorithms::llcs::LLCS2_SR_MQ;
using ::lcs_solver::algorithms::llcs::LLCS2_SR_RMQ;
using ::lcs_solver::algorithms::llcs::LLCS2_SA_MQ;
using ::lcs_solver::algorithms::llcs::LLCS2_SA_RMQ;

using ::lcs_solver::constraints::local::Constraint_MC_1C;
using ::lcs_solver::constraints::local::Constraint_MC_O1C_SYNC;

using ::lcs_solver::util::ParamGenerator;
using ::lcs_solver::algorithms::BaseSolution;
using ConstraintType = ::lcs_solver::constraints::ConstraintType;

/*******************************************************************************
 * Benchmarks LLCS2_MC_k
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_MC_k(::benchmark::State &state) {
  const auto m = 200;
  const auto n = 200;
  const auto k = static_cast<size_t>(state.range(0));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    auto param = ParamGenerator::genWithUniqSol(
        ConstraintType::MC,
        seed, // Seed used to init random number generators
        1, // Number of AlgoParams to generate
        {{m, m}, {n, n}}, // Bounds on string lengths
        {k, k}, // Length of the LCS
        {'x', 'a', 'b'} // Alphabet for strings
    );

    const auto &spv = param.back().pointers;
    auto &map = param.back().map;
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
 * Benchmarks BM_LLCS2_MC_O1_SYNC with MC_1C Constraint
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_MC_O1_SYNC_1c(::benchmark::State &state) {
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
    auto gaps = lcs_solver::constraints::Get<Constraint_MC_1C>(map)->GetGapVector();
    map[ConstraintType::MC_O1C_SYNC] = std::make_shared<Constraint_MC_O1C_SYNC>(gaps);
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

std::u32string MakeUniqueStr(const std::size_t sig, const std::size_t len) {
  std::u32string result;
  result.reserve(sig);
  char32_t c = U'A'; // Startzeichen
  for (std::size_t i = 0; i < sig; ++i) {
    result.push_back(c++);
  }
  while (result.length() < len) {
    const auto symb = result.back();
    result.push_back(symb);
  }

  return result;
}

/*******************************************************************************
 * Benchmarks LLCS2_SR_MQ
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_SR_MQ_sig(::benchmark::State &state) {
  const auto m = 200;
  const auto n = 200;
  const auto sig =  static_cast<size_t>(state.range(0));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    lcs_solver::util::String s1 = MakeUniqueStr(sig, m);
    lcs_solver::util::String s2 = MakeUniqueStr(sig, m);
    std::vector<lcs_solver::util::String> vec = {s1,s2};

    lcs_solver::util::StrPtrVector spv;
    spv.emplace_back(std::make_shared<lcs_solver::util::String>(s1));
    spv.emplace_back(std::make_shared<lcs_solver::util::String>(s2));
    auto map = ParamGenerator::genRelaxedToStdMap(ConstraintType::SIGMA_R, vec);

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
 * Benchmarks BM_LLCS2_SA_MQ_sig
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_SA_MQ_sig(::benchmark::State &state) {
  const auto m = 200;
  const auto n = 200;
  const auto sig =  static_cast<size_t>(state.range(0));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    lcs_solver::util::String s1 = MakeUniqueStr(sig, m);
    lcs_solver::util::String s2 = MakeUniqueStr(sig, m);
    std::vector<lcs_solver::util::String> vec = {s1,s2};

    lcs_solver::util::StrPtrVector spv;
    spv.emplace_back(std::make_shared<lcs_solver::util::String>(s1));
    spv.emplace_back(std::make_shared<lcs_solver::util::String>(s2));
    auto map = ParamGenerator::genRelaxedToStdMap(ConstraintType::SIGMA, vec);

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
 * Benchmarks BM_LLCS2_SA_RMQ_sig
 * @param state Reference to the benchmark state object provided by Benchmark
 ******************************************************************************/
void BM_LLCS2_SA_RMQ_sig(::benchmark::State &state) {
  const auto m = 200;
  const auto n = 200;
  const auto sig =  static_cast<size_t>(state.range(0));
  size_t seed = std::random_device{}();
  for (auto _ : state) {
    //=== Setup of Benchmark====================================================
    state.PauseTiming();
    seed++;
    lcs_solver::util::String s1 = MakeUniqueStr(sig, m);
    lcs_solver::util::String s2 = MakeUniqueStr(sig, m);
    std::vector<lcs_solver::util::String> vec = {s1,s2};

    lcs_solver::util::StrPtrVector spv;
    spv.emplace_back(std::make_shared<lcs_solver::util::String>(s1));
    spv.emplace_back(std::make_shared<lcs_solver::util::String>(s2));
    auto map = ParamGenerator::genRelaxedToStdMap(ConstraintType::SIGMA, vec);

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

static void CustomArguments(benchmark::internal::Benchmark *b) {
  for (int i = 1; i <= 110; i += 10)
      b->Args({i});
}

static void CustomArguments2(benchmark::internal::Benchmark *b) {
  for (int i = 8; i <= 256; i *= 2)
    for (int j = 8; j <= 256; j *= 2)
      b->Args({i, j});
}

static void CustomArguments3(benchmark::internal::Benchmark *b) {
  for (int i = 2; i <= 20; i += 1)
    b->Args({i});
}

BENCHMARK(BM_LLCS2_MC_k)->Apply(CustomArguments);
 BENCHMARK(BM_LLCS2_MC_O1_SYNC_1c)->Apply(CustomArguments2);
BENCHMARK(BM_LLCS2_SR_MQ_sig)->Apply(CustomArguments3);
BENCHMARK(BM_LLCS2_SA_MQ_sig)->Apply(CustomArguments3);
BENCHMARK(BM_LLCS2_SA_RMQ_sig)->Apply(CustomArguments3);

BENCHMARK_MAIN();