## === Configure Benchmark =====================================================
find_package(Benchmark QUIET)
if (Benchmark_FOUND)
    message(STATUS "Found installed Google Benchmark library")
elseif (NOT FETCH_GTEST)
    message(FATAL_ERROR "Google Benchmark not found")
else ()
    message(STATUS "Fetching Google Benchmark library from GitHub (v1.9.2)")
    set(BENCHMARK_ENABLE_TESTING off) # suppress benchmark internal tests
    FetchContent_Declare(
            google_benchmark
            GIT_REPOSITORY https://github.com/google/benchmark.git
            GIT_TAG v1.9.2
    )
    FetchContent_MakeAvailable(google_benchmark)
endif ()
set(benchmark_LIB benchmark::benchmark)

##=== Benchmark Executables ====================================================

# Helper function to add a benchmark executable
function(add_benchmark target source_file)
    add_executable(${target} ${CMAKE_CURRENT_SOURCE_DIR}/${source_file})
    target_include_directories(${target} PRIVATE ${CMAKE_SOURCE_DIR}/include/lcs_solver)
    target_link_libraries(${target} PRIVATE ${benchmark_LIB} ${PROJECT_NAME})
endfunction()

# Add benchmark executables
add_benchmark(Bench_RMQ bm_rmq.cc)
add_benchmark(Bench_LLCS bm_llcs.cc)
add_benchmark(Bench_SegTree bm_segtree.cc)
add_benchmark(Bench_Fine bm_llcs_mc.cc)
