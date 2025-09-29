cmake_minimum_required(VERSION 3.28.3)
project(lcs_solver VERSION 1.0.0 LANGUAGES C CXX)

# Add cmake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(FETCH_JSON "Use FetchContent to provide json" ON)
option(FETCH_GTEST "Use FetchContent to provide json" ON)
option(FETCH_BENCH "Use FetchContent to provide benchmark" ON)

option(BUILD_STATIC_LIB "Whether to build lcs_solver as STATIC or SHARED" ON)
option(BUILD_TESTS "Build tests of the libary" ON)
option(BUILD_BENCH "Build benchmarks of the libary" ON)
option(BUILD_DEMOS "Build demo executables" ON)
option(BUILD_DOCS "Build Doxygen documentation" ON)

include(GNUInstallDirs)

if(FETCH_JSON OR FETCH_GTEST OR FETCH_BENCH)
    include(FetchContent)
endif ()

find_package(nlohmann_json CONFIG QUIET)
if(nlohmann_json_FOUND)
    message(STATUS "Found installed nlohmann_json")
elseif (NOT FETCH_JSON)
    message(FATAL_ERROR "nlohmann_json not found")
else()
    message(STATUS "Fetching nlohmann_json library from GitHub (v3.11.3)")
    include(FetchContent)
    FetchContent_Declare(
            nlohmann_json
            GIT_REPOSITORY https://github.com/nlohmann/json.git
            GIT_TAG v3.11.3
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif ()

add_subdirectory(lib)

if (BUILD_DEMOS)
    add_subdirectory(demo)      # diff, dna and other demos
endif ()

if (BUILD_TESTS)
    enable_testing()
    include(CTest)              # required for gtest_discover_tests
    add_subdirectory(tests)     # tests the library
endif ()

if (BUILD_BENCH)
    add_subdirectory(benchmark) # benchmarks the libary
endif ()

if(BUILD_DOCS)
    add_subdirectory(docs)      # doxygen-docs target
endif ()