##=== Find/Download GTest ======================================================
find_package(GTest QUIET)

set(NEED_FETCH_GTEST FALSE)
if (GTest_FOUND)
    if (NOT TARGET GTest::gmock OR NOT TARGET GTest::gmock_main)
        message(STATUS "GTest found but GMock is missing.")
        set(NEED_FETCH_GTEST TRUE)
    else()
        message(STATUS "Found installed GTest with GMock")
    endif()
else ()
    message(STATUS "GTest not found")
    set(NEED_FETCH_GTEST TRUE)
endif ()

if(NEED_FETCH_GTEST)
    if (NOT FETCH_GTEST)
        message(FATAL_ERROR "GTest or GMock not found and FETCH_GTEST is OFF")
    endif()

    message(STATUS "Fetching Google Test library from GitHub (v1.16.0)")
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.16.0
    )
    # Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif ()
set(gtest_LIB GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main)

##=== Setup Tests ==============================================================
include(GoogleTest)

add_executable(LCS_Test)
target_include_directories(LCS_Test PRIVATE
        ${CMAKE_SOURCE_DIR}/include/lcs_solver
        ${CMAKE_SOURCE_DIR}
)

add_subdirectory(integration)
add_subdirectory(LLCS)
add_subdirectory(LCS)
add_subdirectory(structures)
add_subdirectory(util)

target_link_libraries(LCS_Test
        ${gtest_LIB}
        ${PROJECT_NAME}
        dna_lib
        diffgc_lib
        solver_lib
)

gtest_discover_tests(LCS_Test
        PROPERTIES TIMEOUT 120
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        EXTRA_ARGS --gtest_color=yes
        TEST_LIST all_tests
        DISCOVERY_MODE PRE_TEST
        #PROPERTIES ENVIRONMENT "GTEST_COLOR=yes"
)