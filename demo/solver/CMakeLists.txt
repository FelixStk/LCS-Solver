# Demo: Run Algorithms via the Terminal and save Problems

set(EXECUTABLE_SOLVER solver)

# Create a library for reusable logic and tests
add_library(solver_lib STATIC
        solver.cc
)
target_compile_options(solver_lib PRIVATE -Wall)
target_include_directories(solver_lib PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(solver_lib PUBLIC ${PROJECT_NAME})

# Define the solver executable target
add_executable(${EXECUTABLE_SOLVER} main.cc)
target_compile_options(${EXECUTABLE_SOLVER} PRIVATE -Wall)
target_link_libraries(${EXECUTABLE_SOLVER} PRIVATE solver_lib)

if(NOT BUILD_LIB_STATIC)
    copy_runtime_dependencies(${EXECUTABLE_SOLVER})
endif ()