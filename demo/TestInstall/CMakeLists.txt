# TestProject/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(TestInstall VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find the installed package
find_package(lcs_solver REQUIRED)

# Create a test executable
add_executable(test_install main.cc)

# Link against the found library
target_link_libraries(test_install PRIVATE lcs_solver::lcs_solver)
