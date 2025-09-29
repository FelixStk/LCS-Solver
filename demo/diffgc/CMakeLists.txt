# Demo: Diff with Gap Constraints

set(EXECUTABLE_DIFF diffgc)

# Create a library for reusable Diff logic
add_library(diffgc_lib STATIC
        diffgc.cc
        FileAdapter.cc
)
target_compile_options(diffgc_lib PRIVATE -Wall)
target_include_directories(diffgc_lib PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(diffgc_lib PUBLIC ${PROJECT_NAME})

# Define the diff executable
add_executable(${EXECUTABLE_DIFF} main.cc)
target_compile_options(${EXECUTABLE_DIFF} PRIVATE -Wall) # Maybe add -Werror
target_link_libraries(${EXECUTABLE_DIFF} PRIVATE diffgc_lib)

configure_file(res/file1.cc ${CMAKE_CURRENT_BINARY_DIR}/res/file1.cc COPYONLY)
configure_file(res/file2.cc ${CMAKE_CURRENT_BINARY_DIR}/res/file2.cc COPYONLY)
configure_file(res/constraint.json ${CMAKE_CURRENT_BINARY_DIR}/res/constraint.json COPYONLY)
configure_file(res/config_mc.json ${CMAKE_CURRENT_BINARY_DIR}/res/config_mc.json COPYONLY)

if (NOT BUILD_LIB_STATIC)
    copy_runtime_dependencies(${EXECUTABLE_DIFF})
endif ()