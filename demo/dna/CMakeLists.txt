# Demo: Reconstruct of evolutionary trees based on protein sequences

set(EXECUTABLE_DNA dna)

# Download benchmark data
include(DownloadAndUnzip)
download_and_unzip("swisstree.zip" "https://afproject.org/media/genetree/swisstree/dataset/swisstree.zip")

# Create a library for reusable DNA logic
add_library(dna_lib STATIC
        dna.cc
        OutputParser.cc
)
target_compile_options(dna_lib PRIVATE -Wall)
target_include_directories(dna_lib PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(dna_lib PUBLIC ${PROJECT_NAME})

# Define the main target
add_executable(${EXECUTABLE_DNA} main.cc)
target_compile_options(${EXECUTABLE_DNA} PRIVATE -Wall)
target_link_libraries(${EXECUTABLE_DNA} PRIVATE dna_lib)

if (NOT BUILD_LIB_STATIC)
    copy_runtime_dependencies(${EXECUTABLE_DNA})
endif ()