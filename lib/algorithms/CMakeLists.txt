add_subdirectory(LCS)
add_subdirectory(LLCS)
add_subdirectory(solutions)

target_sources(${PROJECT_NAME} PRIVATE
        AlgoFactory.cc
        BaseSolution.cc
        BaseAlgorithm.cc
)