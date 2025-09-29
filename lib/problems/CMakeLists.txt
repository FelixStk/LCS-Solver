add_subdirectory(AdamsonEtAl)

target_sources(${PROJECT_NAME} PRIVATE
        BaseProblem.cpp
        LCS_Classic.cpp
        ProblemFactory.cc
)