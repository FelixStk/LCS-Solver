add_subdirectory(local)
add_subdirectory(global)
add_subdirectory(input)

target_sources(${PROJECT_NAME} PRIVATE
        ConstraintType.cc
        BaseConstraint.cc
        ConstraintFactory.cc
        ConstraintMap.cc
)