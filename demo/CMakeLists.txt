# Provide copy_runtime_dependencies(target)
if(NOT BUILD_LIB_STATIC)
    include(CopyRuntimeDeps)
endif ()

# Add example demos
add_subdirectory(dev)
add_subdirectory(diffgc)
add_subdirectory(dna)
add_subdirectory(solver)