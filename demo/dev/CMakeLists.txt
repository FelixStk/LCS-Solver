# devStuff runs running arbitrary code for testing and development

add_executable(devStuff main.cc)
#target_compile_options(devStuff PRIVATE -Wall -Werror)
target_link_libraries(devStuff PRIVATE
        lcs_solver::lcs_solver
        nlohmann_json::nlohmann_json
)

if(NOT BUILD_LIB_STATIC)
    copy_runtime_dependencies(devStuff)
endif ()