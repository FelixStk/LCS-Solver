# Enable doxygen doc target

set(DOXYFILE "${PROJECT_SOURCE_DIR}/docs/Doxyfile")

find_package(Doxygen REQUIRED OPTIONAL_COMPONENTS dot mscgen dia)

doxygen_add_docs(
        doxygen-docs ALL
        ${PROJECT_SOURCE_DIR}/lib
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/README.md
        COMMENT "Generating documentation - see ${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
        CONFIG_FILE "${DOXYFILE}"
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
)