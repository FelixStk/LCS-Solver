# Build lcs_solver libary
if(BUILD_STATIC_LIB)
    add_library(${PROJECT_NAME} STATIC)
else ()
    add_library(${PROJECT_NAME} SHARED)
endif ()

add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror)

target_link_libraries(${PROJECT_NAME} nlohmann_json)

target_include_directories(${PROJECT_NAME}
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/lcs_solver>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/lcs_solver>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/lib
)

# Add libary implementation directories
add_subdirectory(algorithms)
add_subdirectory(constraints)
add_subdirectory(problems)
add_subdirectory(structures)
add_subdirectory(util)

# Add Header With Version Number
configure_file(${PROJECT_SOURCE_DIR}/cmake/LibVersion.h.in
        ${PROJECT_BINARY_DIR}/include/lcs_solver/LibVersion.h
)

# Install Dependency if it is not installed
if(NOT nlohmann_json_FOUND AND FETCH_JSON)
    install(TARGETS nlohmann_json
            EXPORT "${PROJECT_NAME}Targets"
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

# Install lib
install(TARGETS ${PROJECT_NAME}
        EXPORT "${PROJECT_NAME}Targets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # for STATIC libaries
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # for SHARED libaries
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # for executables
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install version headers
install(FILES ${PROJECT_BINARY_DIR}/include/lcs_solver/LibVersion.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install exported targets (for find_package)
install(EXPORT "${PROJECT_NAME}Targets"
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME} # or CMAKE_INSTALL_DATAROOTDIR?
)

# Create and install config file
include(CMakePackageConfigHelpers)
configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME} # or CMAKE_INSTALL_DATAROOTDIR?
)

# Install Licenses
install(DIRECTORY "${CMAKE_SOURCE_DIR}/LICENSES"
        DESTINATION "${CMAKE_INSTALL_DATADIR}"
        FILES_MATCHING PATTERN "*")
install(FILES "${CMAKE_SOURCE_DIR}/LICENSE"
        DESTINATION "${CMAKE_INSTALL_DATADIR}")

# Combine Licenses
file(WRITE "${CMAKE_BINARY_DIR}/CombinedLicenses.txt" "")
file(READ "${CMAKE_SOURCE_DIR}/LICENSE" MAIN_LICENSE)
file(APPEND "${CMAKE_BINARY_DIR}/CombinedLicenses.txt" "${MAIN_LICENSE}\n\n")

file(GLOB LICENSE_FILES "${CMAKE_SOURCE_DIR}/LICENSES/*")
foreach (LIC ${LICENSE_FILES})
    file(READ "${LIC}" CONTENTS)
    file(RELATIVE_PATH REL_PATH "${CMAKE_SOURCE_DIR}" "${LIC}")
    file(APPEND "${CMAKE_BINARY_DIR}/CombinedLicenses.txt"
            "----- ${REL_PATH} -----\n\n${CONTENTS}\n\n")
endforeach ()

# CPack Stuff
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_BINARY_DIR}/CombinedLicenses.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
include(CPack)