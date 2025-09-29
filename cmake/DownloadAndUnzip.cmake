include(FetchContent)

function(download_and_unzip name url)
    set(download_dir "${CMAKE_BINARY_DIR}/downloads")
    set(zip_file "${download_dir}/${name}")
    set(unzip_dir "${CMAKE_CURRENT_BINARY_DIR}/")

    file(MAKE_DIRECTORY "${download_dir}")
    file(MAKE_DIRECTORY "${unzip_dir}")

    # Download the file if not already present
    if(NOT EXISTS "${zip_file}")
        message(STATUS "Downloading ${name} from ${url}")
        file(DOWNLOAD "${url}" "${zip_file}" SHOW_PROGRESS STATUS status LOG log)
        list(GET status 0 status_code)
        if(NOT status_code EQUAL 0)
            message(FATAL_ERROR "Download failed: ${log}")
        endif()
    endif()

    # Unzip the file (add --strip-components=1 to extract in the current dir)
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf "${zip_file}"
            WORKING_DIRECTORY "${unzip_dir}"
            RESULT_VARIABLE unzip_result
    )

    if(NOT unzip_result EQUAL 0)
        message(FATAL_ERROR "Unzipping ${zip_file} failed with exit code ${unzip_result}")
    else()
        message(STATUS "Unzipped ${zip_file} to ${unzip_dir}")
    endif()
endfunction()
