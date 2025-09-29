function(copy_runtime_dependencies target)
    get_target_property(target_libs ${target} LINK_LIBRARIES)
    foreach (lib ${target_libs})
        if (TARGET ${lib})
            get_target_property(lib_type ${lib} TYPE)
            if (lib_type STREQUAL "SHARED_LIBRARY")
                add_custom_command(TARGET ${target} POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        "$<TARGET_FILE:${lib}>"
                        "$<TARGET_FILE_DIR:${target}>"
                )
            endif ()
        endif ()
    endforeach ()
endfunction()