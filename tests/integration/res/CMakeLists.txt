configure_file(gap.json ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
configure_file(alph.json ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
configure_file(test.fasta ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)

configure_file(file1.cc ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
configure_file(file2.cc ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
configure_file(constraint.json ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
configure_file(config_mc.json ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)