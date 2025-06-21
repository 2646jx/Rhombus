macro(add_module_test _name)
    add_executable(${_name} "test_${_name}.cpp")
    target_include_directories(${_name} PUBLIC ${CMAKE_SOURCE_DIR}/src)
    target_link_libraries(${_name} PUBLIC Rhombus)
endmacro()

add_module_test(matmul)
add_module_test(matvec)
