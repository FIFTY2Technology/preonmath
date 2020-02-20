set(test_src_files
        ${CMAKE_CURRENT_SOURCE_DIR}/test_math_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/test_vec2.cpp
    )

add_executable(preon_math_tests main.cpp ${test_src_files})
target_link_libraries(preon_math_tests PUBLIC preon_math)

add_test(NAME test COMMAND type_safe_test)
