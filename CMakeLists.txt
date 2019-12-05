# Following http://mariobadr.com/creating-a-header-only-library-with-cmake.html
cmake_minimum_required(VERSION 3.1)

project(preon_math)

add_library(preon_math INTERFACE)
set(header_files
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/compile_helper.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/basics.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/euler.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/euler_fwd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/math_misc.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/math_utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix11.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix22.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix33.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix44.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix_fwd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix_simd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix_util.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/quat.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/quat_fwd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/rotation_utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/scalar_avx.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/scalar_simd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/scalar_sse.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/simd_scalar_wrapper.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/simd_scalar_wrapper_std.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/traits.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec2.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec3.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec4.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec_fwd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec_simd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec_utils.h
    )
target_sources(preon_math INTERFACE "$<BUILD_INTERFACE:${header_files}>")

target_include_directories(preon_math INTERFACE include/)

# Tests.
add_subdirectory(test)
