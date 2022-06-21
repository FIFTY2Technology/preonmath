# Following http://mariobadr.com/creating-a-header-only-library-with-cmake.html
cmake_minimum_required(VERSION 3.5)

project(preon_math)

add_library(preon_math INTERFACE)
set(header_files
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/basics.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/cast.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/compile_helper.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/euler_fwd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/euler.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/initializers.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/math_misc.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/math_utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix_fwd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix_simd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix_utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix11.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix22.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix33.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/matrix44.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/quat_fwd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/quat.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/rotation_utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/scalar_avx.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/scalar_simd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/scalar_sse.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/simd_scalar_wrapper_std.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/simd_scalar_wrapper.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/traits.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec_fwd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec_simd.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec_utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec2.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec3.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/preon_math/vec4.h
    )
target_sources(preon_math INTERFACE "$<BUILD_INTERFACE:${header_files}>")

target_include_directories(preon_math INTERFACE include/)
target_compile_features(preon_math INTERFACE cxx_std_17)

# Adapted from https://stackoverflow.com/questions/18233513/cmake-g-set-cxx-flags-02-but-03-is-still-there
macro(remove_cxx_release_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endmacro()
macro(add_cxx_release_flag flag)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${flag}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${flag}")
endmacro()
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /WX")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -msse4")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall")
    add_cxx_release_flag("-O3")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -msse4")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wno-ignored-attributes")
    remove_cxx_release_flag("-O2")  # Force O3 flag, some CMake choose O2 by default.
    add_cxx_release_flag("-O3")
else()
    message(FATAL_ERROR "Compiler unknown and not supported.")
endif()

# Tests.
add_subdirectory(test)
