/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

// #define PREONMATH_USED_IN_PREON_CODE
// #define PREONMATH_QT_INTEGRATION

// This is something from the preon code. We will try to get rid of it in the preon math code in the future.
// Includes or defines macros etc depending on if the library is used
// in the Preon code base.
#ifdef PREONMATH_USED_IN_PREON_CODE
    #include "core/utility/error_handling.h"

// Needed for the "using" below.
namespace Preon
{
namespace Math
{
}  // namespace Math
}  // namespace Preon
using namespace Preon::Math;
#endif  // PREONMATH_USED_IN_PREON_CODE
#ifndef THROW_EXCEPTION
    #define THROW_EXCEPTION(condition, exception) \
        {                                         \
            if (condition)                        \
            {                                     \
                throw exception;                  \
            }                                     \
        }
#endif  // THROW_EXCEPTION

// Define macros for identifying msvc compiler.
// See https://blog.kowalczyk.info/article/j/guide-to-predefined-macros-in-c-compilers-gcc-clang-msvc-etc..html
#if defined(_MSC_VER) && !defined(__clang__)
    #define PREONMATH_COMPILER_MSVC
#endif

#ifndef PREONMATH_FORCEINLINE
    #ifdef PREONMATH_COMPILER_MSVC
        #define PREONMATH_FORCEINLINE __forceinline
    #else
        #define PREONMATH_FORCEINLINE __attribute__((always_inline)) inline
    #endif  // PREONMATH_COMPILER_MSVC
#endif  // PREONMATH_FORCEINLINE

#ifndef PREONMATH_NOINLINE
    #ifdef PREONMATH_COMPILER_MSVC
        #define PREONMATH_NOINLINE __declspec(noinline)
    #else
        #define PREONMATH_NOINLINE __attribute__((noinline))
    #endif  // PREONMATH_COMPILER_MSVC
#endif  // PREONMATH_NOINLINE

#ifdef PREONMATH_QT_INTEGRATION
    #include <QtGlobal>
#endif  // PREONMATH_QT_INTEGRATION
