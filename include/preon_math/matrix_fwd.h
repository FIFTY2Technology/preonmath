/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"
#include "scalar_simd.h"

#include <cstddef>  // For PrMathSize

namespace Preon
{
namespace Math
{
    template<PrMathSize M, PrMathSize N, typename T>
    class matrix;

    // Matrix 1x1
    template<typename T>
    using matrix11 = matrix<1, 1, T>;
    using matrix11f = matrix11<float>;
    using matrix11d = matrix11<double>;
#ifdef PREONMATH_ENABLE_SIMD
    using matrix11_Simd = matrix<1, 1, float_simd>;
#endif

    // Matrix 2x2
    template<typename T>
    using matrix22 = matrix<2, 2, T>;
    using matrix22f = matrix22<float>;
    using matrix22d = matrix22<double>;
#ifdef PREONMATH_ENABLE_SIMD
    using matrix22_Simd = matrix<2, 2, float_simd>;
#endif

    // Matrix 3x3
    template<typename T>
    using matrix33 = matrix<3, 3, T>;
    using matrix33f = matrix33<float>;
    using matrix33d = matrix33<double>;
#ifdef PREONMATH_ENABLE_SIMD
    using matrix33_Simd = matrix33<float_simd>;
#endif

    // Matrix 4x4
    template<typename T>
    using matrix44 = matrix<4, 4, T>;
    using matrix44f = matrix44<float>;
    using matrix44d = matrix44<double>;
}  // namespace Math
}  // namespace Preon
