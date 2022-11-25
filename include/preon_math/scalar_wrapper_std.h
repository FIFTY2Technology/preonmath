/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "scalar_simd.h"

#include <cmath>

namespace Preon
{
namespace Math
{
    // sqrt
#ifdef PREONMATH_CUDA
    template<typename T>
    PREONMATH_DEVICE T sqrt(T val)
    {
        return ::sqrt(val);
    }
    template<typename T>
    PREONMATH_DEVICE T abs(T val)
    {
        return ::abs(val);
    }
    template<typename T>
    PREONMATH_DEVICE T min(T val)
    {
        return ::min(val);
    }
    template<typename T>
    PREONMATH_DEVICE T max(T val)
    {
        return ::max(val);
    }
#else
#ifdef PREONMATH_ENABLE_SIMD
    using Simd::sqrt;
#endif // PREONMATH_ENABLE_SIMD
    using std::max;
    using std::min;
    using std::sqrt;
#endif

    // add other std functions that are implemented in Simd...
}  // namespace Math
}  // namespace Preon
