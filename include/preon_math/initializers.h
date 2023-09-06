/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "basics.h"
#include "scalar_simd.h"

#ifdef PREONMATH_UNITS_INTEGRATION
    #include "core/utility/units_math_scalar.h"  // TODO: Including a file in "core" is not so nice and needs to be changed when releasing a new version of PreonMath.
#endif  // PREONMATH_UNITS_INTEGRATION

namespace Preon
{
namespace Math
{
    template<typename T>
    PREONMATH_FORCEINLINE T scalar_zero()
    {
#ifdef PREONMATH_ENABLE_SIMD
        using Type = decltype(PrimitiveType(std::declval<T>()));
        if constexpr (is_simd_scalar<Type>::value)
            return T(Simd::zero<typename Simd::Scalar<Type>::type>());
        else
#endif  // PREONMATH_ENABLE_SIMD
            return T(0);
    }
}  // namespace Math
}  // namespace Preon
