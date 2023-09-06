/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"
#include "vec.h"
#include "matrix.h"
#include "simd_scalar_wrapper_definitions.h"

namespace Preon
{
namespace Math
{
    namespace SimdWrapper
    {
        template<>
        struct Scalar<float>
        {
            using type = float;
        };
        template<>
        struct Scalar<double>
        {
            using type = double;
        };

        // scalar make wrapper
        template<>
        PREONMATH_FORCEINLINE float make<float>(float val)
        {
            return val;
        }

        template<>
        PREONMATH_FORCEINLINE double make<double>(float val)
        {
            return val;
        }
        template<>
        PREONMATH_FORCEINLINE double make<double>(double val)
        {
            return val;
        }

        // also overload the following for convenience
        template<>
        PREONMATH_FORCEINLINE float make<float>(double val)
        {
            return val;
        }

        // vector make wrapper
        template<>
        PREONMATH_FORCEINLINE vec3f make<vec3f>(const vec3f& val)
        {
            return val;
        }

        template<>
        PREONMATH_FORCEINLINE float zero<float>()
        {
            return 0.0f;
        }
        template<>
        PREONMATH_FORCEINLINE double zero<double>()
        {
            return 0.0;
        }

        // scalar hSum wrapper
        PREONMATH_FORCEINLINE float hSum(const float& val)
        {
            return val;
        }

        // vector hSum wrapper
        PREONMATH_FORCEINLINE vec3f hSum(const vec3f& val)
        {
            return val;
        }

        // matrix hSum wrapper
        PREONMATH_FORCEINLINE matrix33f hSum(const matrix33f& val)
        {
            return val;
        }

        // choose
        PREONMATH_FORCEINLINE float choose(const bool mask, const float v1, const float v2)
        {
            return mask ? v1 : v2;
        }

        template<PrMathSize D>
        PREONMATH_FORCEINLINE vec<D, float> choose(const bool mask, const vec<D, float>& v1, const vec<D, float>& v2)
        {
            return mask ? v1 : v2;
        }

        // AND mask
        PREONMATH_FORCEINLINE float and_mask(bool mask, float val)
        {
            return choose(mask, val, 0);
        }
        PREONMATH_FORCEINLINE float and_mask(float val, bool mask)
        {
            return choose(mask, val, 0);
        }

        // max
        PREONMATH_FORCEINLINE float max(float a, float b)
        {
            return std::max(a, b);
        }

        // min
        PREONMATH_FORCEINLINE float min(float a, float b)
        {
            return std::min(a, b);
        }

        // cosine
        PREONMATH_FORCEINLINE float cos(float val)
        {
            return std::cos(val);
        }

        // exp
        PREONMATH_FORCEINLINE float exp(float val)
        {
            return std::exp(val);
        }
    }  // namespace SimdWrapper
}  // namespace Math
}  // namespace Preon

#undef SIMDWRAPPER_FORCEINLINE