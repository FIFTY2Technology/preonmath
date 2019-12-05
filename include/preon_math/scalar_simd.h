/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_SCALAR_SIMD_H
#define PREONMATH_SCALAR_SIMD_H

#include "compile_helper.h"

// Unfortunately, AVX seems to offer no benefit for us at the moment.
// #define USE_AVX

#ifdef USE_AVX
#include "scalar_avx.h"
#else
#include "scalar_sse.h"
#endif

#include <algorithm>

namespace Preon
{
    namespace Math
    {
        template<typename T>
        struct is_simd_scalar
        {
            static const bool value = false;
        };
        template<>
        struct is_simd_scalar <float_simd>
        {
            static const bool value = true;
        };
        template<>
        struct is_simd_scalar <double_simd>
        {
            static const bool value = true;
        };

        template <typename T, typename E>
        using enable_if_simd = typename std::enable_if<is_simd_scalar<T>::value, E>::type;

        template <typename T, typename E>
        using enable_if_non_simd = typename std::enable_if<!is_simd_scalar<T>::value, E>::type;

        namespace Simd
        {
            template<typename T>
            struct Scalar {};
            template<>
            struct Scalar<float_simd>
            {
                using type = float;
            };
            template<>
            struct Scalar<double_simd>
            {
                using type = double;
            };

            //! scalar_cast converts a non-SIMD scalar (float, double) to type T.
            //! If T is a scalar type, nothing will happen, but if T is a SIMD type, the scalar will be written to all slots of the SIMD register.
            template<typename A, typename B>
            PREONMATH_FORCEINLINE static B scalar_cast(enable_if_non_simd<B, A> in) { return static_cast<B>(in); }
            template<typename A, typename B>
            PREONMATH_FORCEINLINE static B scalar_cast(enable_if_simd<B, A> in) { return Simd::make<A>(in); }

            PREONMATH_FORCEINLINE float_simd choose(const float_simd mask, const float_simd v1, const float_simd v2)
            {
                return blendv(v2, v1, mask);
            }

            PREONMATH_FORCEINLINE double_simd choose(const double_simd mask, const double_simd v1, const double_simd v2)
            {
                return blendv(v2, v1, mask);
            }

            PREONMATH_FORCEINLINE float hMax(float_simd v)
            {
                float values[Simd::Register<float>::size];
                store(v, values);
                float out = values[0];
                for (int i = 1; i < Simd::Register<float>::size; i++)
                    out = std::max(out, values[i]);
                return out;
            }

            PREONMATH_FORCEINLINE float hMin(float_simd v)
            {
                float values[Simd::Register<float>::size];
                store(v, values);
                float out = values[0];
                for (int i = 1; i < Simd::Register<float>::size; i++)
                    out = std::min(out, values[i]);
                return out;
            }

            //! Retrieves the next size (<=) which aligns with the register size.
            template<typename T>
            PREONMATH_FORCEINLINE unsigned int alignedSize(unsigned int size)
            {
                return (size >> Register<T>::expo) << Register<T>::expo;
            }
        }
    }  // namespace Math
}  // namespace Preon

#endif  // PREONMATH_SCALAR_SIMD_H
