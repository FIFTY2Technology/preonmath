/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_SCALAR_SSE_H
#define PREONMATH_SCALAR_SSE_H

#include "compile_helper.h"

#include "basics.h"

#ifdef PREONMATH_COMPILER_MSVC
#include <intrin.h>
#else  // GCC and Clang
#include <pmmintrin.h>
#include <smmintrin.h>
#endif  // COMPILER

#include <type_traits>

namespace Preon
{
    namespace Math
    {
        using float_simd = __m128;
        using double_simd = __m128d;

        namespace Simd
        {
            template<typename scalar>
            using simd_if_float = typename std::enable_if<!std::is_same<scalar, double>::value, float_simd>::type;
            template<typename scalar>
            using simd_if_double = typename std::enable_if<std::is_same<scalar, double>::value, double_simd>::type;

            template<typename T>
            struct Register {};
            template<>
            struct Register<float>
            {
                enum { size = 4 };
                enum { expo = 2 };
                using type = float_simd;
            };
            template<>
            struct Register<double>
            {
                enum { size = 2 };
                enum { expo = 1 };
                using type = double_simd;
            };

            // We demand explicit float / double typing here, because otherwise it would be too easy to create the wrong type by accident.
            template<typename scalar>
            PREONMATH_FORCEINLINE simd_if_float<scalar> make(float a) { return _mm_set1_ps(a); }
            template<typename scalar>
            PREONMATH_FORCEINLINE simd_if_double<scalar> make(double a) { return _mm_set1_pd(a); }

            PREONMATH_FORCEINLINE float_simd load(const float* values) { return _mm_loadu_ps(values); }
            PREONMATH_FORCEINLINE double_simd load(const double* values) { return _mm_loadu_pd(values); }

            PREONMATH_FORCEINLINE void store(float_simd input, float* values) { _mm_storeu_ps(values, input); }
            PREONMATH_FORCEINLINE void store(double_simd input, double* values) { _mm_storeu_pd(values, input); }

            template <typename scalar>
            PREONMATH_FORCEINLINE simd_if_float<scalar> zero() { return _mm_setzero_ps(); }

            template <typename scalar>
            PREONMATH_FORCEINLINE simd_if_double<scalar> zero() { return _mm_setzero_pd(); }

            PREONMATH_FORCEINLINE float_simd and_mask(float_simd a, float_simd b) { return _mm_and_ps(a, b); }
            PREONMATH_FORCEINLINE float_simd or_mask(float_simd a, float_simd b) { return _mm_or_ps(a, b); }
            PREONMATH_FORCEINLINE float_simd xor_mask(float_simd a, float_simd b) { return _mm_xor_ps(a, b); }
            PREONMATH_FORCEINLINE float_simd blendv(float_simd a, float_simd b, float_simd mask) { return _mm_blendv_ps(a, b, mask); }

            PREONMATH_FORCEINLINE double_simd and_mask(double_simd a, double_simd b) { return _mm_and_pd(a, b); }
            PREONMATH_FORCEINLINE double_simd or_mask(double_simd a, double_simd b) { return _mm_or_pd(a, b); }
            PREONMATH_FORCEINLINE double_simd xor_mask(double_simd a, double_simd b) { return _mm_xor_pd(a, b); }
            PREONMATH_FORCEINLINE double_simd blendv(double_simd a, double_simd b, double_simd mask) { return _mm_blendv_pd(a, b, mask); }

            template<class Getter>
            PREONMATH_FORCEINLINE void setScalars(float_simd* out, const Getter& func)
            {
                (*out) = _mm_setr_ps(func(0), func(1), func(2), func(3));
            }
            template<class Getter>
            PREONMATH_FORCEINLINE void setScalars(double_simd* out, const Getter& func)
            {
                (*out) = _mm_setr_pd(func(0), func(1));
            }

            template<class Getter>
            PREONMATH_FORCEINLINE void setScalars(float_simd* out, const Getter& func, unsigned int numElements, float fillValue)
            {
                (*out) = _mm_setr_ps(func(0),
                        (numElements > 1) ? func(1) : fillValue,
                        (numElements > 2) ? func(2) : fillValue,
                        (numElements > 3) ? func(3) : fillValue);
            }
            template<class Getter>
            PREONMATH_FORCEINLINE void setScalars(double_simd* out, const Getter& func,  unsigned int numElements, double fillValue)
            {
                (*out) = _mm_setr_pd(func(0),
                        (numElements > 1) ? func(1) : fillValue);
            }

            template<class Getter>
            PREONMATH_FORCEINLINE void storeScalars(float_simd input, const Getter& func)
            {
                float values[4];
                store(input, values);
                func(0) = values[0];
                func(1) = values[1];
                func(2) = values[2];
                func(3) = values[3];
            }
            template<class Getter>
            PREONMATH_FORCEINLINE void storeScalars(double_simd input, const Getter& func)
            {
                double values[2];
                store(input, values);
                func(0) = values[0];
                func(1) = values[1];
            }

            PREONMATH_FORCEINLINE float hSum(float_simd v)
            {
                float out;
                float_simd x = _mm_hadd_ps(v, v);
                x = _mm_hadd_ps(x, x);
                _mm_store_ss(&out, x);
                return out;
            }
            PREONMATH_FORCEINLINE double hSum(double_simd v)
            {
                double out;
                double_simd x = _mm_hadd_pd(v, v);
                _mm_store_sd(&out, x);
                return out;
            }

            PREONMATH_FORCEINLINE bool _testAllZeroFloat(const float_simd x)
            {
                return _mm_test_all_zeros(_mm_castps_si128(x), _mm_castps_si128(x)) != 0;
            }

            template<typename T>
            PREONMATH_FORCEINLINE bool testAllZero(const T x)
            {
                return _mm_test_all_zeros(_mm_castps_si128((float_simd)x), _mm_castps_si128((float_simd)x)) != 0;
            }

            PREONMATH_FORCEINLINE float_simd round(const float_simd x)
            {
                return _mm_round_ps(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
            }

            PREONMATH_FORCEINLINE double_simd round(const double_simd x)
            {
                return _mm_round_pd(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
            }

            PREONMATH_FORCEINLINE float_simd max(float_simd a, float_simd b)
            {
                return _mm_max_ps(a, b);
            }

            PREONMATH_FORCEINLINE double_simd max(double_simd a, double_simd b)
            {
                return _mm_max_pd(a, b);
            }

            PREONMATH_FORCEINLINE float_simd min(float_simd a, float_simd b)
            {
                return _mm_min_ps(a, b);
            }

            PREONMATH_FORCEINLINE double_simd min(double_simd a, double_simd b)
            {
                return _mm_min_pd(a, b);
            }

            PREONMATH_FORCEINLINE float_simd sqrt(float_simd a)
            {
                return _mm_sqrt_ps(a);
            }

            PREONMATH_FORCEINLINE double_simd sqrt(double_simd a)
            {
                return _mm_sqrt_pd(a);
            }

            PREONMATH_FORCEINLINE float_simd abs(float_simd a)
            {
                const float_simd signMask = _mm_set1_ps(-0.f);
                return _mm_andnot_ps(signMask, a);
            }

            PREONMATH_FORCEINLINE double_simd abs(double_simd a)
            {
                const double_simd signMask = _mm_set1_pd(-0.0);
                return _mm_andnot_pd(signMask, a);
            }
        }
        using preal_simd = Simd::Register<preal>::type;
    }  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_COMPILER_MSVC    // Operators are already defined on gcc / clang
PREONMATH_FORCEINLINE Preon::Math::float_simd operator+(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_add_ps(a, b); }
PREONMATH_FORCEINLINE Preon::Math::float_simd operator-(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_sub_ps(a, b); }
PREONMATH_FORCEINLINE Preon::Math::float_simd operator*(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_mul_ps(a, b); }
PREONMATH_FORCEINLINE Preon::Math::float_simd operator/(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_div_ps(a, b); }

PREONMATH_FORCEINLINE Preon::Math::float_simd& operator+=(Preon::Math::float_simd& a, Preon::Math::float_simd b) { return a = a + b; }
PREONMATH_FORCEINLINE Preon::Math::float_simd& operator-=(Preon::Math::float_simd& a, Preon::Math::float_simd b) { return a = a - b; }
PREONMATH_FORCEINLINE Preon::Math::float_simd& operator*=(Preon::Math::float_simd& a, Preon::Math::float_simd b) { return a = a * b; }
PREONMATH_FORCEINLINE Preon::Math::float_simd& operator/=(Preon::Math::float_simd& a, Preon::Math::float_simd b) { return a = a / b; }

PREONMATH_FORCEINLINE Preon::Math::float_simd operator>(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_cmpgt_ps(a, b); }
PREONMATH_FORCEINLINE Preon::Math::float_simd operator>=(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_cmpge_ps(a, b); }
PREONMATH_FORCEINLINE Preon::Math::float_simd operator<(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_cmplt_ps(a, b); }
PREONMATH_FORCEINLINE Preon::Math::float_simd operator<=(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_cmple_ps(a, b); }
PREONMATH_FORCEINLINE Preon::Math::float_simd operator==(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_cmpeq_ps(a, b); }
PREONMATH_FORCEINLINE Preon::Math::float_simd operator!=(Preon::Math::float_simd a, Preon::Math::float_simd b) { return _mm_cmpneq_ps(a, b); }

PREONMATH_FORCEINLINE Preon::Math::float_simd operator-(Preon::Math::float_simd a) { return Preon::Math::Simd::xor_mask(a, Preon::Math::Simd::make<float>(-0.f)); }

PREONMATH_FORCEINLINE Preon::Math::double_simd operator+(Preon::Math::double_simd a, Preon::Math::double_simd b) { return _mm_add_pd(a, b); }
PREONMATH_FORCEINLINE Preon::Math::double_simd operator-(Preon::Math::double_simd a, Preon::Math::double_simd b) { return _mm_sub_pd(a, b); }
PREONMATH_FORCEINLINE Preon::Math::double_simd operator*(Preon::Math::double_simd a, Preon::Math::double_simd b) { return _mm_mul_pd(a, b); }
PREONMATH_FORCEINLINE Preon::Math::double_simd operator/(Preon::Math::double_simd a, Preon::Math::double_simd b) { return _mm_div_pd(a, b); }

PREONMATH_FORCEINLINE Preon::Math::double_simd& operator+=(Preon::Math::double_simd& a, Preon::Math::double_simd b) { return a = a + b; }
PREONMATH_FORCEINLINE Preon::Math::double_simd& operator-=(Preon::Math::double_simd& a, Preon::Math::double_simd b) { return a = a - b; }
PREONMATH_FORCEINLINE Preon::Math::double_simd& operator*=(Preon::Math::double_simd& a, Preon::Math::double_simd b) { return a = a * b; }
PREONMATH_FORCEINLINE Preon::Math::double_simd& operator/=(Preon::Math::double_simd& a, Preon::Math::double_simd b) { return a = a / b; }

PREONMATH_FORCEINLINE Preon::Math::double_simd operator>(Preon::Math::double_simd a, Preon::Math::double_simd b) { return _mm_cmpgt_pd(a, b); }
PREONMATH_FORCEINLINE Preon::Math::double_simd operator>=(Preon::Math::double_simd a, Preon::Math::double_simd b) { return _mm_cmpge_pd(a, b); }
PREONMATH_FORCEINLINE Preon::Math::double_simd operator<(Preon::Math::double_simd a, Preon::Math::double_simd b) { return _mm_cmplt_pd(a, b); }
PREONMATH_FORCEINLINE Preon::Math::double_simd operator<=(Preon::Math::double_simd a, Preon::Math::double_simd b) { return _mm_cmple_pd(a, b); }
PREONMATH_FORCEINLINE Preon::Math::double_simd operator==(Preon::Math::double_simd a, Preon::Math::double_simd b) { return _mm_cmpeq_pd(a, b); }

PREONMATH_FORCEINLINE Preon::Math::double_simd operator-(Preon::Math::double_simd a) { return Preon::Math::Simd::xor_mask(a, Preon::Math::Simd::make<double>(-0.0)); }
#endif  // PREONMATH_COMPILER_MSVC

#endif  // PREONMATH_SCALAR_SSE_H
