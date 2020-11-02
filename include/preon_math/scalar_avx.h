/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "basics.h"

#ifdef PREONMATH_COMPILER_MSVC
    #include <intrin.h>
#else
    #include <pmmintrin.h>
#endif

namespace Preon
{
namespace Math
{
    using float_simd = __m256;
    using double_simd = __m256d;

    namespace Simd
    {
        template<typename scalar>
        using simd_if_float = typename std::enable_if<!std::is_same<scalar, double>::value, float_simd>::type;
        template<typename scalar>
        using simd_if_double = typename std::enable_if<std::is_same<scalar, double>::value, double_simd>::type;

        template<typename T>
        struct Register
        {
        };
        template<>
        struct Register<float>
        {
            enum
            {
                size = 8
            };
            enum
            {
                expo = 3
            };
            using type = float_simd;
        };
        template<>
        struct Register<double>
        {
            enum
            {
                size = 4
            };
            enum
            {
                expo = 2
            };
            using type = double_simd;
        };

        // We demand explicit float / double typing here, because otherwise it would be too easy to create the wrong type by accident.
        template<typename scalar>
        PREONMATH_FORCEINLINE simd_if_float<scalar> make(float a)
        {
            return _mm256_set1_ps(a);
        }
        template<typename scalar>
        PREONMATH_FORCEINLINE simd_if_double<scalar> make(double a)
        {
            return _mm256_set1_pd(a);
        }

        PREONMATH_FORCEINLINE float_simd load(const float* values) { return _mm256_loadu_ps(values); }
        PREONMATH_FORCEINLINE double_simd load(const double* values) { return _mm256_loadu_pd(values); }

        PREONMATH_FORCEINLINE void store(float_simd input, float* values) { _mm256_storeu_ps(values, input); }
        PREONMATH_FORCEINLINE void store(double_simd input, double* values) { _mm256_storeu_pd(values, input); }

        template<typename scalar>
        PREONMATH_FORCEINLINE simd_if_float<scalar> zero()
        {
            return _mm256_setzero_ps();
        }

        template<typename scalar>
        PREONMATH_FORCEINLINE simd_if_double<scalar> zero()
        {
            return _mm256_setzero_pd();
        }

        PREONMATH_FORCEINLINE float_simd and_mask(float_simd a, float_simd b) { return _mm256_and_ps(a, b); }
        PREONMATH_FORCEINLINE float_simd or_mask(float_simd a, float_simd b) { return _mm256_or_ps(a, b); }
        PREONMATH_FORCEINLINE float_simd xor_mask(float_simd a, float_simd b) { return _mm256_xor_ps(a, b); }
        PREONMATH_FORCEINLINE float_simd blendv(float_simd a, float_simd b, float_simd mask) { return _mm256_blendv_ps(a, b, mask); }

        PREONMATH_FORCEINLINE double_simd and_mask(double_simd a, double_simd b) { return _mm256_and_pd(a, b); }
        PREONMATH_FORCEINLINE double_simd or_mask(double_simd a, double_simd b) { return _mm256_or_pd(a, b); }
        PREONMATH_FORCEINLINE double_simd xor_mask(double_simd a, double_simd b) { return _mm256_xor_pd(a, b); }
        PREONMATH_FORCEINLINE double_simd blendv(double_simd a, double_simd b, double_simd mask) { return _mm256_blendv_pd(a, b, mask); }

        template<class Getter>
        PREONMATH_FORCEINLINE void setScalars(float_simd* out, const Getter& func)
        {
            (*out) = _mm256_setr_ps(func(0), func(1), func(2), func(3), func(4), func(5), func(6), func(7));
        }
        template<class Getter>
        PREONMATH_FORCEINLINE void setScalars(double_simd* out, const Getter& func)
        {
            (*out) = _mm256_setr_pd(func(0), func(1), func(2), func(3));
        }

        template<class Getter>
        PREONMATH_FORCEINLINE void setScalars(float_simd* out, const Getter& func, uint numElements, float fillValue)
        {
            (*out) = _mm256_setr_ps(
                func(0),
                (numElements > 1) ? func(1) : fillValue,
                (numElements > 2) ? func(2) : fillValue,
                (numElements > 3) ? func(3) : fillValue,
                (numElements > 4) ? func(4) : fillValue,
                (numElements > 5) ? func(5) : fillValue,
                (numElements > 6) ? func(6) : fillValue,
                (numElements > 7) ? func(7) : fillValue);
        }
        template<class Getter>
        PREONMATH_FORCEINLINE void setScalars(double_simd* out, const Getter& func, uint numElements, double fillValue)
        {
            (*out) = _mm256_setr_pd(func(0), (numElements > 1) ? func(1) : fillValue, (numElements > 2) ? func(2) : fillValue, (numElements > 3) ? func(3) : fillValue);
        }

        template<class Getter>
        PREONMATH_FORCEINLINE void storeScalars(float_simd input, const Getter& func)
        {
            float values[8];
            store(input, values);
            for (int i = 0; i < 8; i++)
                func(i) = values[i];
        }
        template<class Getter>
        PREONMATH_FORCEINLINE void storeScalars(double_simd input, const Getter& func)
        {
            double values[4];
            store(input, values);
            for (int i = 0; i < 4; i++)
                func(i) = values[i];
        }

        PREONMATH_FORCEINLINE float hSum(float_simd v)
        {
            float out;
            __m128 x128 = _mm_add_ps(_mm256_extractf128_ps(v, 1), _mm256_castps256_ps128(v));
            x128 = _mm_hadd_ps(x128, x128);
            x128 = _mm_hadd_ps(x128, x128);
            _mm_store_ss(&out, x128);
            return out;
        }
        PREONMATH_FORCEINLINE double hSum(double_simd v)
        {
            double out;
            double_simd x = _mm_hadd_pd(v, v);
            x = _mm_hadd_pd(x, x);
            _mm_store_sd(&out, x);
            return out;
        }

        PREONMATH_FORCEINLINE bool _testAllZeroFloat(const float_simd x) { return static_cast<bool>(_mm256_test_all_zeros(_mm256_castps_si256(x), _mm256_castps_si256(x))); }

        template<typename T>
        PREONMATH_FORCEINLINE bool testAllZero(const T x)
        {
            return _testAllZeroFloat((float_simd)x);
        }

        PREONMATH_FORCEINLINE float_simd round(const float_simd x) { return _mm256_round_ps(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); }
        PREONMATH_FORCEINLINE double_simd round(const double_simd x) { return _mm256_round_pd(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); }
    }  // namespace Simd

    using preal_simd = typename Simd::Register<preal>::type;

}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_COMPILER_MSVC  // Operators are already defined on gcc / clang
PREONMATH_FORCEINLINE float_simd operator+(float_simd a, float_simd b)
{
    return _mm256_add_ps(a, b);
}
PREONMATH_FORCEINLINE float_simd operator-(float_simd a, float_simd b)
{
    return _mm256_sub_ps(a, b);
}
PREONMATH_FORCEINLINE float_simd operator*(float_simd a, float_simd b)
{
    return _mm256_mul_ps(a, b);
}
PREONMATH_FORCEINLINE float_simd operator/(float_simd a, float_simd b)
{
    return _mm256_div_ps(a, b);
}

PREONMATH_FORCEINLINE float_simd& operator+=(float_simd& a, float_simd b)
{
    return a = a + b;
}
PREONMATH_FORCEINLINE float_simd& operator-=(float_simd& a, float_simd b)
{
    return a = a - b;
}
PREONMATH_FORCEINLINE float_simd& operator*=(float_simd& a, float_simd b)
{
    return a = a * b;
}
PREONMATH_FORCEINLINE float_simd& operator/=(float_simd& a, float_simd b)
{
    return a = a / b;
}

PREONMATH_FORCEINLINE float_simd operator>(float_simd a, float_simd b)
{
    return _mm256_cmpgt_ps(a, b);
}
PREONMATH_FORCEINLINE float_simd operator>=(float_simd a, float_simd b)
{
    return _mm256_cmpge_ps(a, b);
}
PREONMATH_FORCEINLINE float_simd operator<(float_simd a, float_simd b)
{
    return _mm256_cmplt_ps(a, b);
}
PREONMATH_FORCEINLINE float_simd operator<=(float_simd a, float_simd b)
{
    return _mm256_cmple_ps(a, b);
}
PREONMATH_FORCEINLINE float_simd operator==(float_simd a, float_simd b)
{
    return _mm256_cmpeq_ps(a, b);
}

PREONMATH_FORCEINLINE float_simd operator-(float_simd a)
{
    return Simd::xor_mask(a, Simd::make<float>(-0.f));
}

PREONMATH_FORCEINLINE double_simd operator+(double_simd a, double_simd b)
{
    return _mm256_add_pd(a, b);
}
PREONMATH_FORCEINLINE double_simd operator-(double_simd a, double_simd b)
{
    return _mm256_sub_pd(a, b);
}
PREONMATH_FORCEINLINE double_simd operator*(double_simd a, double_simd b)
{
    return _mm256_mul_pd(a, b);
}
PREONMATH_FORCEINLINE double_simd operator/(double_simd a, double_simd b)
{
    return _mm256_div_pd(a, b);
}

PREONMATH_FORCEINLINE double_simd& operator+=(double_simd& a, double_simd b)
{
    return a = a + b;
}
PREONMATH_FORCEINLINE double_simd& operator-=(double_simd& a, double_simd b)
{
    return a = a - b;
}
PREONMATH_FORCEINLINE double_simd& operator*=(double_simd& a, double_simd b)
{
    return a = a * b;
}
PREONMATH_FORCEINLINE double_simd& operator/=(double_simd& a, double_simd b)
{
    return a = a / b;
}

PREONMATH_FORCEINLINE double_simd operator>(double_simd a, double_simd b)
{
    return _mm256_cmpgt_pd(a, b);
}
PREONMATH_FORCEINLINE double_simd operator>=(double_simd a, double_simd b)
{
    return _mm256_cmpge_pd(a, b);
}
PREONMATH_FORCEINLINE double_simd operator<(double_simd a, double_simd b)
{
    return _mm256_cmplt_pd(a, b);
}
PREONMATH_FORCEINLINE double_simd operator<=(double_simd a, double_simd b)
{
    return _mm256_cmple_pd(a, b);
}
PREONMATH_FORCEINLINE double_simd operator==(double_simd a, double_simd b)
{
    return _mm256_cmpeq_pd(a, b);
}

PREONMATH_FORCEINLINE double_simd operator-(double_simd a)
{
    return Simd::xor_mask(a, Simd::make<double>(-0.0));
}
#endif

namespace std
{
PREONMATH_FORCEINLINE float_simd max(float_simd a, float_simd b)
{
    return _mm256_max_ps(a, b);
}
PREONMATH_FORCEINLINE float_simd min(float_simd a, float_simd b)
{
    return _mm256_min_ps(a, b);
}
PREONMATH_FORCEINLINE float_simd sqrt(float_simd a)
{
    return _mm256_sqrt_ps(a);
}
PREONMATH_FORCEINLINE float_simd abs(float_simd a)
{
    const float_simd signMask = _mm256_set1_ps(-0.f);
    return _mm256_andnot_ps(signMask, a);
}

PREONMATH_FORCEINLINE double_simd max(double_simd a, double_simd b)
{
    return _mm256_max_pd(a, b);
}
PREONMATH_FORCEINLINE double_simd min(double_simd a, double_simd b)
{
    return _mm256_min_pd(a, b);
}
PREONMATH_FORCEINLINE double_simd sqrt(double_simd a)
{
    return _mm256_sqrt_pd(a);
}
PREONMATH_FORCEINLINE double_simd abs(double_simd a)
{
    const double_simd signMask = _mm256_set1_pd(-0.0);
    return _mm256_andnot_pd(signMask, a);
}
}  // namespace std
