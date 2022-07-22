/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "vec.h"
#include "scalar_simd.h"
#include "vec.h"
#include "vec3.h"
#include "static_for.h"

#ifdef PREONMATH_COMPILER_MSVC
    #include <intrin.h>
#else
    #include <immintrin.h>
#endif

#include <utility>

namespace Preon
{
namespace Math
{
    template<size_t D, typename T = float>
    using vec_simd = vec<D, typename Simd::Register<T>::type>;

    namespace Simd
    {
#ifdef USE_AVX
        // =======================================================================================
        // ================================== AVX specific =======================================
        // =======================================================================================

        // =======================================================================================
        // ==================================== Load =============================================
        // =======================================================================================

        template<size_t D>
        PREONMATH_FORCEINLINE void setVecs(
            typename std::enable_if<D != 3, vec<3, float>>::type* out,
            const vec<D, float>& v0,
            const vec<D, float>& v1,
            const vec<D, float>& v2,
            const vec<D, float>& v3,
            const vec<D, float>& v4,
            const vec<D, float>& v5,
            const vec<D, float>& v6,
            const vec<D, float>& v7)
        {
            StaticFor<0, D>([&](size_t d) { (*out)[d] = _mm256_setr_ps(v0[d], v1[d], v2[d], v3[d], v4[d], v5[d], v6[d], v7[d]); });
        }

        //! Vec3f optimized version.
        PREONMATH_FORCEINLINE void setVecs(
            vec_simd<3>* out,
            const vec<3, float>& v0,
            const vec<3, float>& v1,
            const vec<3, float>& v2,
            const vec<3, float>& v3,
            const vec<3, float>& v4,
            const vec<3, float>& v5,
            const vec<3, float>& v6,
            const vec<3, float>& v7)
        {
            __m128 xy0, xy1;
            xy0 = xy1 = _mm_setzero_ps();
            xy0 = _mm_loadl_pi(xy0, (const __m64*)&v0[0]);
            xy1 = _mm_loadl_pi(xy1, (const __m64*)&v2[0]);
            xy0 = _mm_loadh_pi(xy0, (const __m64*)&v1[0]);
            xy1 = _mm_loadh_pi(xy1, (const __m64*)&v3[0]);
            (*out)[0] = _mm256_castps128_ps256(_mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(2, 0, 2, 0)));
            (*out)[1] = _mm256_castps128_ps256(_mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(3, 1, 3, 1)));

            xy0 = _mm_loadl_pi(xy0, (const __m64*)&v4[0]);
            xy1 = _mm_loadl_pi(xy1, (const __m64*)&v6[0]);
            xy0 = _mm_loadh_pi(xy0, (const __m64*)&v5[0]);
            xy1 = _mm_loadh_pi(xy1, (const __m64*)&v7[0]);
            (*out)[0] = _mm256_insertf128_ps((*out)[0], _mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(2, 0, 2, 0)), 1);
            (*out)[1] = _mm256_insertf128_ps((*out)[1], _mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(3, 1, 3, 1)), 1);

            (*out)[2] = _mm256_setr_ps(v0[2], v1[2], v2[2], v3[2], v4[2], v5[2], v6[2], v7[2]);
        }

        PREONMATH_FORCEINLINE void setVecs_1(vec_simd<3>* out, const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3)
        {
            __m128 xy0, xy1;
            xy0 = xy1 = _mm_setzero_ps();
            xy0 = _mm_loadl_pi(xy0, (const __m64*)&v0[0]);
            xy1 = _mm_loadl_pi(xy1, (const __m64*)&v2[0]);
            xy0 = _mm_loadh_pi(xy0, (const __m64*)&v1[0]);
            xy1 = _mm_loadh_pi(xy1, (const __m64*)&v3[0]);
            (*out)[0] = _mm256_castps128_ps256(_mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(2, 0, 2, 0)));
            (*out)[1] = _mm256_castps128_ps256(_mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(3, 1, 3, 1)));

            (*out)[2] = _mm256_castps128_ps256(_mm_setr_ps(v0[2], v1[2], v2[2], v3[2]));
        }

        PREONMATH_FORCEINLINE void setVecs_2(vec_simd<3>* out, const vec3f& v4, const vec3f& v5, const vec3f& v6, const vec3f& v7)
        {
            __m128 xy0, xy1;
            xy0 = xy1 = _mm_setzero_ps();
            xy0 = _mm_loadl_pi(xy0, (const __m64*)&v4[0]);
            xy1 = _mm_loadl_pi(xy1, (const __m64*)&v6[0]);
            xy0 = _mm_loadh_pi(xy0, (const __m64*)&v5[0]);
            xy1 = _mm_loadh_pi(xy1, (const __m64*)&v7[0]);
            (*out)[0] = _mm256_insertf128_ps((*out)[0], _mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(2, 0, 2, 0)), 1);
            (*out)[1] = _mm256_insertf128_ps((*out)[1], _mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(3, 1, 3, 1)), 1);

            (*out)[2] = _mm256_insertf128_ps((*out)[2], _mm_setr_ps(v4[2], v5[2], v6[2], v7[2]), 1);
        }

        template<class Getter>
        PREONMATH_FORCEINLINE void setVecs(vec_simd<3>* out, const Getter& func)
        {
            setVecs(out, func(0), func(1), func(2), func(3), func(4), func(5), func(6), func(7));
        }

        template<class Getter>
        PREONMATH_FORCEINLINE void setVecs(vec_simd<3>* out, const Getter& func, uint32_t numElements, const vec<3, float>& fillValue)
        {
            if (numElements > 4)
            {
                setVecs_1(out, func(0), func(1), func(2), func(3));
                setVecs_2(out, func(4), (numElements > 5) ? func(5) : fillValue, (numElements > 6) ? func(6) : fillValue, (numElements > 7) ? func(7) : fillValue);
            }
            else
            {
                setVecs_1(out, func(0), (numElements > 1) ? func(1) : fillValue, (numElements > 2) ? func(2) : fillValue, (numElements > 3) ? func(3) : fillValue);
                StaticFor<0, 3>([&](size_t d) { (*out)[d] = _mm256_insertf128_ps((*out)[d], _mm_broadcast_ss(&fillValue[d]), 1); });
            }
        }

        // See https://software.intel.com/en-us/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx.
        PREONMATH_FORCEINLINE vec_simd<3, float> load(const float* vecs)
        {
            __m128* m = (__m128*)vecs;
            __m256 m03;
            __m256 m14;
            __m256 m25;
            m03 = _mm256_castps128_ps256(m[0]);  // load lower halves
            m14 = _mm256_castps128_ps256(m[1]);
            m25 = _mm256_castps128_ps256(m[2]);
            m03 = _mm256_insertf128_ps(m03, m[3], 1);  // load upper halves
            m14 = _mm256_insertf128_ps(m14, m[4], 1);
            m25 = _mm256_insertf128_ps(m25, m[5], 1);

            __m256 xy = _mm256_shuffle_ps(m14, m25, _MM_SHUFFLE(2, 1, 3, 2));  // upper x's and y's
            __m256 yz = _mm256_shuffle_ps(m03, m14, _MM_SHUFFLE(1, 0, 2, 1));  // lower y's and z's
            vec_simd<3, float> out;
            (*out)[0] = _mm256_shuffle_ps(m03, xy, _MM_SHUFFLE(2, 0, 3, 0));
            (*out)[1] = _mm256_shuffle_ps(yz, xy, _MM_SHUFFLE(3, 1, 2, 0));
            (*out)[2] = _mm256_shuffle_ps(yz, m25, _MM_SHUFFLE(3, 0, 3, 1));
            return out;
        }

        // =======================================================================================
        // ==================================== Store ============================================
        // =======================================================================================

        template<size_t D>
        PREONMATH_FORCEINLINE void getVecs(
            const vec_simd<D>& input,
            vec<D, float>& v0,
            vec<D, float>& v1,
            vec<D, float>& v2,
            vec<D, float>& v3,
            vec<D, float>& v4,
            vec<D, float>& v5,
            vec<D, float>& v6,
            vec<D, float>& v7)
        {
            float buffer[8];
            StaticFor<0, D>([&](size_t d) {
                _mm256_storeu_ps(buffer, input[d]);
                v0[d] = buffer[0];
                v1[d] = buffer[1];
                v2[d] = buffer[2];
                v3[d] = buffer[3];
                v4[d] = buffer[4];
                v5[d] = buffer[5];
                v6[d] = buffer[6];
                v7[d] = buffer[7];
            });
        }

        template<size_t D, class Getter>
        PREONMATH_FORCEINLINE void storeVecs(const vec_simd<D>& input, const Getter& func)
        {
            getVecs<D>(input, func(0), func(1), func(2), func(3), func(4), func(5), func(6), func(7));
        }

        // =======================================================================================
        // ==================================== Horizontal operators =============================
        // =======================================================================================

        //! For each register, the 4 values are summed up.
        template<size_t D>
        PREONMATH_FORCEINLINE typename std::enable_if<D != 3, vec<D, float>>::type hSum(const vec_simd<D>& v)
        {
            vec<D, float> out;
            StaticFor<0, D>([&](int d) { out[d] = hSum(v[d]); });
            return out;
        }

        //! Optimized version of horizontal add for 3D-vectors. Performs four _mm_hadd_ps instead of six.
        template<size_t D>
        PREONMATH_FORCEINLINE typename std::enable_if<D == 3, vec3f>::type hSum(const vec_simd<D>& v)
        {
            vec3f out;

            // extraxt x and y component
            __m128 tmp1 = _mm_add_ps(_mm256_extractf128_ps(v[0], 1), _mm256_castps256_ps128(v[0]));
            __m128 tmp2 = _mm_add_ps(_mm256_extractf128_ps(v[1], 1), _mm256_castps256_ps128(v[1]));
            tmp1 = _mm_hadd_ps(tmp1, tmp2);
            tmp1 = _mm_hadd_ps(tmp1, tmp1);
            _mm_storel_pi((__m64*)&out[0], tmp1);

            // extract z
            tmp2 = _mm_add_ps(_mm256_extractf128_ps(v[2], 1), _mm256_castps256_ps128(v[2]));
            tmp2 = _mm_hadd_ps(tmp2, tmp2);
            tmp2 = _mm_hadd_ps(tmp2, tmp2);
            _mm_store_ss(&out[2], tmp2);

            return out;
        }
#else
        // =======================================================================================
        // ================================== SSE specific =======================================
        // =======================================================================================

        // =======================================================================================
        // ==================================== Load =============================================
        // =======================================================================================

        template<>
        struct Register<vec3f>  // Todo: Is it possible to to make this generic for other dimensions?
        {
            using type = vec_simd<3, float>;
        };
        template<>
        struct Register<vec3d>
        {
            using type = vec_simd<3, double>;
        };

        template<typename V>
        PREONMATH_FORCEINLINE typename std::enable_if<std::is_same<V, vec3d>::value, vec<3, double_simd>>::type zero()
        {
            return vec<3, double_simd>::zero();
        }

        template<typename V>
        PREONMATH_FORCEINLINE typename std::enable_if<std::is_same<V, vec3f>::value, vec<3, float_simd>>::type zero()
        {
            return vec<3, float_simd>::zero();
        }

        template<size_t D>
        PREONMATH_FORCEINLINE void
        setVecs(typename std::enable_if<D != 3, vec_simd<D, float>>::type* out, const vec<D, float>& v0, const vec<D, float>& v1, const vec<D, float>& v2, const vec<D, float>& v3)
        {
            StaticFor<0, D>([&](size_t d) { (*out)[d] = _mm_setr_ps(v0[d], v1[d], v2[d], v3[d]); });
        }

        //! Vec3f optimized version.
        PREONMATH_FORCEINLINE void setVecs(vec_simd<3, float>* out, const vec<3, float>& v0, const vec<3, float>& v1, const vec<3, float>& v2, const vec<3, float>& v3)
        {
            __m128 xy0 = _mm_castsi128_ps(_mm_loadl_epi64((const __m128i*)&v0[0]));
            __m128 xy1 = _mm_castsi128_ps(_mm_loadl_epi64((const __m128i*)&v2[0]));
            xy0 = _mm_loadh_pi(xy0, (const __m64*)&v1[0]);
            xy1 = _mm_loadh_pi(xy1, (const __m64*)&v3[0]);
            (*out)[0] = _mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(2, 0, 2, 0));
            (*out)[1] = _mm_shuffle_ps(xy0, xy1, _MM_SHUFFLE(3, 1, 3, 1));
            (*out)[2] = _mm_setr_ps(v0[2], v1[2], v2[2], v3[2]);
        }

        //! Vec4f optimized version.
        PREONMATH_FORCEINLINE void setVecs(vec_simd<4, float>* out, const vec<4, float>& v0, const vec<4, float>& v1, const vec<4, float>& v2, const vec<4, float>& v3)
        {
            (*out)[0] = _mm_loadu_ps((const float*)&v0);
            (*out)[1] = _mm_loadu_ps((const float*)&v1);
            (*out)[2] = _mm_loadu_ps((const float*)&v2);
            (*out)[3] = _mm_loadu_ps((const float*)&v3);

            _MM_TRANSPOSE4_PS((*out)[0], (*out)[1], (*out)[2], (*out)[3]);
        }

        template<class Getter, size_t D>
        PREONMATH_FORCEINLINE void setVecs(vec_simd<D, float>* out, const Getter& func)
        {
            setVecs(out, func(0), func(1), func(2), func(3));
        }

        template<class Getter, size_t D>
        PREONMATH_FORCEINLINE void setVecs(vec_simd<D, float>* out, const Getter& func, uint32_t numElements, const vec<D, float>& fillValue)
        {
            setVecs(out, func(0), (numElements > 1) ? func(1) : fillValue, (numElements > 2) ? func(2) : fillValue, (numElements > 3) ? func(3) : fillValue);
        }

        // See https://software.intel.com/en-us/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx.
        PREONMATH_FORCEINLINE vec_simd<3, float> load(const vec3f* vecs)
        {
            __m128 x0y0z0x1 = _mm_loadu_ps((float*)vecs);
            __m128 y1z1x2y2 = _mm_loadu_ps(((float*)vecs) + 4);
            __m128 z2x3y3z3 = _mm_loadu_ps(((float*)vecs) + 8);
            __m128 x2y2x3y3 = _mm_shuffle_ps(y1z1x2y2, z2x3y3z3, _MM_SHUFFLE(2, 1, 3, 2));
            __m128 y0z0y1z1 = _mm_shuffle_ps(x0y0z0x1, y1z1x2y2, _MM_SHUFFLE(1, 0, 2, 1));
            return vec_simd<3, float>(
                _mm_shuffle_ps(x0y0z0x1, x2y2x3y3, _MM_SHUFFLE(2, 0, 3, 0)),
                _mm_shuffle_ps(y0z0y1z1, x2y2x3y3, _MM_SHUFFLE(3, 1, 2, 0)),
                _mm_shuffle_ps(y0z0y1z1, z2x3y3z3, _MM_SHUFFLE(3, 0, 3, 1)));
        }

        //! Double loaders. No fancy optimizations here for now.
        template<size_t D>
        PREONMATH_FORCEINLINE void setVecs(vec_simd<3, double>* out, const vec<D, double>& v0, const vec<D, double>& v1)
        {
            StaticFor<0, D>([&](size_t d) { (*out)[d] = _mm_setr_pd(v0[d], v1[d]); });
        }

        template<class Getter>
        PREONMATH_FORCEINLINE void setVecs(vec_simd<3, double>* out, const Getter& func)
        {
            setVecs(out, func(0), func(1));
        }

        template<class Getter>
        PREONMATH_FORCEINLINE void setVecs(vec_simd<3, double>* out, const Getter& func, uint32_t numElements, const vec<3, double>& fillValue)
        {
            setVecs(out, func(0), (numElements > 1) ? func(1) : fillValue);
        }

        PREONMATH_FORCEINLINE vec_simd<3, double> load(const vec3d* vecs)
        {
            vec_simd<3, double> out;
            setVecs(&out, [vecs](int i) { return vecs[i]; });
            return out;
        }

        // =======================================================================================
        // ==================================== Store ============================================
        // =======================================================================================

        template<size_t D>
        PREONMATH_FORCEINLINE void getVecs(const vec_simd<D, float>& input, vec<D, float>& v0, vec<D, float>& v1, vec<D, float>& v2, vec<D, float>& v3)
        {
            float buffer[4];
            StaticFor<0, D>([&](size_t d) {
                _mm_storeu_ps(buffer, input[d]);
                v0[d] = buffer[0];
                v1[d] = buffer[1];
                v2[d] = buffer[2];
                v3[d] = buffer[3];
            });
        }

        template<size_t D, class Getter>
        PREONMATH_FORCEINLINE void storeVecs(const vec_simd<D, float>& input, const Getter& func)
        {
            getVecs<D>(input, func(0), func(1), func(2), func(3));
        }

        template<size_t D>
        PREONMATH_FORCEINLINE void getVecs(const vec_simd<D, double>& input, vec<D, double>& v0, vec<D, double>& v1)
        {
            double buffer[2];
            for (size_t d = 0; d < D; ++d)
            {
                _mm_storeu_pd(buffer, input[d]);
                v0[d] = buffer[0];
                v1[d] = buffer[1];
            }
        }

        template<size_t D, class Getter>
        PREONMATH_FORCEINLINE void storeVecs(const vec_simd<D, double>& input, const Getter& func)
        {
            getVecs<D>(input, func(0), func(1));
        }

        // =======================================================================================
        // ==================================== Horizontal operators =============================
        // =======================================================================================

        //! For each register, the 4 values are summed up.
        template<size_t D>
        PREONMATH_FORCEINLINE typename std::enable_if<D != 3, vec<D, float>>::type hSum(const vec_simd<D, float>& v)
        {
            vec<D, float> out;
            StaticFor<0, D>([&](int d) { out[d] = hSum(v[d]); });
            return out;
        }

        //! Optimized version of horizontal add for 3D-vectors. Performs four _mm_hadd_ps instead of six.
        template<size_t D>
        PREONMATH_FORCEINLINE typename std::enable_if<D == 3, vec3f>::type hSum(const vec_simd<D, float>& v)
        {
            vec3f out;
            __m128 xy = _mm_hadd_ps(v[0], v[1]);
            __m128 z = _mm_hadd_ps(v[2], v[2]);
            xy = _mm_hadd_ps(xy, xy);
            z = _mm_hadd_ps(z, z);
            _mm_storel_pi((__m64*)&out[0], xy);
            _mm_store_ss(&out[2], z);
            return out;
        }

        template<size_t D>
        PREONMATH_FORCEINLINE vec<D, double> hSum(const vec_simd<D, double>& v)
        {
            vec<D, double> out;
            for (size_t d = 0; d < D; ++d)
                out[d] = hSum(v[d]);
            return out;
        }

#endif
        // =======================================================================================
        // ==================================== SSE / AVX agnostic ===============================
        // =======================================================================================

        template<size_t D, typename T>
        PREONMATH_FORCEINLINE void load(vec_simd<D, T>* out, const T* values, int d)
        {
            (*out)[d] = load(values);
        }

        template<size_t D, typename T>
        PREONMATH_FORCEINLINE void store(const vec_simd<D, T>& input, T* values, int d)
        {
            store(values, input[d]);
        }

        template<typename T>
        PREONMATH_FORCEINLINE vec<Simd::Register<T>::size, T> simdScalarToVec(typename Simd::Register<T>::type input)
        {
            vec<Simd::Register<T>::size, T> out;
            store(input, &out[0]);
            return out;
        }

        //! Performs mask[d] ? v1[d] : v2[d] for each component and returns the result.
        template<size_t D, typename T>
        PREONMATH_FORCEINLINE vec_simd<D, T> choose(const vec_simd<D, T>& mask, const vec_simd<D, T>& v1, const vec_simd<D, T>& v2)
        {
            vec_simd<D, T> out;
            StaticFor<0, D>([&mask, &v1, &v2, &out](int d) { out[d] = blendv(v2[d], v1[d], mask[d]); });
            return out;
        }

        //! Performs mask[d] ? v1[d] : v2[d] for each component and returns the result.
        template<size_t D, typename T>
        PREONMATH_FORCEINLINE vec_simd<D, T> choose(T mask, const vec_simd<D, T>& v1, const vec_simd<D, T>& v2)
        {
            vec_simd<D, T> out;
            StaticFor<0, D>([&mask, &v1, &v2, &out](int d) { out[d] = blendv(v2[d], v1[d], mask); });
            return out;
        }

        //! Performs mask[d] ? vec[d] : 0 for each component and returns the result.
        template<size_t D>
        PREONMATH_FORCEINLINE vec_simd<D, float> and_mask(const vec_simd<D, float>& vec, float_simd mask)
        {
            return vec_simd<D, float>([&](size_t d) { return Simd::and_mask(vec[d], mask); });
        }

        // Todo: Implement more efficient horizontal max and min.
        template<size_t D, typename T>
        PREONMATH_FORCEINLINE vec<D, T> hMax(const vec_simd<D, T>& v)
        {
            vec<D, T> out;
            StaticFor<0, D>([&](int d) { out[d] = hMax(v[d]); });
            return out;
        }

        template<size_t D, typename T>
        PREONMATH_FORCEINLINE vec<D, T> hMin(const vec_simd<D, T>& v)
        {
            vec<D, T> out;
            StaticFor<0, D>([&](int d) { out[d] = hMin(v[d]); });
            return out;
        }
    }  // namespace Simd
}  // namespace Math
}  // namespace Preon
