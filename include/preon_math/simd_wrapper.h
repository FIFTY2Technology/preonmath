/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"
#include "scalar_simd.h"
#include "vec_simd.h"
#include "vec.h"
#include "matrix.h"
#include "matrix_simd.h"
#include "simd_scalar_wrapper_definitions.h"

#include <type_traits>
namespace Preon
{
namespace Math
{
    namespace SimdWrapper
    {
        template<typename T>
        struct Scalar
        {
            using type = typename Simd::Scalar<T>::type;
        };

        // scalar make wrapper
        template<>
        inline float_simd make<float_simd>(float val)
        {
            return Simd::make<float>(val);
        }

        template<>
        inline double_simd make<double_simd>(double val)
        {
            return Simd::make<double>(val);
        }
        // also overload the following for convenience
        template<>
        inline float_simd make<float_simd>(double val)
        {
            return Simd::make<float>(val);
        }

        // vector make wrapper
        template<>
        inline vec_simd<3> make<vec_simd<3>>(const vec3f& val)
        {
            return vec_simd<3>(val);
        }

        template<>
        inline float_simd zero<float_simd>()
        {
            return Simd::zero<float>();
        }
        template<>
        inline double_simd zero<double_simd>()
        {
            return Simd::zero<double>();
        }

        // scalar hSum wrapper
        inline float hSum(const float_simd& val)
        {
            return Simd::hSum(val);
        }

        // vector hSum wrapper
        inline vec3f hSum(const vec<3, float_simd>& val)
        {
            return Simd::hSum(val);
        }

        // matrix hSum wrapper
        inline matrix33f hSum(const matrix_simd<3, 3, float>& val)
        {
            return Simd::hSum(val);
        }

        // choose
        inline float_simd choose(const float_simd mask, const float_simd v1, const float_simd v2)
        {
            return Simd::choose(mask, v1, v2);
        }

        template<PrMathSize D>
        inline vec<D, float_simd> choose(const float_simd mask, const vec<D, float_simd>& v1, const vec<D, float_simd>& v2)
        {
            return vec_simd<D>::choose(mask, v1, v2);
        }

        // AND mask
        inline float_simd and_mask(float_simd a, float_simd b)
        {
            return Simd::and_mask(a, b);
        }

        // max
        inline float_simd max(float_simd a, float_simd b)
        {
            return Simd::max(a, b);
        }

        // min
        inline float_simd min(float_simd a, float_simd b)
        {
            return Simd::min(a, b);
        }

        // cosine
        inline float_simd cos(float_simd val)
        {
            // TODO: this is not very efficient
            std::array<float, Simd::Register<float>::size> vals;
            Simd::store(val, vals.data());
            std::for_each(vals.begin(), vals.end(), [](float& v) { v = cos(v); });
            return Simd::load(vals.data());
        }
        // exp
        inline float_simd exp(float_simd val)
        {
            // TODO: this is not very efficient
            std::array<float, Simd::Register<float>::size> vals;
            Simd::store(val, vals.data());
            std::for_each(vals.begin(), vals.end(), [](float& v) { v = exp(v); });
            return Simd::load(vals.data());
        }
    }  // namespace SimdWrapper
}  // namespace Math
}  // namespace Preon