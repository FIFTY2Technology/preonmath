/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_VEC_UTILS_H
#define PREONMATH_VEC_UTILS_H

#include "vec.h"
#include "matrix.h"
#include "math_utils.h"

namespace Preon
{
    namespace Math
    {
        namespace VecUtils
        {
            template <size_t D, typename T>
            static vec<D, T> round(vec<D, T> v)
            {
                for (size_t i = 0; i < D; i++)
                    v[i] = MathUtils::round(v[i]);
                return v;
            }

            template <size_t D, typename T>
            static vec<D, T> floor(vec<D, T> v)
            {
                for (size_t i = 0; i < D; i++)
                    v[i] = std::floor(v[i]);
                return v;
            }

            template <size_t D, typename T>
            static vec<D, T> ceil(vec<D, T> v)
            {
                for (size_t i = 0; i < D; i++)
                    v[i] = std::ceil(v[i]);
                return v;
            }

            template <size_t D, typename T>
            static bool isNan(const vec<D, T>& v)
            {
                for (size_t i = 0; i < D; i++)
                    if (std::isnan(v[i])) return true;
                return false;
            }

            template <size_t D, typename T>
            static matrix<D, D, T> outerProduct(const vec<D, T>& v1, const vec<D, T>& v2)
            {
                return matrix<D, D, T>([&] (size_t r, size_t c) { return v1[r] * v2[c]; } );
            }
        };
    }  // namespace Math
}  // namespace Preon

#endif  // PREONMATH_VEC_UTILS_H
