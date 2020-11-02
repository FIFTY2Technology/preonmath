/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "vec.h"
#include "matrix.h"
#include "math_utils.h"

namespace Preon
{
namespace Math
{
    namespace VecUtils
    {
        template<size_t D, typename T>
        static vec<D, T> round(vec<D, T> v)
        {
            for (size_t i = 0; i < D; i++)
                v[i] = MathUtils::round(v[i]);
            return v;
        }

        template<size_t D, typename T>
        static vec<D, T> floor(vec<D, T> v)
        {
            for (size_t i = 0; i < D; i++)
                v[i] = std::floor(v[i]);
            return v;
        }

        template<size_t D, typename T>
        static vec<D, T> ceil(vec<D, T> v)
        {
            for (size_t i = 0; i < D; i++)
                v[i] = std::ceil(v[i]);
            return v;
        }

        template<size_t D, typename T>
        static bool isNan(const vec<D, T>& v)
        {
            for (size_t i = 0; i < D; i++)
                if (std::isnan(v[i]))
                    return true;
            return false;
        }

        template<size_t D, typename T>
        static matrix<D, D, T> outerProduct(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return matrix<D, D, T>([&](size_t r, size_t c) { return v1[r] * v2[c]; });
        }

        template<size_t D, typename T>
        static vec<D, T> project(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return (v1 * v2) / v2.lengthSquared() * v2;
        }

        //! Projects the vector \arg \v into the direction pointed by \arg \n. \arg \n has to be normalized.
        template<size_t D, typename T>
        static vec<D, T> projectToDirection(const vec<D, T>& v, const vec<D, T>& n)
        {
            return (v * n) * n;
        }

        template<size_t D, typename T>
        static vec<D, T> reject(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return v1 - project(v1, v2);
        }

        //! Rejects the vector \arg \v from the direction pointed by \arg \n. \arg \n has to be normalized.
        template<size_t D, typename T>
        static vec<D, T> rejectFromDirection(const vec<D, T>& v, const vec<D, T>& n)
        {
            return v - projectToDirection(v, n);
        }
    }  // namespace VecUtils
}  // namespace Math
}  // namespace Preon
