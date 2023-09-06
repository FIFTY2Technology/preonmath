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
        template<PrMathSize D, typename T>
        static vec<D, T> round(vec<D, T> v)
        {
            for (PrMathSize i = 0; i < D; i++)
                v[i] = MathUtils::round(v[i]);
            return v;
        }

        template<PrMathSize D, typename T>
        static vec<D, T> floor(vec<D, T> v)
        {
            for (PrMathSize i = 0; i < D; i++)
                v[i] = std::floor(v[i]);
            return v;
        }

        template<PrMathSize D, typename T>
        static vec<D, T> ceil(vec<D, T> v)
        {
            for (PrMathSize i = 0; i < D; i++)
                v[i] = std::ceil(v[i]);
            return v;
        }

        template<PrMathSize D, typename T>
        static bool isNan(const vec<D, T>& v)
        {
            for (PrMathSize i = 0; i < D; i++)
                if (std::isnan(v[i]))
                    return true;
            return false;
        }

        template<PrMathSize D, typename T1, typename T2>
        PREONMATH_DEVICE static matrix<D, D, product_t<T1, T2>> outerProduct(const vec<D, T1>& v1, const vec<D, T2>& v2)
        {
            return matrix<D, D, product_t<T1, T2>>([&](PrMathSize r, PrMathSize c) { return v1[r] * v2[c]; });
        }

        template<PrMathSize D, typename T>
        PREONMATH_DEVICE static vec<D, T> project(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return (v1 * v2) / v2.lengthSquared() * v2;
        }

        //! Projects the vector \arg \v into the direction pointed by \arg \n. \arg \n has to be normalized.
        template<PrMathSize D, typename Tv, typename Tn>
        PREONMATH_DEVICE static vec<D, Tv> projectToDirection(const vec<D, Tv>& v, const vec<D, Tn>& n)
        {
            return (v * n) * n;
        }

        template<PrMathSize D, typename T>
        PREONMATH_DEVICE static vec<D, T> reject(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return v1 - project(v1, v2);
        }

        //! Rejects the vector \arg \v from the direction pointed by \arg \n. \arg \n has to be normalized.
        template<PrMathSize D, typename Tv, typename Tn>
        PREONMATH_DEVICE static vec<D, Tv> rejectFromDirection(const vec<D, Tv>& v, const vec<D, Tn>& n)
        {
            return v - projectToDirection(v, n);
        }
    }  // namespace VecUtils
}  // namespace Math
}  // namespace Preon
