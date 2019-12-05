/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_TRAITS_H
#define PREONMATH_TRAITS_H

#include "vec_fwd.h"

namespace Preon
{
    namespace Math
    {
        /* Dimension of scalar and vectors */
        struct Dimension
        {
            template <typename T>
            static size_t get();
        };

        template <>
        inline size_t Dimension::get<float>() { return 1; }

        template <>
        inline size_t Dimension::get<vec3f>() { return 3; }


        template <typename T>
        struct Wrapper
        {
            template <typename Lambda>
            static void forAll(const T& val, const Lambda& lambda);
        };

        template <> template <typename Lambda>
        inline void Wrapper<float>::forAll(const float& val, const Lambda& lambda) { lambda(val); }

        template <> template <typename Lambda>
        inline void Wrapper<vec3f>::forAll(const vec3f& val, const Lambda& lambda) { val.forAll([&](size_t i) { lambda(val[i]); }); }
    }  // namespace Math
}  // namespace Preon

#endif  // PREONMATH_TRAITS_H
