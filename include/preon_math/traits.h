/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "vec_fwd.h"
#include "euler.h"

#include <type_traits>

namespace Preon
{
namespace Math
{
    /* Dimension of scalar and vectors */
    struct Dimension
    {
        template<typename T>
        static PrMathSize get();
    };

    template<>
    inline PrMathSize Dimension::get<float>()
    {
        return 1;
    }

    template<>
    inline PrMathSize Dimension::get<vec3f>()
    {
        return 3;
    }

    template<typename T>
    struct Wrapper
    {
        template<typename Lambda>
        static void forAll(const T& val, const Lambda& lambda);
    };

    template<>
    template<typename Lambda>
    inline void Wrapper<float>::forAll(const float& val, const Lambda& lambda)
    {
        lambda(val);
    }

    template<>
    template<typename Lambda>
    inline void Wrapper<vec3f>::forAll(const vec3f& val, const Lambda& lambda)
    {
        val.forAll([&](PrMathSize i) { lambda(val[i]); });
    }

    template<typename T>
    struct IsEuler : std::false_type
    {
    };
    template<typename T>
    struct IsEuler<euler<T>> : std::true_type
    {
    };
    template<typename T>
    inline constexpr bool IsEulerV = IsEuler<T>::value;

}  // namespace Math
}  // namespace Preon
