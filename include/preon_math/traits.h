/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"
#include "vec_fwd.h"
#include "euler.h"

#ifdef PREONMATH_UNITS_INTEGRATION
    #include <units/core.h>
#endif

#include <tuple>
#include <type_traits>

namespace Preon
{
namespace Math
{
    template<typename T>
    inline constexpr bool isSupportedScalar = std::is_arithmetic_v<T>
#ifdef PREONMATH_UNITS_INTEGRATION
        || units::traits::is_unit_v<T>
#endif
        ;

    /* Dimension of scalar and vectors */
    template<typename T>
    struct Dimension
    {
        static_assert(isSupportedScalar<T>, "No valid specialization found for the provided data type");

        static constexpr PrMathSize get() { return 1; }
    };

    template<PrMathSize D, typename T>
    struct Dimension<vec<D, T>>
    {
        static constexpr PrMathSize get() { return D; }
    };

    template<typename T>
    struct DataSize
    {
        static constexpr PrMathSize get() { return std::tuple_size<std::remove_reference_t<decltype(std::declval<T>().data())>>::value; }
    };

// Unfortunatelly, our CUDA version does not support C++20 yet. This can be included once we upgrade to CUDA 12.
#if __cplusplus >= 202002L
    template<typename T>
    requires(isSupportedScalar<T>) struct DataSize<T>
    {
        static constexpr PrMathSize get() { return 1; }
    };
#endif

    template<typename T>
    struct Wrapper
    {
        static_assert(isSupportedScalar<T>, "No valid specialization found for the provided data type");

        template<typename Lambda>
        constexpr static void forAll(const T& val, const Lambda& lambda)
        {
            lambda(val);
        }
    };

    template<PrMathSize D, typename T>
    struct Wrapper<vec<D, T>>
    {
        template<typename Lambda>
        constexpr static void forAll(const vec<D, T>& val, const Lambda& lambda)
        {
            val.forAll([&](PrMathSize i) { lambda(val[i]); });
        }
    };

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
