/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_STATIC_FOR_H
#define PREONMATH_STATIC_FOR_H

#include "compile_helper.h"
#include <cstddef>
#include <type_traits>

namespace Preon
{
    namespace Math
    {
        template <size_t First, size_t Last, size_t StepSize = 1, typename Lambda>
        PREONMATH_FORCEINLINE typename std::enable_if<First == Last, void>::type StaticFor(const Lambda&) {}

        // With C++17 we could use fold expressions to avoid recursive calls.
        template <size_t First, size_t Last, size_t StepSize = 1, typename Lambda>
        PREONMATH_FORCEINLINE typename std::enable_if<First != Last, void>::type StaticFor(const Lambda& func)
        {
            func(First);
            StaticFor<First + StepSize, Last, StepSize, Lambda>(func);
        }
    }  // namespace Math
}  // namespace Preon

#endif  // PREONMATH_STATIC_FOR_H
