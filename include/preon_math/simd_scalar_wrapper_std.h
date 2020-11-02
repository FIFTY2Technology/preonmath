/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#include "compile_helper.h"

namespace Preon
{
namespace Math
{
    namespace SimdWrapper
    {
        // sqrt
        using Simd::sqrt;
        using std::sqrt;

        // add other std functions that are implemented in Simd...
    }  // namespace SimdWrapper
}  // namespace Math
}  // namespace Preon
