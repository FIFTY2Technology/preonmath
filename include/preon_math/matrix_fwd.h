/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"
#include "scalar_simd.h"

#include <cstddef>  // For PrMathSize

namespace Preon
{
namespace Math
{
    template<PrMathSize M, PrMathSize N, typename T>
    class matrix;
    // Matrix 3x3
    typedef matrix<3, 3, float> matrix33f;
    typedef matrix<3, 3, double> matrix33d;
    typedef matrix<3, 3, float_simd> matrix33_Simd;

    // Matrix 4x4
    typedef matrix<4, 4, float> matrix44f;
    typedef matrix<4, 4, double> matrix44d;
}  // namespace Math
}  // namespace Preon
