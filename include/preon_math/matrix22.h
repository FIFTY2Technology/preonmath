/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "matrix.h"
#include "scalar_simd.h"

namespace Preon
{
namespace Math
{
    // Matrix 2x2
    typedef matrix<2, 2, float> matrix22f;
    typedef matrix<2, 2, double> matrix22d;
    typedef matrix<2, 2, float_simd> matrix22_Simd;
}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
Q_DECLARE_TYPEINFO(matrix22f, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(matrix22d, Q_MOVABLE_TYPE);

Q_DECLARE_METATYPE(matrix22f)
Q_DECLARE_METATYPE(matrix22d)
#endif  // PREONMATH_QT_INTEGRATION
