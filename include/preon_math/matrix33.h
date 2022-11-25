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
    // Matrix 3x3
    typedef matrix<3, 3, float> matrix33f;
    typedef matrix<3, 3, double> matrix33d;
#ifdef PREONMATH_ENABLE_SIMD
    typedef matrix<3, 3, float_simd> matrix33_Simd;
#endif
}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
Q_DECLARE_TYPEINFO(matrix33f, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(matrix33d, Q_MOVABLE_TYPE);

Q_DECLARE_METATYPE(matrix33f)
Q_DECLARE_METATYPE(matrix33d)
#endif  // PREONMATH_QT_INTEGRATION
