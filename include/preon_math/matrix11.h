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
    // Matrix 1x1
    typedef matrix<1, 1, float> matrix11f;
    typedef matrix<1, 1, double> matrix11d;
#ifdef PREONMATH_ENABLE_SIMD
    typedef matrix<1, 1, float_simd> matrix11_Simd;
#endif
}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
Q_DECLARE_TYPEINFO(matrix11f, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(matrix11d, Q_MOVABLE_TYPE);

Q_DECLARE_METATYPE(matrix11f)
Q_DECLARE_METATYPE(matrix11d)
#endif  // PREONMATH_QT_INTEGRATION
