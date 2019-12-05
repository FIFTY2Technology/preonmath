/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_MATRIX44_H
#define PREONMATH_MATRIX44_H

#include "compile_helper.h"

#include "matrix.h"
#include "scalar_simd.h"

namespace Preon
{
    namespace Math
    {
        // Matrix 4x4
        typedef matrix<4, 4, float> matrix44f;
        typedef matrix<4, 4, double> matrix44d;
    }  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
Q_DECLARE_TYPEINFO(matrix44f, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(matrix44d, Q_MOVABLE_TYPE);

Q_DECLARE_METATYPE(matrix44f)
Q_DECLARE_METATYPE(matrix44d)
#endif  // PREONMATH_QT_INTEGRATION

#endif  // PREONMATH_MATRIX44_H
