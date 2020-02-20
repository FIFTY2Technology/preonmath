/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_VEC4_H
#define PREONMATH_VEC4_H

#include "compile_helper.h"

#include "vec.h"

namespace Preon
{
    namespace Math
    {
        // Vec4.
        typedef vec<4, float> vec4f;
        typedef vec<4, double> vec4d;
        typedef vec<4, int> vec4i;
    }  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
Q_DECLARE_TYPEINFO(vec4f, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(vec4d, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(vec4i, Q_MOVABLE_TYPE);
Q_DECLARE_METATYPE(vec4f)
Q_DECLARE_METATYPE(vec4d)
#endif  // PREONMATH_QT_INTEGRATION

#endif  // PREONMATH_VEC4_H
