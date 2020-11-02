/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "vec.h"

namespace Preon
{
namespace Math
{
    // Vec2.
    typedef vec<2, float> vec2f;
    typedef vec<2, double> vec2d;
    typedef vec<2, int> vec2i;
}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
Q_DECLARE_TYPEINFO(vec2f, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(vec2d, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(vec2i, Q_MOVABLE_TYPE);
Q_DECLARE_METATYPE(vec2f)
Q_DECLARE_METATYPE(vec2d)
#endif  // PREONMATH_QT_INTEGRATION
