/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "vec.h"
#include "matrix.h"

#include "ImathVec.h"
#include "ImathMatrix.h"

namespace Preon
{
namespace Math
{
    // Vec2.
    template<typename T>
    inline const vec<2, T>& fromImath(const Imath::Vec2<T>& v)
    {
        return *reinterpret_cast<const vec<2, T>*>(&v);
    }
    template<typename T>
    inline const Imath::Vec2<T>& toImath(const vec<2, T>& v)
    {
        return *reinterpret_cast<const Imath::Vec2<T>*>(&v);
    }

    // Vec3.
    template<typename T>
    inline const vec<3, T>& fromImath(const Imath::Vec3<T>& v)
    {
        return *reinterpret_cast<const vec<3, T>*>(&v);
    }
    template<typename T>
    inline const Imath::Vec3<T>& toImath(const vec<3, T>& v)
    {
        return *reinterpret_cast<const Imath::Vec3<T>*>(&v);
    }

    // Matrix44.
    template<typename TIn, typename TOut = TIn>
    inline Imath::Matrix44<TOut> toImath(const matrix<4, 4, TIn>& m)
    {
        return Imath::Matrix44<TOut>(
            (TOut)m(0, 0),
            (TOut)m(1, 0),
            (TOut)m(2, 0),
            (TOut)m(3, 0),
            (TOut)m(0, 1),
            (TOut)m(1, 1),
            (TOut)m(2, 1),
            (TOut)m(3, 1),
            (TOut)m(0, 2),
            (TOut)m(1, 2),
            (TOut)m(2, 2),
            (TOut)m(3, 2),
            (TOut)m(0, 3),
            (TOut)m(1, 3),
            (TOut)m(2, 3),
            (TOut)m(3, 3));
    }
}  // namespace Math
}  // namespace Preon
