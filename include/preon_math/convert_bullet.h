/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_CONVERT_BULLET_H
#define PREONMATH_CONVERT_BULLET_H

#include "compile_helper.h"

#include "vec3.h"
#include "matrix33.h"
#include "quat.h"

#include <bullet/LinearMath/btQuaternion.h>
#include <bullet/LinearMath/btVector3.h>
#include <bullet/LinearMath/btMatrix3x3.h>

namespace Preon
{
    namespace Math
    {
        inline vec3f fromBullet(const btVector3& vec)
        {
            return vec3f(vec.getX(), vec.getY(), vec.getZ());
        }

        template <typename T>
        inline btVector3 toBullet(const vec<3, T>& vec)
        {
            return btVector3(vec.x(), vec.y(), vec.z());
        }

        inline quatf fromBullet(const btQuaternion& quat)
        {
            return quatf(quat.getW(), quat.getX(), quat.getY(), quat.getZ());
        }

        inline btQuaternion toBullet(const quatf& quat)
        {
            return btQuaternion(quat.x(), quat.y(), quat.z(), quat.scalar());
        }

        inline btMatrix3x3 toBullet(const matrix33f& mat)
        {
            return btMatrix3x3(mat(0, 0),
                               mat(0, 1),
                               mat(0, 2),
                               mat(1, 0),
                               mat(1, 1),
                               mat(1, 2),
                               mat(2, 0),
                               mat(2, 1),
                               mat(2, 2));
        }
    }  // namespace Math
}  // namespace Preon

#endif  // PREONMATH_CONVERT_BULLET_H
