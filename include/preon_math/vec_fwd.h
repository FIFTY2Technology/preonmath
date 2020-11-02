/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include <cstddef>  // For size_t

namespace Preon
{
namespace Math
{
    template<std::size_t D, typename T>
    class vec;

    // Vec2.
    typedef vec<2, float> vec2f;
    typedef vec<2, double> vec2d;
    typedef vec<2, int> vec2i;

    // Vec3.
    typedef vec<3, float> vec3f;
    typedef vec<3, double> vec3d;
    typedef vec<3, int> vec3i;

    // Vec4.
    typedef vec<4, float> vec4f;
    typedef vec<4, double> vec4d;
    typedef vec<4, int> vec4i;
}  // namespace Math
}  // namespace Preon
