/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

namespace Preon
{
namespace Math
{
    template<PrMathSize D, typename T>
    class vec;

    // Vec2.
    template<typename T>
    using vec2 = vec<2, T>;
    using vec2f = vec2<float>;
    using vec2d = vec2<double>;
    using vec2i = vec2<int>;

    // Vec3.
    template<typename T>
    using vec3 = vec<3, T>;
    using vec3f = vec3<float>;
    using vec3d = vec3<double>;
    using vec3i = vec3<int>;
    using vec3us = vec3<unsigned short>;
    using vec3uc = vec3<unsigned char>;
    using vec3ui = vec3<unsigned int>;

    // Vec4.
    template<typename T>
    using vec4 = vec<4, T>;
    using vec4f = vec4<float>;
    using vec4d = vec4<double>;
    using vec4i = vec4<int>;
}  // namespace Math
}  // namespace Preon
