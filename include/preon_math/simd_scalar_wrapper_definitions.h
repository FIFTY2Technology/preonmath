/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"
#include "vec_fwd.h"

namespace Preon::Math::SimdWrapper
{
template<typename T>
struct Scalar;

template<typename T>
PREONMATH_FORCEINLINE T make(float val);
template<typename T>
PREONMATH_FORCEINLINE T make(double val);

template<typename T>
PREONMATH_FORCEINLINE T make(const vec3f& val);

template<typename T>
PREONMATH_FORCEINLINE T zero();

}  // namespace Preon::Math::SimdWrapper
