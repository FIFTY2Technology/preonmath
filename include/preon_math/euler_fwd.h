/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

namespace Preon
{
namespace Math
{
    template<typename Real>
    class euler;
    typedef euler<float> eulerf;
    typedef euler<double> eulerd;
}  // namespace Math
}  // namespace Preon
