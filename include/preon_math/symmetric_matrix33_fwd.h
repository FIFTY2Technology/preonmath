/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

namespace Preon
{
namespace Math
{
    template<typename T>
    class SymmetricMatrix33;

    using SymmetricMatrix33f = SymmetricMatrix33<float>;
    using SymmetricMatrix33d = SymmetricMatrix33<double>;
}  // namespace Math
}  // namespace Preon
