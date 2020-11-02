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

    typedef SymmetricMatrix33<float> SymmetricMatrix33f;
    typedef SymmetricMatrix33<double> SymmetricMatrix33d;
}  // namespace Math
}  // namespace Preon
