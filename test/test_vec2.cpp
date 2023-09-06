/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include "catch_amalgamated.hpp"

#include "preon_math/vec2.h"

#include <cmath>

using namespace Preon::Math;

TEST_CASE("vec2f::length", "[vec]")
{
    {
        vec2f test(0.0f, 0.0f);
        CHECK(test.length() == 0.0f);
    }
    {
        vec2f test(1.0f, -1.0f);
        CHECK(test.length() == std::sqrt(2.0f));
    }
    {
        vec2f test(-1.0f, -1.0f);
        CHECK(test.length() == std::sqrt(2.0f));
    }
}

TEST_CASE("vec2f::lengthSquared", "[vec]")
{
    {
        vec2f test(0.0f, 0.0f);
        CHECK(test.length() == 0.0f);
    }
    {
        vec2f test(1.0f, -1.0f);
        CHECK(test.length() == std::sqrt(2.0f));
    }
    {
        vec2f test(-1.0f, -1.0f);
        CHECK(test.length() == std::sqrt(2.0f));
    }
}
