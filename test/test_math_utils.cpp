/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include "catch_amalgamated.hpp"

#include "preon_math/math_utils.h"

using namespace Preon::Math;

TEST_CASE("MathUtils::testIsZero", "[math_utils]")
{
    // Float.
    CHECK(IsZero(0.0f));
    CHECK(IsZero(-0.0f));
    CHECK(IsZero(static_cast<float>(DEps)));
    CHECK(IsZero(-static_cast<float>(DEps)));
    CHECK(IsZero(std::numeric_limits<float>::min()));
    CHECK(IsZero(-std::numeric_limits<float>::min()));

    CHECK(!IsZero(FEps));
    CHECK(!IsZero(-FEps));
    CHECK(!IsZero(std::numeric_limits<float>::max()));
    CHECK(!IsZero(-std::numeric_limits<float>::max()));

    // Double.
    CHECK(IsZero(0.0));
    CHECK(IsZero(-0.0));
    CHECK(IsZero(std::numeric_limits<double>::min()));
    CHECK(IsZero(-std::numeric_limits<double>::min()));

    CHECK(!IsZero(DEps));
    CHECK(!IsZero(-DEps));
    CHECK(!IsZero(std::numeric_limits<double>::max()));
    CHECK(!IsZero(-std::numeric_limits<double>::max()));

    // Int.
    CHECK(IsZero(0));
    CHECK(IsZero(-0));

    CHECK(!IsZero(1));
    CHECK(!IsZero(-1));
}

TEST_CASE("MathUtils::testIsExactlyZero", "[math_utils]")
{
    // Float.
    CHECK(IsExactlyZero(0.0f));
    CHECK(IsExactlyZero(-0.0f));

    CHECK(!IsExactlyZero(std::numeric_limits<float>::min()));
    CHECK(!IsExactlyZero(-std::numeric_limits<float>::min()));
    CHECK(!IsExactlyZero(FEps));
    CHECK(!IsExactlyZero(-FEps));

    // Double.
    CHECK(IsExactlyZero(0.0));
    CHECK(IsExactlyZero(-0.0));

    CHECK(!IsExactlyZero(std::numeric_limits<double>::min()));
    CHECK(!IsExactlyZero(-std::numeric_limits<double>::min()));
    CHECK(!IsExactlyZero(DEps));
    CHECK(!IsExactlyZero(-DEps));

    // Int.
    CHECK(IsExactlyZero(0));
    CHECK(IsExactlyZero(-0));

    CHECK(!IsExactlyZero(1));
    CHECK(!IsExactlyZero(-1));
}

TEST_CASE("MathUtils::testAreEqual", "[math_utils]")
{
    // Float.
    CHECK(AreEqual(0.0f, 0.0f));
    CHECK(AreEqual(0.0f, -0.0f));
    CHECK(AreEqual(3.14159f, 3.14159f));
    CHECK(AreEqual(-15418.45f, -15418.45f + std::numeric_limits<float>::min()));
    CHECK(AreEqual(std::numeric_limits<float>::min(), std::numeric_limits<float>::min()));
    CHECK(AreEqual(0.0f, std::numeric_limits<float>::min()));

    CHECK(!AreEqual(12310.0f, 1256460.0f));
    CHECK(!AreEqual(-10.0f, 10.0f));

    // Double.
    CHECK(AreEqual(0.0, 0.0));
    CHECK(AreEqual(0.0, -0.0));
    CHECK(AreEqual(3.14159, 3.14159));
    CHECK(AreEqual(-15418.45, -15418.45 + std::numeric_limits<double>::min()));
    CHECK(AreEqual(std::numeric_limits<double>::min(), std::numeric_limits<double>::min()));
    CHECK(AreEqual(0.0, std::numeric_limits<double>::min()));

    CHECK(!AreEqual(12310.0, 1256460.0));
    CHECK(!AreEqual(-10.0, 10.0));

    // Int.
    CHECK(AreEqual(1, 1));
    CHECK(AreEqual(0, -0));
    CHECK(!AreEqual(1, -1));
}

TEST_CASE("MathUtils::testAreExactlyEqual", "[math_utils]")
{
    // Float.
    CHECK(AreExactlyEqual(0.0f, 0.0f));
    CHECK(AreExactlyEqual(0.0f, -0.0f));
    CHECK(AreExactlyEqual(1231.787673611123357857f, 1231.787673611123357857f));

    CHECK(!AreExactlyEqual(0.0f, std::numeric_limits<float>::min()));

    // Double,
    CHECK(AreExactlyEqual(0.0, 0.0));
    CHECK(AreExactlyEqual(0.0, -0.0));
    CHECK(AreExactlyEqual(1231.787673611123357857, 1231.787673611123357857));

    CHECK(!AreExactlyEqual(0.0, std::numeric_limits<double>::min()));

    // Int.
    CHECK(AreExactlyEqual(1, 1));
    CHECK(AreExactlyEqual(-0, 0));

    CHECK(!AreExactlyEqual(11276567, 0));
    CHECK(!AreExactlyEqual(11276567, 11276568));
}

TEST_CASE("MathUtils::testSignum", "[math_utils]")
{
    // Float.
    CHECK(Signum(0.0f) == 0);
    CHECK(Signum(-0.0f) == 0);
    CHECK(Signum(112435.125537568f) == 1);
    CHECK(Signum(-112435.125537568f) == -1);
    CHECK(Signum(-std::numeric_limits<float>::min()) == -1);
    CHECK(Signum(std::numeric_limits<float>::min()) == 1);
    CHECK(Signum(-std::numeric_limits<float>::max()) == -1);
    CHECK(Signum(std::numeric_limits<float>::max()) == 1);

    // Double.
    CHECK(Signum(0.0) == 0);
    CHECK(Signum(-0.0) == 0);
    CHECK(Signum(112435.125537568) == 1);
    CHECK(Signum(-112435.125537568) == -1);
    CHECK(Signum(-std::numeric_limits<double>::min()) == -1);
    CHECK(Signum(std::numeric_limits<double>::min()) == 1);
    CHECK(Signum(-std::numeric_limits<double>::max()) == -1);
    CHECK(Signum(std::numeric_limits<double>::max()) == 1);

    // Int.
    CHECK(Signum(0) == 0);
    CHECK(Signum(-0) == 0);
    CHECK(Signum(-789780) == -1);
    CHECK(Signum(-std::numeric_limits<int>::max()) == -1);
    CHECK(Signum(std::numeric_limits<int>::max()) == 1);
}

TEST_CASE("MathUtils::testDegToRad", "[math_utils]")
{
    // Float.
    CHECK(DegToRad(0.0f) == 0.0f);
    CHECK(DegToRad(360.0f) == 2.0f * FPi);
    CHECK(DegToRad(90.0f) == 0.5f * FPi);
    CHECK(DegToRad(-180.0f) == -FPi);

    // Double.
    CHECK(DegToRad(0.0) == 0.0);
    CHECK(DegToRad(360.0) == 2.0 * DPi);
    CHECK(DegToRad(90.0) == 0.5 * DPi);
    CHECK(DegToRad(-180.0) == -DPi);
}

TEST_CASE("MathUtils::testRadToDeg", "[math_utils]")
{
    // Float.
    CHECK(0.0f == DegToRad(0.0f));
    CHECK(2.0f * FPi == DegToRad(360.0f));
    CHECK(0.5f * FPi == DegToRad(90.0f));
    CHECK(-FPi == DegToRad(-180.0f));

    // Double.
    CHECK(0.0 == DegToRad(0.0));
    CHECK(2.0 * DPi == DegToRad(360.0));
    CHECK(0.5 * DPi == DegToRad(90.0));
    CHECK(-DPi == DegToRad(-180.0));
}

TEST_CASE("MathUtils::nextPowerOfTwo", "[math_utils]")
{
    CHECK(MathUtils::nextPowerOfTwo(0) == 1);
    CHECK(MathUtils::nextPowerOfTwo(1) == 1);
    CHECK(MathUtils::nextPowerOfTwo(2) == 2);
    CHECK(MathUtils::nextPowerOfTwo(3) == 4);
    CHECK(MathUtils::nextPowerOfTwo(64) == 64);
    CHECK(MathUtils::nextPowerOfTwo(12364) == 16384);

    REQUIRE_THROWS_AS(MathUtils::nextPowerOfTwo(-1), std::exception);
    REQUIRE_THROWS_AS(MathUtils::nextPowerOfTwo(-510), std::exception);
}
