/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"
#include "scalar_wrapper_std.h"

#include <cmath>
#include <cstring>

// RadToDeg: 180 / Pi
#define FRadiansToDegree 57.295779513082320876798154814105f
#define DRadiansToDegree 57.295779513082320876798154814105
// DegToRad: Pi / 180
#define FDegreeToRadians 0.01745329251994329576923690768489f
#define DDegreeToRadians 0.01745329251994329576923690768489

// Math constants.
#define FPi 3.1415926535897932384626433832795f
#define FPiInv 0.31830988618f
#define FPiHalf 1.5707963267948966192313216916398f
#define FPiOver180 0.01745329251994329576923690768489f
#define FEps 1.0e-6f
#define FTwoPi 2.0f * FPi
#define DEps 2.2204460492503131e-015
#define D6Eps 2.2204460492503131e-06
#define DPi 3.1415926535897932384626433832795
#define DPiHalf 1.5707963267948966192313216916398
#define DTwoPi 2.0 * DPi

#define F_INFINITY 9999999.0f

namespace Preon
{
//! Copies the content of input into the returned value byte by byte.
//! Both types must have the same size.
template<class TOut, class TIn>
TOut interpretAs(const TIn& input)
{
    static_assert(sizeof(TOut) == sizeof(TIn), "TOut and TIn must have same byte size!");
    TOut out;
    std::memcpy((char*)&out, (char*)&input, sizeof(TIn));
    return out;  // Could also implement this with *reinterpret_cast<TOut*>(&input) but apparently this wouldn't be 100% because pointer alignment requirements may be different for different types.
}

namespace Math
{
    template<typename T>
    struct PreonReal
    {
    };
    template<>
    struct PreonReal<float>
    {
        PREONMATH_FORCEINLINE static float eps() { return FEps; }
    };
    template<>
    struct PreonReal<double>
    {
        PREONMATH_FORCEINLINE static double eps() { return DEps; }
    };

    // IsZero.
    // We do not define it generically,
    // because we want an error if it is called for a non-instantiated class.
    template<typename T>
    PREONMATH_DEVICE inline bool IsZero(T val);
    template<>
    PREONMATH_DEVICE inline bool IsZero<float>(float val)
    {
        return abs(val) < FEps;
    }
    template<>
    constexpr inline bool IsZero<double>(double val)
    {
        return abs(val) < DEps;
    }
    template<>
    PREONMATH_DEVICE inline bool IsZero<int>(int val)
    {
        return val == 0;
    }

    // IsExactlyZero.
    template<typename T>
    PREONMATH_DEVICE inline bool IsExactlyZero(T val);
    template<>
    PREONMATH_DEVICE inline bool IsExactlyZero<float>(float val)
    {
        return abs(val) == 0.0f;
    }
    template<>
    PREONMATH_DEVICE inline bool IsExactlyZero<double>(double val)
    {
        return abs(val) == 0.0;
    }
    template<>
    PREONMATH_DEVICE inline bool IsExactlyZero<int>(int val)
    {
        return val == 0;
    }

    // AreEqual.
    template<typename T>
    PREONMATH_DEVICE inline bool AreEqual(T a, T b)
    {
        return IsZero<T>(a - b);
    }

    // AreExactlyEqual.
    template<typename T>
    PREONMATH_DEVICE inline bool AreExactlyEqual(T a, T b)
    {
        return IsExactlyZero<T>(a - b);
    }

    // Signum function.
    template<typename T>
    PREONMATH_DEVICE inline int Signum(T val);
    template<>
    PREONMATH_DEVICE inline int Signum<float>(float val)
    {
        return (val > 0.0f) - (val < 0.0f);
    }
    template<>
    PREONMATH_DEVICE inline int Signum<double>(double val)
    {
        return (val > 0.0) - (val < 0.0);
    }
    template<>
    PREONMATH_DEVICE inline int Signum<int>(int val)
    {
        return (val > 0) - (val < 0);
    }

    // DegToRad.
    template<typename T>
    PREONMATH_DEVICE inline T DegToRad(T val);
    template<>
    PREONMATH_DEVICE inline float DegToRad<float>(float degree)
    {
        return FDegreeToRadians * degree;
    }
    template<>
    PREONMATH_DEVICE inline double DegToRad<double>(double degree)
    {
        return DDegreeToRadians * degree;
    }

    // RadToDeg.
    template<typename T>
    PREONMATH_DEVICE inline T RadToDeg(T val);
    template<>
    PREONMATH_DEVICE inline float RadToDeg<float>(float radians)
    {
        return FRadiansToDegree * radians;
    }
    template<>
    PREONMATH_DEVICE inline double RadToDeg<double>(double radians)
    {
        return DRadiansToDegree * radians;
    }

    // Returns the value of the given argument as primitive type. Can be specialized for custom data types.
    template<typename T>
    PREONMATH_DEVICE T PrimitiveType(const T& val)
    {
        return val;
    }

    // Resulting type for multiplication.
    template<typename T1, typename T2>
    using product_t = decltype(std::declval<T1>() * std::declval<T2>());

    // Resulting type for division.
    template<typename T1, typename T2>
    using quotient_t = decltype(std::declval<T1>() / std::declval<T2>());

    // Depracted, just here for existing code, please do not use anymore!
    inline bool FIsZero(float a) { return IsZero<float>(a); }
    inline bool FIsExactlyZero(float a) { return IsExactlyZero<float>(a); }
    inline bool FAreEqual(float a, float b) { return AreEqual<float>(a, b); }
    inline bool FAreExactlyEqual(float a, float b) { return AreExactlyEqual<float>(a, b); }
    inline int FSignum(float a) { return Signum<float>(a); }
    inline float FDegToRad(float deg) { return DegToRad<float>(deg); }
    inline float FRadToDeg(float radians) { return RadToDeg<float>(radians); }

    inline bool DIsZero(double a) { return IsZero<double>(a); }
    inline bool DIsExactlyZero(double a) { return IsExactlyZero<double>(a); }
    inline bool DAreEqual(double a, double b) { return AreEqual<double>(a, b); }
    inline bool DAreExactlyEqual(double a, double b) { return AreExactlyEqual<double>(a, b); }
    inline double DDegToRad(double deg) { return DegToRad<double>(deg); }
    inline double DRadToDeg(double radians) { return RadToDeg<double>(radians); }
}  // namespace Math
}  // namespace Preon
