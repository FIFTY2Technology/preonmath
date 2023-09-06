/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"
#include "scalar_wrapper.h"

#ifdef PREONMATH_ENABLE_SIMD
    #include "simd_wrapper.h"

namespace Preon::Math::SimdWrapper
{
/* Base class for transparently storing scalar and SIMD versions of a value. */
template<typename T_Scalar, typename T_Simd>
class Base
{
public:
    using ScalarType = T_Scalar;
    using SimdType = T_Simd;

protected:
    Base(T_Scalar value, T_Simd valueSimd)
        : m_Value(value), m_ValueSimd(valueSimd) {}

    T_Scalar m_Value;
    alignas(16) T_Simd m_ValueSimd;
};

/* This class transparently stores the same value as scalar and SIMD floating point value.
 * It offers the get() method to return the value in the needed type. */
template<typename T_Scalar, typename T_Simd = typename Simd::Register<T_Scalar>::type>
class FloatingPoint : public Base<T_Scalar, T_Simd>
{
public:
    FloatingPoint()
        : FloatingPoint(0) {}
    FloatingPoint(T_Scalar value)
        : Base<T_Scalar, T_Simd>(value, Simd::make<T_Scalar>(value)) {}

    template<typename T>
    T get() const;
};
using Float = FloatingPoint<float>;
using Double = FloatingPoint<double>;

template<>
template<>
inline float Float::get() const
{
    return m_Value;
}
template<>
template<>
inline float_simd Float::get() const
{
    return m_ValueSimd;
}
template<>
template<>
inline double Double::get() const
{
    return m_Value;
}
template<>
template<>
inline double_simd Double::get() const
{
    return m_ValueSimd;
}

/* This class transparently stores the same value as vec3f and as SIMD version.
 * It offers the get() method to return the value in the needed type. */
class Vec3f : public Base<vec<3, float>, vec<3, float_simd>>
{
public:
    Vec3f(vec<3, float> value)
        : Base(value, vec<3, float_simd>(value)) {}

    template<typename T>
    T get() const;
};

template<>
inline vec<3, float> Vec3f::get() const
{
    return m_Value;
}

template<>
inline vec<3, float_simd> Vec3f::get() const
{
    return m_ValueSimd;
}

// Operators
template<
    template<typename, typename>
    class WrapperType,
    typename T,
    typename = std::enable_if_t<std::is_convertible<WrapperType<T, typename Simd::Register<T>::type>, Base<T, typename Simd::Register<T>::type>>::value>>
T operator*(const WrapperType<T, typename Simd::Register<T>::type>& a, const WrapperType<T, typename Simd::Register<T>::type>& b)
{
    return a.template get<T>() * b.template get<T>();
}

template<
    template<typename, typename>
    class WrapperType,
    typename T,
    typename = std::enable_if_t<std::is_convertible<WrapperType<T, typename Simd::Register<T>::type>, Base<T, typename Simd::Register<T>::type>>::value>,
    typename = std::enable_if_t<std::is_arithmetic<T>::value>>
T operator*(const WrapperType<T, typename Simd::Register<T>::type>& wrapped, const T scalar)
{
    return wrapped.template get<T>() * scalar;
}

template<
    template<typename, typename>
    class WrapperType,
    typename T,
    typename = std::enable_if_t<std::is_convertible<WrapperType<T, typename Simd::Register<T>::type>, Base<T, typename Simd::Register<T>::type>>::value>,
    typename = std::enable_if_t<std::is_arithmetic<T>::value>>
T operator*(const T scalar, const WrapperType<T, typename Simd::Register<T>::type>& wrapped)
{
    return scalar * wrapped.template get<T>();
}

template<
    template<typename, typename>
    class WrapperType,
    typename T,
    typename = std::enable_if_t<std::is_convertible<WrapperType<T, typename Simd::Register<T>::type>, Base<T, typename Simd::Register<T>::type>>::value>,
    typename = std::enable_if_t<std::is_arithmetic<T>::value>>
T operator/(const WrapperType<T, typename Simd::Register<T>::type>& wrapped, const T scalar)
{
    return wrapped.template get<T>() / scalar;
}

}  // namespace Preon::Math::SimdWrapper

#endif