/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "scalar_simd.h"
#include "vec_simd.h"
#include "vec3.h"
#include "matrix33.h"
#include "matrix_simd.h"

#include <type_traits>

namespace Preon
{
namespace Math
{
    namespace SimdWrapper
    {
        /* Base class for transparently storing scalar and SIMD versions of a value. */
        template<typename T_Scalar, typename T_Simd>
        class Base
        {
        protected:
            Base(T_Scalar value, T_Simd valueSimd)
                : m_Value(value), m_ValueSimd(valueSimd) {}

            T_Scalar m_Value;
            T_Simd m_ValueSimd;
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

        /** Overloaded wrapper methods for transparent handling of scalar and SIMD functions **/

        // scalar make wrapper
        template<typename T>
        T make(float val);
        template<typename T>
        T make(double val);

        template<>
        inline float make<float>(float val)
        {
            return val;
        }
        template<>
        inline float_simd make<float_simd>(float val)
        {
            return Simd::make<float>(val);
        }
        template<>
        inline double make<double>(double val)
        {
            return val;
        }
        template<>
        inline double_simd make<double_simd>(double val)
        {
            return Simd::make<double>(val);
        }
        // also overload the following for convenience
        template<>
        inline float make<float>(double val)
        {
            return val;
        }
        template<>
        inline float_simd make<float_simd>(double val)
        {
            return Simd::make<float>(val);
        }

        // vector make wrapper
        template<typename T>
        T make(const vec3f& val);

        template<>
        inline vec3f make<vec3f>(const vec3f& val)
        {
            return val;
        }

        template<>
        inline vec_simd<3> make<vec_simd<3>>(const vec3f& val)
        {
            return vec_simd<3>(val);
        }

        template<typename T>
        T zero();

        template<>
        inline float zero<float>()
        {
            return 0.0f;
        }
        template<>
        inline double zero<double>()
        {
            return 0.0;
        }
        template<>
        inline float_simd zero<float_simd>()
        {
            return Simd::zero<float>();
        }
        template<>
        inline double_simd zero<double_simd>()
        {
            return Simd::zero<double>();
        }

        // scalar hSum wrapper
        inline float hSum(const float& val) { return val; }
        inline float hSum(const float_simd& val) { return Simd::hSum(val); }

        // vector hSum wrapper
        inline vec3f hSum(const vec3f& val) { return val; }
        inline vec3f hSum(const vec<3, float_simd>& val) { return Simd::hSum(val); }

        // matrix hSum wrapper
        inline matrix33f hSum(const matrix33f& val) { return val; }
        inline matrix33f hSum(const matrix_simd<3, 3, float>& val) { return Simd::hSum(val); }

        // choose
        inline float choose(const bool mask, const float v1, const float v2) { return mask ? v1 : v2; }
        inline float_simd choose(const float_simd mask, const float_simd v1, const float_simd v2) { return Simd::choose(mask, v1, v2); }

        template<size_t D>
        inline vec<D, float> choose(const bool mask, const vec<D, float>& v1, const vec<D, float>& v2)
        {
            return mask ? v1 : v2;
        }
        template<size_t D>
        inline vec<D, float_simd> choose(const float_simd mask, const vec<D, float_simd>& v1, const vec<D, float_simd>& v2)
        {
            return vec_simd<D>::choose(mask, v1, v2);
        }

        // AND mask
        inline float and_mask(bool mask, float val) { return choose(mask, val, 0); }
        inline float and_mask(float val, bool mask) { return choose(mask, val, 0); }
        inline float_simd and_mask(float_simd a, float_simd b) { return Simd::and_mask(a, b); }

        // max
        inline float max(float a, float b) { return std::max(a, b); }
        inline float_simd max(float_simd a, float_simd b) { return Simd::max(a, b); }

        // min
        inline float min(float a, float b) { return std::min(a, b); }
        inline float_simd min(float_simd a, float_simd b) { return Simd::min(a, b); }

        // cosine
        inline float cos(float val) { return std::cos(val); }
        inline float_simd cos(float_simd val)
        {
            // TODO: this is not very efficient
            std::array<float, Simd::Register<float>::size> vals;
            Simd::store(val, vals.data());
            std::for_each(vals.begin(), vals.end(), [](float& v) { v = cos(v); });
            return Simd::load(vals.data());
        }
    }  // namespace SimdWrapper
}  // namespace Math
}  // namespace Preon
