/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "basics.h"
#include "initializers.h"
#include "scalar_simd.h"
#include "scalar_wrapper_std.h"
#include "static_for.h"

#ifdef PREONMATH_QT_INTEGRATION
    #include <QDataStream>
    #include <QDebug>
    #include <QString>
#endif  // PREONMATH_QT_INTEGRATION

#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <type_traits>
#include <sstream>
#include <iostream>

namespace Preon
{
namespace Math
{
    namespace
    {
        template<typename T, typename U>
        using enable_if_not_convertible = typename std::enable_if<!std::is_convertible<T, U>::value>::type*;

        template<typename T>
        using enable_if_scalar = std::enable_if_t<std::is_arithmetic<T>::value || is_simd_scalar<T>::value>;
    }  // namespace

// Replace with std::isgreater when we use c++14.
#define ISGREATER(a, b) (a > b)

    template<PrMathSize D, typename T>
    class vec
    {
    private:
        //! legacyInitZeroFunc initializes the components of the vector to zero for non-SIMD types.
        //! This is necessary for legacy reasons, but in the future, we want to move away from it.
        template<typename E>
        PREONMATH_DEVICE void legacyInitZeroFunc(enable_if_non_simd<E, void*> = nullptr)
        {
#ifdef PREONMATH_DEFAULTINIT
            setZero();
#endif
        }

        //! Do nothing for SIMD types here.
        template<typename E>
        PREONMATH_DEVICE void legacyInitZeroFunc(enable_if_simd<E, void*> = nullptr)
        {
        }

        template<PrMathSize Begin, PrMathSize End>
        struct AssignLoopUnroller
        {
            template<typename Lambda>
            PREONMATH_FORCEINLINE static void step(vec<D, T>* v, const Lambda& func)
            {
                (*v)[Begin] = func(Begin);
                AssignLoopUnroller<Begin + 1, End>::step(v, func);
            }
        };
        template<PrMathSize End>
        struct AssignLoopUnroller<End, End>
        {
            template<typename Lambda>
            PREONMATH_FORCEINLINE static void step(vec<D, T>*, const Lambda&)
            {
            }
        };

    public:
        /** Constructors **/
        PREONMATH_FORCEINLINE vec()
        {
            legacyInitZeroFunc<T>();
        }
        PREONMATH_FORCEINLINE explicit vec(T fillValue)
        {
            setAll([&](PrMathSize) { return fillValue; });
        }

        template<PrMathSize E = D>
        PREONMATH_FORCEINLINE vec(T x, T y, typename std::enable_if<E == 2, T>::type* = 0)
        {
            set(x, y);
        }
        template<PrMathSize E = D>
        PREONMATH_FORCEINLINE vec(T x, T y, T z, typename std::enable_if<E == 3, T>::type* = 0)
        {
            set(x, y, z);
        }
        template<PrMathSize E = D>
        PREONMATH_FORCEINLINE vec(T x, T y, T z, T w, typename std::enable_if<E == 4, T>::type* = 0)
        {
            set(x, y, z, w);
        }
        // MSVC currently has problems with this variadic args cnstructor, so we only implement the ones we need currently.
#ifdef PREONMATH_COMPILER_MSVC
        template<PrMathSize E = D>
        PREONMATH_FORCEINLINE vec(T a, T b, T c, T d, T e, typename std::enable_if<E == 5, T>::type* = 0)
        {
            set(a, b, c, d, e);
        }
        template<PrMathSize E = D>
        PREONMATH_FORCEINLINE vec(T a, T b, T c, T d, T e, T f, typename std::enable_if<E == 6, T>::type* = 0)
        {
            set(a, b, c, d, e, f);
        }
#else
        template<typename... Args, typename = typename std::enable_if<sizeof...(Args) == D && ISGREATER(D, 4), void>::type>
        PREONMATH_DEVICE vec(Args&&... vals)
            : m_Data({{vals...}})
        {
        }
#endif  // PREONMATH_COMPILER_MSVC

        //! Builds a vector of dimension D by taking two inputs: a vector of dimension D-1 and the missing component.
        template<PrMathSize E>
        PREONMATH_FORCEINLINE explicit vec(const vec<E, T>& v, T component, typename std::enable_if<E == D - 1, T>::type* = 0)
        {
            StaticFor<0, E>([&](PrMathSize d) { m_Data[d] = v[d]; });
            m_Data[E] = component;
        }

        //! Conversion between vectors with same dimension but different types.
        template<typename S>
        PREONMATH_FORCEINLINE explicit vec(const vec<D, S>& vec)
        {
            setAll([&](PrMathSize d) { return Simd::scalar_cast<S, T>(vec[d]); });
        }

        //! Builds the vector by executing the given function f(PrMathSize d) -> T on all elements of the vector, indexed by d.
        template<typename Lambda>
        PREONMATH_FORCEINLINE explicit vec(const Lambda& f, enable_if_not_convertible<Lambda, T> = 0)
        {
            setAll(f);
        }

        PREONMATH_DEVICE void setZero()
        {
            // Calling the AssignLoopUnroller directly as setAll results in a compile error for vec<3, matrix<2, 3, float>>.
            AssignLoopUnroller<0, D>::step(this, [](PrMathSize) { return scalar_zero<T>(); });
            // setAll([&] (PrMathSize) { return scalar_zero<T>(); });
        }

        /** Getters **/

        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 1, T>::type x() const
        {
            return element(0);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 1, T>::type& x()
        {
            return element(0);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 2, T>::type y() const
        {
            return element(1);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 2, T>::type& y()
        {
            return element(1);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 3, T>::type z() const
        {
            return element(2);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 3, T>::type& z()
        {
            return element(2);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 4, T>::type w() const
        {
            return element(3);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 4, T>::type& w()
        {
            return element(3);
        }

        PREONMATH_FORCEINLINE const T& operator()(PrMathSize d) const
        {
            return m_Data[d];
        }
        PREONMATH_FORCEINLINE T& operator()(PrMathSize d)
        {
            return m_Data[d];
        }

        PREONMATH_FORCEINLINE T& operator[](PrMathSize d)
        {
            return m_Data[d];
        }
        PREONMATH_FORCEINLINE const T& operator[](PrMathSize d) const
        {
            return m_Data[d];
        }

        PREONMATH_FORCEINLINE const T& element(PrMathSize d) const
        {
            return m_Data[d];
        }
        PREONMATH_FORCEINLINE T& element(PrMathSize d)
        {
            return m_Data[d];
        }

        const std::array<T, D>& data() const
        {
            return m_Data;
        }

        //! Returns a vector of dimension E, copying the data. E must be lesser or equal D.
        template<PrMathSize E>
        PREONMATH_DEVICE typename std::enable_if<E <= D, vec<E, T>>::type subVector() const
        {
            return vec<E, T>([this](int i) -> const T& { return m_Data[i]; });
        }

        /** Setters **/

        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 1, void>::type setX(T x)
        {
            setElement(0, x);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 2, void>::type setY(T y)
        {
            setElement(1, y);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 3, void>::type setZ(T z)
        {
            setElement(2, z);
        }
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E >= 4, void>::type setW(T w)
        {
            setElement(3, w);
        }

        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E == 2, void>::type set(T x, T y)
        {
            setX(x);
            setY(y);
        }  // TODO
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E == 3, void>::type set(T x, T y, T z)
        {
            setX(x);
            setY(y);
            setZ(z);
        }  // TODO
        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E == 4, void>::type set(T x, T y, T z, T w)
        {
            setX(x);
            setY(y);
            setZ(z);
            setW(w);
        }  // TODO
        template<typename... Args, PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E == sizeof...(Args) && ISGREATER(E, 4), void>::type set(Args&&... vals)
        {
            m_Data = {{vals...}};
        }

        PREONMATH_DEVICE void setElement(PrMathSize d, T value)
        {
            m_Data[d] = value;
        }

        //! Sets all the elements to the given fill value
        PREONMATH_DEVICE void setAll(T fillValue)
        {
            return setAll([&](PrMathSize) { return fillValue; });
        }

        //! Sets the elements by evaluating the given lambda expression f(PrMathSize d) -> T for all elements of the vector, indexed by d.
        template<typename Lambda>
        PREONMATH_FORCEINLINE void setAll(const Lambda& f, enable_if_not_convertible<Lambda, T> = 0)
        {
            AssignLoopUnroller<0, D>::step(this, f);
        }

        //! Executes the given lambda expression f(PrMathSize d) on all elements of the vector, indexed by d.
        template<typename Lambda>
        PREONMATH_FORCEINLINE void forAll(const Lambda& f) const
        {
            StaticFor<0, D>(f);
        }

        /** Static predefined vectors **/

        PREONMATH_DEVICE static vec<D, T> zero()
        {
            vec<D, T> v;
            v.setZero();
            return v;
        }

        /** Length **/

        PREONMATH_DEVICE T length() const
        {
            return Preon::Math::sqrt(lengthSquared());
        }
        PREONMATH_DEVICE product_t<T, T> lengthSquared() const
        {
            return dotProduct(*this, *this);
        }
        PREONMATH_DEVICE product_t<T, T> squaredDistance(const vec<D, T>& other) const
        {
            return ((*this) - other).lengthSquared();
        }  // TODO: move it to an auxiliary class

        /** Normalization **/

        template<typename E = T>
        PREONMATH_DEVICE void normalize(enable_if_non_simd<E, void*> = nullptr)
        {
            T length = this->length();
            THROW_EXCEPTION(IsZero(length), std::domain_error("Normalizing a zero vector is undefined"))
            (*this) /= length;
        }

#ifdef PREONMATH_ENABLE_SIMD
        template<typename E = T>
        PREONMATH_DEVICE void pseudoNormalize(enable_if_simd<E, void*> = nullptr)
        {
            T len = length();
            (*this) /= len;

            // This code ensures that zero vectors stay zero.
            T mask = (len > Simd::zero<typename Simd::Scalar<T>::type>());
            setAll([&](PrMathSize d) { return Simd::and_mask(m_Data[d], mask); });
        }
#endif

        PREONMATH_DEVICE vec<D, T> normalized() const
        {
            vec<D, T> out = (*this);
            out.normalize();
            return out;
        }

        //! Normalizes the vector, but does not throw an exception if the length of the vector is zero. Zero vectors are not changed by this operation.
        template<typename E = T>
        PREONMATH_DEVICE void pseudoNormalize(enable_if_non_simd<E, void*> = nullptr)
        {
            T len = length();
            if (len > scalar_zero<T>())
                (*this) /= len;
        }

        //! Returns the normalized vector, but does not throw an exception if the length of the vector is zero. In this case, a zero vector is returned.
        PREONMATH_DEVICE vec<D, T> pseudoNormalized() const
        {
            vec<D, T> v = *this;
            v.pseudoNormalize();
            return v;
        }

        /** Operators **/

        PREONMATH_DEVICE vec<D, T>& operator+=(const vec<D, T>& vector)
        {
            setAll([&](PrMathSize d) { return m_Data[d] + vector[d]; });
            return *this;
        }

        PREONMATH_DEVICE vec<D, T>& operator-=(const vec<D, T>& vector)
        {
            setAll([&](PrMathSize d) { return m_Data[d] - vector[d]; });
            return *this;
        }

        template<typename F, typename = enable_if_scalar<F>>
        PREONMATH_DEVICE vec<D, T>& operator*=(F scalar)
        {
            setAll([&](PrMathSize d) { return m_Data[d] * scalar; });
            return *this;
        }

        template<typename F, typename = enable_if_scalar<F>>
        PREONMATH_DEVICE vec<D, T>& operator/=(F scalar)
        {
            setAll([&](PrMathSize d) { return m_Data[d] / scalar; });
            return *this;
        }

        PREONMATH_DEVICE operator T() const
        {
            static_assert(D == 1, "Scalar operator only available for 1-dimensional vectors");
            return m_Data[0];
        }

        template<typename U = T, typename V = T>
        PREONMATH_DEVICE vec<D, U> elementWiseMultiplication(const vec<D, V>& factor) const
        {
            return vec<D, U>([&](PrMathSize d) { return (*this)[d] * factor[d]; });
        }

        template<typename U = T, typename V = T>
        PREONMATH_DEVICE vec<D, U> elementWiseDivision(const vec<D, V>& divisor) const
        {
            return vec<D, U>([&](PrMathSize d) { return (*this)[d] / divisor[d]; });
        }

        template<typename T_Vec1, typename T_Vec2>
        static PREONMATH_FORCEINLINE product_t<T_Vec1, T_Vec2> dotProduct(const vec<D, T_Vec1>& v1, const vec<D, T_Vec2>& v2)
        {
            product_t<T_Vec1, T_Vec2> out = v1[0] * v2[0];
            StaticFor<1, D>([&](PrMathSize d) { out = out + product_t<T_Vec1, T_Vec2>(v1[d] * v2[d]); });
            return out;
        }

        static PREONMATH_FORCEINLINE product_t<T, T> dotProduct(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return dotProduct<T, T>(v1, v2);
        }

        /** Cross product **/

        template<PrMathSize E = D>
        static PREONMATH_DEVICE typename std::enable_if<E == 3, vec<D, T>>::type crossProduct(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return crossProduct<T, T, D>(v1, v2);
        }

        template<typename Out, typename In, PrMathSize E = D>
        static PREONMATH_DEVICE typename std::enable_if<E == 3, vec<D, T>>::type crossProduct(const vec<D, T>& v1, const vec<D, In>& v2)
        {
            return vec<D, Out>(v1.y() * v2.z() - v1.z() * v2.y(), v1.z() * v2.x() - v1.x() * v2.z(), v1.x() * v2.y() - v1.y() * v2.x());
        }

        /** Orthogonal direction **/

        template<PrMathSize E = D>
        PREONMATH_DEVICE typename std::enable_if<E == 3 && !is_simd_scalar<T>::value, vec<D, T>>::type orthogonalDirection() const
        {
            vec<D, T> helperVec(0, 0, 0);
            helperVec[(indexOfMaxAbsComponent() + 1) % 3] = 1;
            vec<D, T> orthogonalDir = vec<D, T>::crossProduct(*this, helperVec);
            orthogonalDir.normalize();
            return orthogonalDir;
        }

        /** Elementwise operations **/

        //! Calls the "value" member function on each element.
        PREONMATH_DEVICE auto value() const
        {
            return vec<D, decltype(std::declval<T>().value())>([this](PrMathSize d) { return element(d).value(); });
        }

        //! Only divides those elements whose divisor is different from 0. Else, set the result to 0.
        template<typename E = T>
        PREONMATH_DEVICE vec<D, T> elementWisePseudoDivision(const enable_if_non_simd<E, vec<D, T>>& divisor) const
        {
            return vec<D, T>([&](PrMathSize d) {
                // This looks so complicated to prevent MSVC from throwing a divide by 0 warning.
                bool iszero = IsZero<T>(divisor[d]);
                return (iszero ? T(0) : m_Data[d]) / (iszero ? T(1) : divisor[d]);
            });
        }

        //! The element wise minimum of two vectors.
        static PREONMATH_DEVICE vec<D, T> min(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return vec<D, T>([&](PrMathSize d) { return std::min(v1[d], v2[d]); });
        }

        //! The element wise maximum of two vectors.
        static PREONMATH_DEVICE vec<D, T> max(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return vec<D, T>([&](PrMathSize d) { return std::max(v1[d], v2[d]); });
        }

        //! The element wise absolute value.
        static PREONMATH_DEVICE vec<D, T> abs(const vec<D, T>& vector)
        {
            return vec<D, T>([&](PrMathSize d) { return std::abs(vector[d]); });
        }

        //! The element wise square root.
        static PREONMATH_DEVICE vec<D, T> sqrt(const vec<D, T>& vector)
        {
            return vec<D, T>([&](PrMathSize d) { return std::sqrt(vector[d]); });
        }

        //! The element wise power.
        template<typename E = T>
        static PREONMATH_DEVICE vec<D, T> pow(const enable_if_non_simd<E, vec<D, E>>& vector, T exponent)
        {
            return vec<D, T>([&](PrMathSize d) { return std::pow(vector[d], exponent); });
        }

        //! The element wise reciprocal.
        template<typename E = T>
        static PREONMATH_DEVICE vec<D, T> reciprocal(const enable_if_non_simd<E, vec<D, E>>& vector)
        {
            return vec<D, T>([&](PrMathSize d) { return 1 / vector[d]; });
        }

#ifdef PREONMATH_ENABLE_SIMD
        //! SIMD specific masking.
        template<typename X = T>
        static vec<D, X> and_mask(const vec<D, X>& v1, const vec<D, X>& v2)
        {
            return vec<D, X>([&](PrMathSize d) { return Simd::and_mask(v1[d], v2[d]); });
        }

        template<typename X = T>
        static vec<D, X> or_mask(const vec<D, X>& v1, const vec<D, X>& v2)
        {
            return vec<D, X>([&](PrMathSize d) { return Simd::or_mask(v1[d], v2[d]); });
        }
#endif

        /** Other stuff **/

        template<typename E = T>
        PrMathSize PREONMATH_DEVICE indexOfMinComponent(enable_if_non_simd<E, void*> = nullptr) const
        {
            PrMathSize index = 0;
            T currentMin = element(0);
            for (PrMathSize d = 1; d < D; ++d)
            {
                if (element(d) < currentMin)
                {
                    index = d;
                    currentMin = element(d);
                }
            }
            return index;
        }
        template<typename E = T>
        PrMathSize PREONMATH_DEVICE indexOfMinAbsComponent(enable_if_non_simd<E, void*> = nullptr) const
        {
            PrMathSize index = 0;
            T currentAbsMin = std::abs(element(0));
            for (PrMathSize d = 1; d < D; ++d)
            {
                if (std::abs(element(d)) < currentAbsMin)
                {
                    index = d;
                    currentAbsMin = std::abs(element(d));
                }
            }
            return index;
        }
        template<typename E = T>
        PrMathSize PREONMATH_DEVICE indexOfMaxComponent(enable_if_non_simd<E, void*> = nullptr) const
        {
            PrMathSize index = 0;
            T currentMax = element(0);
            for (PrMathSize d = 1; d < D; ++d)
            {
                if (element(d) > currentMax)
                {
                    index = d;
                    currentMax = element(d);
                }
            }
            return index;
        }
        template<typename E = T>
        PrMathSize PREONMATH_DEVICE indexOfMaxAbsComponent(enable_if_non_simd<E, void*> = nullptr) const
        {
            PrMathSize index = 0;
            T currentAbsMax = std::abs(element(0));
            for (PrMathSize d = 1; d < D; ++d)
            {
                if (std::abs(element(d)) > currentAbsMax)
                {
                    index = d;
                    currentAbsMax = std::abs(element(d));
                }
            }
            return index;
        }

#ifdef PREONMATH_QT_INTEGRATION
        template<typename E = T>
        QString toUtf8(enable_if_non_simd<E, void*> = nullptr) const
        {
            std::stringstream ss;
            ss << *this;
            return QString::fromStdString(ss.str());
        }
#endif  // PREONMATH_QT_INTEGRATION

    private:
        std::array<T, D> m_Data;
    };

    template<PrMathSize D, typename T>
    using vec_if_non_simd = typename std::enable_if<!is_simd_scalar<T>::value, vec<D, T>>::type;

    template<PrMathSize D, typename T>
    using vec_if_simd = typename std::enable_if<is_simd_scalar<T>::value, vec<D, T>>::type;

    /** Overloaded non-member operators **/

    template<PrMathSize D, typename T>
    PREONMATH_DEVICE bool operator==(const vec_if_non_simd<D, T>& v1, const vec<D, T>& v2)
    {
        bool result = true;
        for (PrMathSize d = 0; d < D; d++)
            result &= v1.element(d) == v2.element(d);
        return result;
    }

    template<PrMathSize D, typename T>
    PREONMATH_DEVICE bool operator!=(const vec_if_non_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return !(v1 == v2);
    }

    template<PrMathSize D, typename T, typename F, typename = enable_if_scalar<F>>
    PREONMATH_DEVICE vec<D, T> operator*(const F scalar, const vec<D, T>& v)
    {
        return vec<D, T>([&](PrMathSize d) { return scalar * v[d]; });
    }

    template<PrMathSize D, typename T, typename F, typename = enable_if_scalar<F>>
    PREONMATH_DEVICE vec<D, T> operator*(const vec<D, T>& v, const F scalar)
    {
        return scalar * v;
    }

    template<PrMathSize D, typename T, typename F, typename = enable_if_scalar<F>>
    PREONMATH_DEVICE vec<D, T> operator/(const vec<D, T>& v, const F scalar)
    {
        return vec<D, T>([&](PrMathSize d) { return v[d] / scalar; });
    }

    template<PrMathSize D, typename T>
    PREONMATH_DEVICE vec<D, T> operator+(const vec<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](PrMathSize d) { return v1[d] + v2[d]; });
    }

    template<PrMathSize D, typename T>
    PREONMATH_DEVICE vec<D, T> operator-(const vec<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](PrMathSize d) { return v1[d] - v2[d]; });
    }

    template<PrMathSize D, typename T_Vec1, typename T_Vec2>
    PREONMATH_DEVICE product_t<T_Vec1, T_Vec2> operator*(const vec<D, T_Vec1>& v1, const vec<D, T_Vec2>& v2)
    {
        return vec<D, product_t<T_Vec1, T_Vec2>>::template dotProduct<T_Vec1, T_Vec2>(v1, v2);
    }

    template<PrMathSize D, typename T>
    PREONMATH_DEVICE product_t<T, T> operator*(const vec<D, T>& v1, const vec<D, T>& v2)
    {
        return operator*<D, T, T>(v1, v2);
    }

    template<PrMathSize D, typename T>
    PREONMATH_DEVICE vec<D, T> operator-(const vec<D, T>& v)
    {
        return vec<D, T>([&](PrMathSize d) { return -v[d]; });
    }

    //! SIMD element wise operators.
    template<PrMathSize D, typename T>
    PREONMATH_DEVICE vec<D, T> operator>(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](PrMathSize d) -> T { return v1[d] > v2[d]; });
    }
    template<PrMathSize D, typename T>
    PREONMATH_DEVICE vec<D, T> operator>=(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](PrMathSize d) -> T { return v1[d] >= v2[d]; });
    }
    template<PrMathSize D, typename T>
    PREONMATH_DEVICE vec<D, T> operator<(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](PrMathSize d) -> T { return v1[d] < v2[d]; });
    }
    template<PrMathSize D, typename T>
    PREONMATH_DEVICE vec<D, T> operator<=(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](PrMathSize d) -> T { return v1[d] <= v2[d]; });
    }
    template<PrMathSize D, typename T>
    PREONMATH_DEVICE vec<D, T> operator==(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](PrMathSize d) -> T { return v1[d] == v2[d]; });
    }

    //! During comparison x is most significant, followed by y etc.
    template<PrMathSize D, typename T>
    PREONMATH_DEVICE bool operator<(const vec_if_non_simd<D, T>& v1, const vec<D, T>& v2)
    {
        for (PrMathSize d = 0; d < D; d++)
        {
            if (v1[d] != v2[d])
                return v1[d] < v2[d];
        }
        return false;
    }

    template<PrMathSize D, typename T>
    std::stringstream& operator<<(std::stringstream& stream, const vec<D, T>& vector)
    {
        stream << "(";
        for (PrMathSize d = 0; d < D - 1; ++d)
        {
            stream << vector.element(d) << ", ";
        }
        stream << vector.element(D - 1) << ")";

        return stream;
    }

    template<PrMathSize D, typename T>
    std::ostream& operator<<(std::ostream& stream, const vec<D, T>& vector)
    {
        std::stringstream ss;
        ss << vector;
        stream << ss.str();
        return stream;
    }

    /** Approximative equality **/

    template<PrMathSize D, typename T>
    PREONMATH_DEVICE bool AreEqual(const vec<D, T>& v1, const vec<D, T>& v2)
    {
        bool result = true;
        for (PrMathSize d = 0; d < D; d++)
            result &= AreEqual(v1.element(d), v2.element(d));
        return result;
    }

#ifdef PREONMATH_QT_INTEGRATION
    /** Qt integration **/

    template<PrMathSize D, typename T>
    QDebug operator<<(QDebug dbg, const vec<D, T>& vector)
    {
        dbg.nospace() << "vec" << D << "(";
        for (PrMathSize d = 0; d < D - 1; d++)
            dbg.nospace() << vector.element(d) << ", ";
        dbg.nospace() << vector.element(D - 1) << ')';
        return dbg.space();
    }

    template<PrMathSize D, typename T>
    QDataStream& operator<<(QDataStream& stream, const vec<D, T>& vector)
    {
        for (PrMathSize d = 0; d < D; d++)
            stream << vector.element(d);
        return stream;
    }

    template<PrMathSize D, typename T>
    QDataStream& operator>>(QDataStream& stream, vec<D, T>& vector)
    {
        for (PrMathSize d = 0; d < D; d++)
            stream >> vector.element(d);
        return stream;
    }
#endif  // PREONMATH_QT_INTEGRATION
}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
template<PrMathSize D, typename T>
bool qFuzzyCompare(const Preon::Math::vec<D, T>& v1, const Preon::Math::vec<D, T>& v2)
{
    bool result = true;
    for (PrMathSize d = 0; d < D; d++)
        result &= qFuzzyCompare(v1.element(d), v2.element(d));
    return result;
}
#endif  // PREONMATH_QT_INTEGRATION
