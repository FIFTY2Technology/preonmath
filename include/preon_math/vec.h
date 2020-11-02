/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "basics.h"
#include "scalar_simd.h"
#include "simd_scalar_wrapper_std.h"
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
    using std::size_t;

    namespace
    {
        template<typename T, typename U>
        using enable_if_not_convertible = typename std::enable_if<!std::is_convertible<T, U>::value>::type*;

        template<typename T>
        using enable_if_scalar = std::enable_if_t<std::is_arithmetic<T>::value || is_simd_scalar<T>::value>;
    }  // namespace

// Replace with std::isgreater when we use c++14.
#define ISGREATER(a, b) (a > b)

    template<size_t D, typename T>
    class vec
    {
    private:
        //! legacyInitZeroFunc initializes the components of the vector to zero for non-SIMD types.
        //! This is necessary for legacy reasons, but in the future, we want to move away from it.
        template<typename E>
        void legacyInitZeroFunc(enable_if_non_simd<E, void*> = nullptr)
        {
            setZero();
        }

        //! Do nothing for SIMD types here.
        template<typename E>
        void legacyInitZeroFunc(enable_if_simd<E, void*> = nullptr)
        {
        }

        template<size_t Begin, size_t End>
        struct AssignLoopUnroller
        {
            template<typename Lambda>
            PREONMATH_FORCEINLINE static void step(vec<D, T>* v, const Lambda& func)
            {
                (*v)[Begin] = func(Begin);
                AssignLoopUnroller<Begin + 1, End>::step(v, func);
            }
        };
        template<size_t End>
        struct AssignLoopUnroller<End, End>
        {
            template<typename Lambda>
            PREONMATH_FORCEINLINE static void step(vec<D, T>*, const Lambda&)
            {
            }
        };

        template<typename E = T>
        PREONMATH_FORCEINLINE static enable_if_non_simd<E, E> scalar_zero()
        {
            return E(0);
        }

        template<typename E = T>
        PREONMATH_FORCEINLINE static enable_if_simd<E, E> scalar_zero()
        {
            return Simd::zero<typename Simd::Scalar<T>::type>();
        }

    public:
        /** Constructors **/
        vec() { legacyInitZeroFunc<T>(); }
        PREONMATH_FORCEINLINE explicit vec(T fillValue)
        {
            setAll([&](size_t) { return fillValue; });
        }

        template<size_t E = D>
        PREONMATH_FORCEINLINE vec(T x, T y, typename std::enable_if<E == 2, T>::type* = 0)
        {
            set(x, y);
        }
        template<size_t E = D>
        PREONMATH_FORCEINLINE vec(T x, T y, T z, typename std::enable_if<E == 3, T>::type* = 0)
        {
            set(x, y, z);
        }
        template<size_t E = D>
        PREONMATH_FORCEINLINE vec(T x, T y, T z, T w, typename std::enable_if<E == 4, T>::type* = 0)
        {
            set(x, y, z, w);
        }
        // MSVC currently has problems with this variadic args cnstructor, so we only implement the ones we need currently.
#ifdef PREONMATH_COMPILER_MSVC
        template<size_t E = D>
        PREONMATH_FORCEINLINE vec(T a, T b, T c, T d, T e, typename std::enable_if<E == 5, T>::type* = 0)
        {
            set(a, b, c, d, e);
        }
        template<size_t E = D>
        PREONMATH_FORCEINLINE vec(T a, T b, T c, T d, T e, T f, typename std::enable_if<E == 6, T>::type* = 0)
        {
            set(a, b, c, d, e, f);
        }
#else
        template<typename... Args, typename = typename std::enable_if<sizeof...(Args) == D && ISGREATER(D, 4), void>::type>
        vec(Args&&... vals)
            : m_Data({{vals...}})
        {
        }
#endif  // PREONMATH_COMPILER_MSVC

        //! Builds a vector of dimension D by taking two inputs: a vector of dimension D-1 and the missing component.
        template<size_t E>
        PREONMATH_FORCEINLINE explicit vec(const vec<E, T>& v, T component, typename std::enable_if<E == D - 1, T>::type* = 0)
        {
            StaticFor<0, E>([&](size_t d) { m_Data[d] = v[d]; });
            m_Data[E] = component;
        }

        //! Conversion between vectors with same dimension but different types.
        template<typename S>
        PREONMATH_FORCEINLINE explicit vec(const vec<D, S>& vec)
        {
            setAll([&](size_t d) { return Simd::scalar_cast<S, T>(vec[d]); });
        }

        //! Builds the vector by executing the given function f(size_t d) -> T on all elements of the vector, indexed by d.
        template<typename Lambda>
        PREONMATH_FORCEINLINE explicit vec(const Lambda& f, enable_if_not_convertible<Lambda, T> = 0)
        {
            setAll(f);
        }

        void setZero()
        {
            // Calling the AssignLoopUnroller directly as setAll results in a compile error for vec<3, matrix<2, 3, float>>.
            AssignLoopUnroller<0, D>::step(this, [](size_t) { return scalar_zero<T>(); });
            // setAll([&] (size_t) { return scalar_zero<T>(); });
        }

        /** Getters **/

        template<size_t E = D>
        typename std::enable_if<E >= 1, T>::type x() const
        {
            return element(0);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 1, T>::type& x()
        {
            return element(0);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 2, T>::type y() const
        {
            return element(1);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 2, T>::type& y()
        {
            return element(1);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 3, T>::type z() const
        {
            return element(2);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 3, T>::type& z()
        {
            return element(2);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 4, T>::type w() const
        {
            return element(3);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 4, T>::type& w()
        {
            return element(3);
        }

        PREONMATH_FORCEINLINE const T& operator()(size_t d) const { return m_Data[d]; }
        PREONMATH_FORCEINLINE T& operator()(size_t d) { return m_Data[d]; }

        PREONMATH_FORCEINLINE T& operator[](size_t d) { return m_Data[d]; }
        PREONMATH_FORCEINLINE const T& operator[](size_t d) const { return m_Data[d]; }

        PREONMATH_FORCEINLINE const T& element(size_t d) const { return m_Data[d]; }
        PREONMATH_FORCEINLINE T& element(size_t d) { return m_Data[d]; }

        //! Returns a vector of dimension E, copying the data. E must be lesser or equal D.
        template<size_t E>
        typename std::enable_if<E <= D, vec<E, T>>::type subVector() const
        {
            return vec<E, T>([this](int i) -> const T& { return m_Data[i]; });
        }

        /** Setters **/

        template<size_t E = D>
        typename std::enable_if<E >= 1, void>::type setX(T x)
        {
            setElement(0, x);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 2, void>::type setY(T y)
        {
            setElement(1, y);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 3, void>::type setZ(T z)
        {
            setElement(2, z);
        }
        template<size_t E = D>
        typename std::enable_if<E >= 4, void>::type setW(T w)
        {
            setElement(3, w);
        }

        template<size_t E = D>
        typename std::enable_if<E == 2, void>::type set(T x, T y)
        {
            setX(x);
            setY(y);
        }  // TODO
        template<size_t E = D>
        typename std::enable_if<E == 3, void>::type set(T x, T y, T z)
        {
            setX(x);
            setY(y);
            setZ(z);
        }  // TODO
        template<size_t E = D>
        typename std::enable_if<E == 4, void>::type set(T x, T y, T z, T w)
        {
            setX(x);
            setY(y);
            setZ(z);
            setW(w);
        }  // TODO
        template<typename... Args, size_t E = D>
        typename std::enable_if<E == sizeof...(Args) && ISGREATER(E, 4), void>::type set(Args&&... vals)
        {
            m_Data = {{vals...}};
        }

        void setElement(size_t d, T value) { m_Data[d] = value; }

        //! Sets all the elements to the given fill value
        void setAll(T fillValue)
        {
            return setAll([&](size_t) { return fillValue; });
        }

        //! Sets the elements by evaluating the given lambda expression f(size_t d) -> T for all elements of the vector, indexed by d.
        template<typename Lambda>
        PREONMATH_FORCEINLINE void setAll(const Lambda& f, enable_if_not_convertible<Lambda, T> = 0)
        {
            AssignLoopUnroller<0, D>::step(this, f);
        }

        //! Executes the given lambda expression f(size_t d) on all elements of the vector, indexed by d.
        template<typename Lambda>
        PREONMATH_FORCEINLINE void forAll(const Lambda& f) const
        {
            StaticFor<0, D>(f);
        }

        /** Static predefined vectors **/

        static vec<D, T> zero()
        {
            vec<D, T> v;
            v.setZero();
            return v;
        }

        /** Length **/

        T length() const { return SimdWrapper::sqrt(lengthSquared()); }
        T lengthSquared() const { return dotProduct(*this, *this); }
        T squaredDistance(const vec<D, T>& other) const { return ((*this) - other).lengthSquared(); }  // TODO: move it to an auxiliary class

        /** Normalization **/

        template<typename E = T>
        void normalize(enable_if_non_simd<E, void*> = nullptr)
        {
            T length = this->length();
            THROW_EXCEPTION(IsZero(length), std::domain_error("Normalizing a zero vector is undefined"))
            (*this) /= length;
        }

        template<typename E = T>
        void pseudoNormalize(enable_if_simd<E, void*> = nullptr)
        {
            T len = length();
            (*this) /= len;

            // This code ensures that zero vectors stay zero.
            T mask = (len > Simd::zero<typename Simd::Scalar<T>::type>());
            setAll([&](size_t d) { return Simd::and_mask(m_Data[d], mask); });
        }

        vec<D, T> normalized() const
        {
            vec<D, T> out = (*this);
            out.normalize();
            return out;
        }

        //! Normalizes the vector, but does not throw an exception if the length of the vector is zero. Zero vectors are not changed by this operation.
        template<typename E = T>
        void pseudoNormalize(enable_if_non_simd<E, void*> = nullptr)
        {
            T len = length();
            if (len > 0)
                (*this) /= len;
        }

        //! Returns the normalized vector, but does not throw an exception if the length of the vector is zero. In this case, a zero vector is returned.
        vec<D, T> pseudoNormalized() const
        {
            vec<D, T> v = *this;
            v.pseudoNormalize();
            return v;
        }

        /** Operators **/

        vec<D, T>& operator+=(const vec<D, T>& vector)
        {
            setAll([&](size_t d) { return m_Data[d] + vector[d]; });
            return *this;
        }

        vec<D, T>& operator-=(const vec<D, T>& vector)
        {
            setAll([&](size_t d) { return m_Data[d] - vector[d]; });
            return *this;
        }

        vec<D, T>& operator*=(T scalar)
        {
            setAll([&](size_t d) { return m_Data[d] * scalar; });
            return *this;
        }

        vec<D, T>& operator/=(T scalar)
        {
            setAll([&](size_t d) { return m_Data[d] / scalar; });
            return *this;
        }

        template<typename U = T, typename V = T>
        vec<D, U> elementWiseMultiplication(const vec<D, V>& factor) const
        {
            return vec<D, U>([&](size_t d) { return (*this)[d] * factor[d]; });
        }

        vec<D, T> elementWiseDivision(const vec<D, T>& divisor) const
        {
            return vec<D, T>([&](size_t d) { return (*this)[d] / divisor[d]; });
        }

        template<typename T_Out = T, typename T_Vec1, typename T_Vec2>
        static PREONMATH_FORCEINLINE T_Out dotProduct(const vec<D, T_Vec1>& v1, const vec<D, T_Vec2>& v2)
        {
            T_Out out = v1[0] * v2[0];
            StaticFor<1, D>([&](size_t d) { out = out + (v1[d] * v2[d]); });
            return out;
        }

        static PREONMATH_FORCEINLINE T dotProduct(const vec<D, T>& v1, const vec<D, T>& v2) { return dotProduct<T, T, T>(v1, v2); }

        /** Cross product **/

        template<size_t E = D>
        static typename std::enable_if<E == 3, vec<D, T>>::type crossProduct(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return crossProduct<T, T, D>(v1, v2);
        }

        template<typename Out, typename In, size_t E = D>
        static typename std::enable_if<E == 3, vec<D, T>>::type crossProduct(const vec<D, T>& v1, const vec<D, In>& v2)
        {
            return vec<D, Out>(v1.y() * v2.z() - v1.z() * v2.y(), v1.z() * v2.x() - v1.x() * v2.z(), v1.x() * v2.y() - v1.y() * v2.x());
        }

        /** Orthogonal direction **/

        template<size_t E = D>
        typename std::enable_if<E == 3 && !is_simd_scalar<T>::value, vec<D, T>>::type orthogonalDirection() const
        {
            vec<D, T> helperVec(0, 0, 0);
            helperVec[(indexOfMaxAbsComponent() + 1) % 3] = 1;
            vec<D, T> orthogonalDir = vec<D, T>::crossProduct(*this, helperVec);
            orthogonalDir.normalize();
            return orthogonalDir;
        }

        /** Elementwise operations **/
        //! Only divides those elements whose divisor is different form 0. Else, set the result to 0.
        template<typename E = T>
        vec<D, T> elementWisePseudoDivision(const enable_if_non_simd<E, vec<D, T>>& divisor) const
        {
            return vec<D, T>([&](size_t d) {
                // This looks so complicated to prevent MSVC from throwing a divide by 0 warning.
                bool iszero = IsZero<T>(divisor[d]);
                return (iszero ? T(0) : m_Data[d]) / (iszero ? T(1) : divisor[d]);
            });
        }

        //! The element wise minimum of two vectors.
        static vec<D, T> min(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return vec<D, T>([&](size_t d) { return std::min(v1[d], v2[d]); });
        }

        //! The element wise maximum of two vectors.
        static vec<D, T> max(const vec<D, T>& v1, const vec<D, T>& v2)
        {
            return vec<D, T>([&](size_t d) { return std::max(v1[d], v2[d]); });
        }

        //! The element wise absolute value.
        static vec<D, T> abs(const vec<D, T>& vector)
        {
            return vec<D, T>([&](size_t d) { return std::abs(vector[d]); });
        }

        //! The element wise square root.
        static vec<D, T> sqrt(const vec<D, T>& vector)
        {
            return vec<D, T>([&](size_t d) { return std::sqrt(vector[d]); });
        }

        //! The element wise power.
        template<typename E = T>
        static vec<D, T> pow(const enable_if_non_simd<E, vec<D, E>>& vector, T exponent)
        {
            return vec<D, T>([&](size_t d) { return std::pow(vector[d], exponent); });
        }

        //! The element wise reciprocal.
        template<typename E = T>
        static vec<D, T> reciprocal(const enable_if_non_simd<E, vec<D, E>>& vector)
        {
            return vec<D, T>([&](size_t d) { return 1 / vector[d]; });
        }

        //! SIMD specific masking.
        template<typename X = T>
        static vec<D, X> and_mask(const vec<D, X>& v1, const vec<D, X>& v2)
        {
            return vec<D, X>([&](size_t d) { return Simd::and_mask(v1[d], v2[d]); });
        }

        template<typename X = T>
        static vec<D, X> or_mask(const vec<D, X>& v1, const vec<D, X>& v2)
        {
            return vec<D, X>([&](size_t d) { return Simd::or_mask(v1[d], v2[d]); });
        }

        /** Other stuff **/

        template<typename E = T>
        size_t indexOfMinComponent(enable_if_non_simd<E, void*> = nullptr) const
        {
            size_t index = 0;
            T currentMin = element(0);
            for (size_t d = 1; d < D; ++d)
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
        size_t indexOfMinAbsComponent(enable_if_non_simd<E, void*> = nullptr) const
        {
            size_t index = 0;
            T currentAbsMin = std::abs(element(0));
            for (size_t d = 1; d < D; ++d)
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
        size_t indexOfMaxComponent(enable_if_non_simd<E, void*> = nullptr) const
        {
            size_t index = 0;
            T currentMax = element(0);
            for (size_t d = 1; d < D; ++d)
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
        size_t indexOfMaxAbsComponent(enable_if_non_simd<E, void*> = nullptr) const
        {
            size_t index = 0;
            T currentAbsMax = std::abs(element(0));
            for (size_t d = 1; d < D; ++d)
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

    template<size_t D, typename T>
    using vec_if_non_simd = typename std::enable_if<!is_simd_scalar<T>::value, vec<D, T>>::type;

    template<size_t D, typename T>
    using vec_if_simd = typename std::enable_if<is_simd_scalar<T>::value, vec<D, T>>::type;

    /** Overloaded non-member operators **/

    template<size_t D, typename T>
    bool operator==(const vec_if_non_simd<D, T>& v1, const vec<D, T>& v2)
    {
        bool result = true;
        for (size_t d = 0; d < D; d++)
            result &= v1.element(d) == v2.element(d);
        return result;
    }

    template<size_t D, typename T>
    bool operator!=(const vec_if_non_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return !(v1 == v2);
    }

    template<size_t D, typename T, typename F, typename = enable_if_scalar<F>>
    vec<D, T> operator*(const F scalar, const vec<D, T>& v)
    {
        return vec<D, T>([&](size_t d) { return scalar * v[d]; });
    }

    template<size_t D, typename T, typename F, typename = enable_if_scalar<F>>
    vec<D, T> operator*(const vec<D, T>& v, const F scalar)
    {
        return scalar * v;
    }

    template<size_t D, typename T, typename F>
    vec<D, T> operator/(const vec<D, T>& v, const F scalar)
    {
        return vec<D, T>([&](size_t d) { return v[d] / scalar; });
    }

    template<size_t D, typename T>
    vec<D, T> operator+(const vec<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](size_t d) { return v1[d] + v2[d]; });
    }

    template<size_t D, typename T>
    vec<D, T> operator-(const vec<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](size_t d) { return v1[d] - v2[d]; });
    }

    template<size_t D, typename T_Out, typename T_Vec1, typename T_Vec2>
    T_Out operator*(const vec<D, T_Vec1>& v1, const vec<D, T_Vec2>& v2)
    {
        return vec<D, T_Out>::template dotProduct<T_Out, T_Vec1, T_Vec2>(v1, v2);
    }

    template<size_t D, typename T>
    T operator*(const vec<D, T>& v1, const vec<D, T>& v2)
    {
        return operator*<D, T, T, T>(v1, v2);
    }

    template<size_t D, typename T>
    vec<D, T> operator-(const vec<D, T>& v)
    {
        return vec<D, T>([&](size_t d) { return -v[d]; });
    }

    //! SIMD element wise operators.
    template<size_t D, typename T>
    vec<D, T> operator>(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](size_t d) -> T { return v1[d] > v2[d]; });
    }
    template<size_t D, typename T>
    vec<D, T> operator>=(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](size_t d) -> T { return v1[d] >= v2[d]; });
    }
    template<size_t D, typename T>
    vec<D, T> operator<(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](size_t d) -> T { return v1[d] < v2[d]; });
    }
    template<size_t D, typename T>
    vec<D, T> operator<=(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](size_t d) -> T { return v1[d] <= v2[d]; });
    }
    template<size_t D, typename T>
    vec<D, T> operator==(const vec_if_simd<D, T>& v1, const vec<D, T>& v2)
    {
        return vec<D, T>([&](size_t d) -> T { return v1[d] == v2[d]; });
    }

    //! During comparison x is most significant, followed by y etc.
    template<size_t D, typename T>
    bool operator<(const vec_if_non_simd<D, T>& v1, const vec<D, T>& v2)
    {
        for (size_t d = 0; d < D; d++)
        {
            if (v1[d] != v2[d])
                return v1[d] < v2[d];
        }
        return false;
    }

    template<size_t D, typename T>
    std::stringstream& operator<<(std::stringstream& stream, const vec<D, T>& vector)
    {
        stream << "(";
        for (size_t d = 0; d < D - 1; ++d)
        {
            stream << vector.element(d) << ", ";
        }
        stream << vector.element(D - 1) << ")";

        return stream;
    }

    template<size_t D, typename T>
    std::ostream& operator<<(std::ostream& stream, const vec<D, T>& vector)
    {
        std::stringstream ss;
        ss << vector;
        stream << ss.str();
        return stream;
    }

    /** Approximative equality **/

    template<size_t D, typename T>
    bool AreEqual(const vec<D, T>& v1, const vec<D, T>& v2)
    {
        bool result = true;
        for (size_t d = 0; d < D; d++)
            result &= AreEqual(v1.element(d), v2.element(d));
        return result;
    }

#ifdef PREONMATH_QT_INTEGRATION
    /** Qt integration **/

    template<size_t D, typename T>
    QDebug operator<<(QDebug dbg, const vec<D, T>& vector)
    {
        dbg.nospace() << "vec" << D << "(";
        for (size_t d = 0; d < D - 1; d++)
            dbg.nospace() << vector.element(d) << ", ";
        dbg.nospace() << vector.element(D - 1) << ')';
        return dbg.space();
    }

    template<size_t D, typename T>
    QDataStream& operator<<(QDataStream& stream, const vec<D, T>& vector)
    {
        for (size_t d = 0; d < D; d++)
            stream << vector.element(d);
        return stream;
    }

    template<size_t D, typename T>
    QDataStream& operator>>(QDataStream& stream, vec<D, T>& vector)
    {
        for (size_t d = 0; d < D; d++)
            stream >> vector.element(d);
        return stream;
    }
#endif  // PREONMATH_QT_INTEGRATION
}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
template<size_t D, typename T>
bool qFuzzyCompare(const Preon::Math::vec<D, T>& v1, const Preon::Math::vec<D, T>& v2)
{
    bool result = true;
    for (size_t d = 0; d < D; d++)
        result &= qFuzzyCompare(v1.element(d), v2.element(d));
    return result;
}
#endif  // PREONMATH_QT_INTEGRATION
