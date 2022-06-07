/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "vec.h"

#include <type_traits>
#include <sstream>
#include <iostream>

namespace Preon
{
namespace Math
{
    using std::size_t;

    template<size_t M, size_t N, typename T>
    class matrix
    {
    public:
        // Constructor.
#ifndef PREONMATH_COMPILER_MSVC
        template<size_t _M = M, size_t _N = N>
        matrix(typename std::enable_if<_M == _N, T>::type* = 0)
        {
            legacyInitIdentityFunc<T>();
        }
        template<size_t _M = M, size_t _N = N>
        matrix(typename std::enable_if<_M != _N, T>::type* = 0)
        {
            legacyInitZeroFunc<T>();
        }
#else
        template<typename = typename std::enable_if<M == N, void>::type>
        matrix()
        {
            legacyInitIdentityFunc<T>();
        }
        template<typename _DUMMY = typename std::enable_if<M != N, void>::type, typename = _DUMMY>
        matrix()
        {
            legacyInitZeroFunc<T>();
        }
#endif  // PREONMATH_COMPILER_MSVC
        explicit matrix(T value)
        {
            fill(value);
        }

        template<typename... Args, typename = typename std::enable_if<sizeof...(Args) == M * N && M * N != 1, void>::type>
        matrix(Args&&... vals)
            : matrix({vals...})
        {
        }

        //! Conversion between matrices with same dimension but different types.
        //! This is especially needed for float -> simd.
        template<typename S>
        PREONMATH_FORCEINLINE explicit matrix(const matrix<M, N, S>& input)
        {
            for (size_t c = 0; c < N; c++)
                for (size_t r = 0; r < M; r++)
                    m_Data[c][r] = Simd::scalar_cast<S, T>(input(r, c));
        }

        template<typename Lambda>
        explicit matrix(const Lambda& lambda, enable_if_not_convertible<Lambda, T> = 0);

        // Column getter & setter
        inline vec<M, T>& column(size_t column);
        const vec<M, T>& column(size_t column) const { return m_Data[column]; }
        inline void setColumn(size_t row, const vec<M, T>& value);

        // Row getter & setter
        inline vec<N, T> row(size_t row) const;
        inline void setRow(size_t row, const vec<N, T>& value);

        // Data access
        T operator()(size_t row, size_t column) const { return element(row, column); }
        T& operator()(size_t row, size_t column) { return element(row, column); }

        T element(size_t row, size_t column) const { return m_Data[column][row]; }
        T& element(size_t row, size_t column) { return m_Data[column][row]; }

        // Filling
        inline void fill(T value);
        inline void setToZero() { fill(scalar_zero<T>()); }
        static matrix<M, N, T> zero()
        {
            matrix<M, N, T> m;
            m.setToZero();
            return m;
        }

        template<size_t _M = M, size_t _N = N>
        inline void setToIdentity(typename std::enable_if<_M == _N, T>::type* = nullptr);
        template<size_t _M = M, size_t _N = N>
        static matrix<M, N, T> identity(typename std::enable_if<_M == _N, T>::type* = nullptr)
        {
            matrix<M, N, T> m;
            m.setToIdentity();
            return m;
        }

        // Transpose
        template<size_t _M = M, size_t _N = N>
        inline void transpose(typename std::enable_if<_M == _N, T>::type* = nullptr);
        inline matrix<N, M, T> transposed() const;

        // Inverse
        template<size_t _M = M, size_t _N = N>
        inline void invert(typename std::enable_if<_M == _N, T>::type* = nullptr)
        {
            *this = inverted();
        }
        template<size_t _M = M, size_t _N = N>
        inline matrix<N, M, T> inverted(typename std::enable_if<_M == _N && _M == 1, T>::type* = nullptr) const;
        template<size_t _M = M, size_t _N = N>
        inline matrix<N, M, T> inverted(typename std::enable_if<_M == _N && _M == 2, T>::type* = nullptr) const;
        template<size_t _M = M, size_t _N = N>
        inline matrix<N, M, T> inverted(typename std::enable_if<_M == _N && _M == 3, T>::type* = nullptr) const;
        template<size_t _M = M, size_t _N = N>
        inline matrix<N, M, T> inverted(typename std::enable_if<_M == _N && _M == 4, T>::type* = nullptr) const;

        // Adjoint and cofactor.
        template<size_t _M = M, size_t _N = N>
        inline matrix<N, M, T> adjugate(typename std::enable_if<_M == _N && _M == 3, T>::type* = nullptr) const;
        template<size_t _M = M, size_t _N = N>
        inline matrix<M, N, T> cofactorMatrix(typename std::enable_if<_M == _N && _M == 3, T>::type* = nullptr) const
        {
            return adjugate().transposed();
        }

        // Reductions
        template<typename E = T>
        inline enable_if_non_simd<E, E> norm(E p = 2) const;

        template<size_t _M = M, size_t _N = N>
        inline T determinant(typename std::enable_if<_M == _N && _M != 2 && _M != 3, T>::type* = nullptr) const;
        template<size_t _M = M, size_t _N = N>
        inline T determinant(typename std::enable_if<_M == 2 && _N == 2, T>::type* = nullptr) const;
        template<size_t _M = M, size_t _N = N>
        inline T determinant(typename std::enable_if<_M == 3 && _N == 3, T>::type* = nullptr) const;

        template<size_t _M = M, size_t _N = N>
        inline T trace(typename std::enable_if<_M == _N, T>::type* = nullptr) const;
        template<size_t _M = M, size_t _N = N>
        vec<M, T> diagonal(typename std::enable_if<_M == _N, T>::type* = nullptr) const;

        // Operators
        inline matrix<M, N, T>& operator+=(const matrix<M, N, T>& other);
        inline matrix<M, N, T>& operator-=(const matrix<M, N, T>& other);
        inline matrix<M, N, T>& operator*=(const matrix<M, N, T>& other);  // TODO: Make it only useable for M == N, else it makes no sense.
        inline matrix<M, N, T>& operator*=(T factor);
        inline matrix<M, N, T>& operator/=(T divisor)
        {
            *this *= (Simd::scalar_cast<float, T>(1) / divisor);
            return *this;
        }

        // Helpers
        template<typename T_Out, typename T_Vec>
        inline T_Out rowVecProduct(size_t row, const vec<N, T_Vec>& vec) const;

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
        std::array<vec<M, T>, N> m_Data;  // Column-major order to match OpenGL.

        // This constructor is private, as the size of the initializer list cannot
        // be checked at compile time, which may lead to unintuitive behavior.
        matrix(const std::initializer_list<T>& values);

        template<typename E>
        void legacyInitZeroFunc(enable_if_non_simd<E, void*> = nullptr)
        {
            setToZero();
        }
        template<typename E>
        void legacyInitZeroFunc(enable_if_simd<E, void*> = nullptr)
        {
        }

        template<typename E>
        void legacyInitIdentityFunc(enable_if_non_simd<E, void*> = nullptr)
        {
            setToIdentity();
        }
        template<typename E>
        void legacyInitIdentityFunc(enable_if_simd<E, void*> = nullptr)
        {
        }
    };

    // Overloaded non-member operators

    template<size_t M, size_t N, typename T>
    inline bool operator==(const matrix<M, N, T>& m1, const matrix<M, N, T>& m2);
    template<size_t M, size_t N, typename T>
    inline bool operator!=(const matrix<M, N, T>& m1, const matrix<M, N, T>& m2)
    {
        return !(m1 == m2);
    }

    template<size_t M, size_t N, typename T>
    inline matrix<M, N, T> operator+(matrix<M, N, T> m1, const matrix<M, N, T>& m2)
    {
        return m1 += m2;
    }
    template<size_t M, size_t N, typename T>
    inline matrix<M, N, T> operator-(matrix<M, N, T> m1, const matrix<M, N, T>& m2)
    {
        return m1 -= m2;
    }
    template<size_t M, size_t N, typename T>
    inline matrix<M, N, T> operator-(const matrix<M, N, T>& m2)
    {
        return matrix<M, N, T>(0) -= m2;
    }

    template<size_t M, size_t N, size_t O, typename T>
    inline matrix<M, O, T> operator*(const matrix<M, N, T>& m1, const matrix<N, O, T>& m2);
    template<size_t M, size_t N, typename T>
    inline matrix<M, N, T> operator*(const T factor, matrix<M, N, T> m)
    {
        return m *= factor;
    }
    template<size_t M, size_t N, typename T>
    inline matrix<M, N, T> operator*(matrix<M, N, T> m, const T factor)
    {
        return factor * m;
    }

    template<typename T, size_t M, size_t N, typename T_Mat, typename T_Vec>
    inline vec<M, T> operator*(const matrix<M, N, T_Mat>& m, const vec<N, T_Vec>& v);
    template<size_t M, size_t N, typename T>
    inline vec<M, T> operator*(const matrix<M, N, T>& m, const vec<N, T>& v)
    {
        return operator*<T, M, N, T, T>(m, v);
    }

    template<size_t M, size_t N, typename T>
    inline matrix<M, N, T> operator/(matrix<M, N, T> m, const T factor)
    {
        return m /= factor;
    }

    template<size_t M, size_t N, typename T>
    inline enable_if_non_simd<T, std::stringstream>& operator<<(std::stringstream&, const matrix<M, N, T>&);

    template<size_t M, size_t N, typename T>
    inline enable_if_non_simd<T, std::ostream>& operator<<(std::ostream&, const matrix<M, N, T>&);

#ifdef PREONMATH_QT_INTEGRATION
    template<size_t M, size_t N, typename T>
    QDebug operator<<(QDebug dbg, const matrix<M, N, T>& m)
    {
        dbg.nospace() << m.toUtf8();
        return dbg.space();
    }
#endif  // PREONMATH_QT_INTEGRATION
}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
template<size_t M, size_t N, typename T>
inline bool qFuzzyCompare(const Preon::Math::matrix<M, N, T>& m1, const Preon::Math::matrix<M, N, T>& m2)
{
    for (size_t c = 0; c < N; c++)
        for (size_t r = 0; r < M; r++)
            if (!qFuzzyCompare(m1(r, c), m2(r, c)))
                return false;
    return true;
}
#endif  // PREONMATH_QT_INTEGRATION

#include "matrix.inl"
