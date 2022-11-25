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
    template<PrMathSize M, PrMathSize N, typename T>
    class matrix
    {
    public:
        PREONMATH_DEVICE matrix()
        {
            if constexpr (M == N)
                legacyInitIdentityFunc<T>();
            else
                legacyInitZeroFunc<T>();
        }

        PREONMATH_DEVICE explicit matrix(T value)
        {
            fill(value);
        }

        template<typename... Args, typename = typename std::enable_if<sizeof...(Args) == M * N && M * N != 1, void>::type>
        PREONMATH_DEVICE matrix(Args&&... vals)
            : matrix({vals...})
        {
        }

        //! Conversion between matrices with same dimension but different types.
        //! This is especially needed for float -> simd.
        template<typename S>
        PREONMATH_FORCEINLINE explicit matrix(const matrix<M, N, S>& input)
        {
            for (PrMathSize c = 0; c < N; c++)
                for (PrMathSize r = 0; r < M; r++)
                    m_Data[c][r] = Simd::scalar_cast<S, T>(input(r, c));
        }

        template<typename Lambda>
        PREONMATH_DEVICE explicit matrix(const Lambda& lambda, enable_if_not_convertible<Lambda, T> = 0);

        // Column getter & setter
        PREONMATH_DEVICE inline vec<M, T>& column(PrMathSize column);
        PREONMATH_DEVICE const vec<M, T>& column(PrMathSize column) const { return m_Data[column]; }
        PREONMATH_DEVICE inline void setColumn(PrMathSize row, const vec<M, T>& value);

        // Row getter & setter
        PREONMATH_DEVICE inline vec<N, T> row(PrMathSize row) const;
        PREONMATH_DEVICE inline void setRow(PrMathSize row, const vec<N, T>& value);

        // Data access
        PREONMATH_DEVICE T operator()(PrMathSize row, PrMathSize column) const { return element(row, column); }
        PREONMATH_DEVICE T& operator()(PrMathSize row, PrMathSize column) { return element(row, column); }

        PREONMATH_DEVICE T element(PrMathSize row, PrMathSize column) const { return m_Data[column][row]; }
        PREONMATH_DEVICE T& element(PrMathSize row, PrMathSize column) { return m_Data[column][row]; }

        // Filling
        PREONMATH_DEVICE inline void fill(T value);
        PREONMATH_DEVICE inline void setToZero() { fill(scalar_zero<T>()); }
        PREONMATH_DEVICE static matrix<M, N, T> zero()
        {
            matrix<M, N, T> m;
            m.setToZero();
            return m;
        }

        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline void setToIdentity(typename std::enable_if<_M == _N, T>::type* = nullptr);
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE static matrix<M, N, T> identity(typename std::enable_if<_M == _N, T>::type* = nullptr)
        {
            matrix<M, N, T> m;
            m.setToIdentity();
            return m;
        }

        // Transpose
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline void transpose(typename std::enable_if<_M == _N, T>::type* = nullptr);
        PREONMATH_DEVICE inline matrix<N, M, T> transposed() const;

        // Inverse
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline void invert(typename std::enable_if<_M == _N, T>::type* = nullptr)
        {
            *this = inverted();
        }
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline matrix<N, M, T> inverted(typename std::enable_if<_M == _N && _M == 1, T>::type* = nullptr) const;
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline matrix<N, M, T> inverted(typename std::enable_if<_M == _N && _M == 2, T>::type* = nullptr) const;
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline matrix<N, M, T> inverted(typename std::enable_if<_M == _N && _M == 3, T>::type* = nullptr) const;
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline matrix<N, M, T> inverted(typename std::enable_if<_M == _N && _M == 4, T>::type* = nullptr) const;

        // Adjoint and cofactor.
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline matrix<N, M, T> adjugate(typename std::enable_if<_M == _N && _M == 3, T>::type* = nullptr) const;
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline matrix<M, N, T> cofactorMatrix(typename std::enable_if<_M == _N && _M == 3, T>::type* = nullptr) const
        {
            return adjugate().transposed();
        }

        // Reductions
        template<typename E = T>
        PREONMATH_DEVICE inline enable_if_non_simd<E, E> norm(E p = 2) const;

        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline T determinant(typename std::enable_if<_M == _N && _M != 2 && _M != 3, T>::type* = nullptr) const;
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline T determinant(typename std::enable_if<_M == 2 && _N == 2, T>::type* = nullptr) const;
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline T determinant(typename std::enable_if<_M == 3 && _N == 3, T>::type* = nullptr) const;

        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE inline T trace(typename std::enable_if<_M == _N, T>::type* = nullptr) const;
        template<PrMathSize _M = M, PrMathSize _N = N>
        PREONMATH_DEVICE vec<M, T> diagonal(typename std::enable_if<_M == _N, T>::type* = nullptr) const;

        // Operators
        PREONMATH_DEVICE inline matrix<M, N, T>& operator+=(const matrix<M, N, T>& other);
        PREONMATH_DEVICE inline matrix<M, N, T>& operator-=(const matrix<M, N, T>& other);
        PREONMATH_DEVICE inline matrix<M, N, T>& operator*=(const matrix<M, N, T>& other);  // TODO: Make it only useable for M == N, else it makes no sense.
        PREONMATH_DEVICE inline matrix<M, N, T>& operator*=(T factor);
        PREONMATH_DEVICE inline matrix<M, N, T>& operator/=(T divisor)
        {
            *this *= (Simd::scalar_cast<float, T>(1) / divisor);
            return *this;
        }

        // Helpers
        template<typename T_Out, typename T_Vec>
        PREONMATH_DEVICE inline T_Out rowVecProduct(PrMathSize row, const vec<N, T_Vec>& vec) const;

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
        PREONMATH_DEVICE void legacyInitZeroFunc(enable_if_non_simd<E, void*> = nullptr)
        {
#ifdef PREONMATH_DEFAULTINIT
            setToZero();
#endif
        }
        template<typename E>
        PREONMATH_DEVICE void legacyInitZeroFunc(enable_if_simd<E, void*> = nullptr)
        {
        }

        template<typename E>
        PREONMATH_DEVICE void legacyInitIdentityFunc(enable_if_non_simd<E, void*> = nullptr)
        {
#ifdef PREONMATH_DEFAULTINIT
            setToIdentity();
#endif
        }
        template<typename E>
        PREONMATH_DEVICE void legacyInitIdentityFunc(enable_if_simd<E, void*> = nullptr)
        {
        }
    };

    // Overloaded non-member operators

    template<PrMathSize M, PrMathSize N, typename T>
    PREONMATH_DEVICE inline bool operator==(const matrix<M, N, T>& m1, const matrix<M, N, T>& m2);
    template<PrMathSize M, PrMathSize N, typename T>
    PREONMATH_DEVICE inline bool operator!=(const matrix<M, N, T>& m1, const matrix<M, N, T>& m2)
    {
        return !(m1 == m2);
    }

    template<PrMathSize M, PrMathSize N, typename T>
    PREONMATH_DEVICE inline matrix<M, N, T> operator+(matrix<M, N, T> m1, const matrix<M, N, T>& m2)
    {
        return m1 += m2;
    }
    template<PrMathSize M, PrMathSize N, typename T>
    PREONMATH_DEVICE inline matrix<M, N, T> operator-(matrix<M, N, T> m1, const matrix<M, N, T>& m2)
    {
        return m1 -= m2;
    }
    template<PrMathSize M, PrMathSize N, typename T>
    PREONMATH_DEVICE inline matrix<M, N, T> operator-(const matrix<M, N, T>& m2)
    {
        return matrix<M, N, T>(0) -= m2;
    }

    template<PrMathSize M, PrMathSize N, PrMathSize O, typename T>
    PREONMATH_DEVICE inline matrix<M, O, T> operator*(const matrix<M, N, T>& m1, const matrix<N, O, T>& m2);
    template<PrMathSize M, PrMathSize N, typename T>
    PREONMATH_DEVICE inline matrix<M, N, T> operator*(const T factor, matrix<M, N, T> m)
    {
        return m *= factor;
    }
    template<PrMathSize M, PrMathSize N, typename T>
    PREONMATH_DEVICE inline matrix<M, N, T> operator*(matrix<M, N, T> m, const T factor)
    {
        return factor * m;
    }

    template<typename T, PrMathSize M, PrMathSize N, typename T_Mat, typename T_Vec>
    PREONMATH_DEVICE inline vec<M, T> operator*(const matrix<M, N, T_Mat>& m, const vec<N, T_Vec>& v);
    template<PrMathSize M, PrMathSize N, typename T>
    PREONMATH_DEVICE inline vec<M, T> operator*(const matrix<M, N, T>& m, const vec<N, T>& v)
    {
        return operator*<T, M, N, T, T>(m, v);
    }

    template<PrMathSize M, PrMathSize N, typename T>
    PREONMATH_DEVICE inline matrix<M, N, T> operator/(matrix<M, N, T> m, const T factor)
    {
        return m /= factor;
    }

    template<PrMathSize M, PrMathSize N, typename T>
    inline enable_if_non_simd<T, std::stringstream>& operator<<(std::stringstream&, const matrix<M, N, T>&);

    template<PrMathSize M, PrMathSize N, typename T>
    inline enable_if_non_simd<T, std::ostream>& operator<<(std::ostream&, const matrix<M, N, T>&);

#ifdef PREONMATH_QT_INTEGRATION
    template<PrMathSize M, PrMathSize N, typename T>
    QDebug operator<<(QDebug dbg, const matrix<M, N, T>& m)
    {
        dbg.nospace() << m.toUtf8();
        return dbg.space();
    }
#endif  // PREONMATH_QT_INTEGRATION
}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
template<PrMathSize M, PrMathSize N, typename T>
inline bool qFuzzyCompare(const Preon::Math::matrix<M, N, T>& m1, const Preon::Math::matrix<M, N, T>& m2)
{
    for (PrMathSize c = 0; c < N; c++)
        for (PrMathSize r = 0; r < M; r++)
            if (!qFuzzyCompare(m1(r, c), m2(r, c)))
                return false;
    return true;
}
#endif  // PREONMATH_QT_INTEGRATION

#include "matrix.inl"
