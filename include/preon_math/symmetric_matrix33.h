/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"
#include "vec.h"
#include "matrix.h"

#include <sstream>
#include <iostream>

namespace Preon
{
namespace Math
{
    template<typename T>
    class SymmetricMatrix33
    {
    public:
        SymmetricMatrix33()
            : m_Data(0) {}

        SymmetricMatrix33(const T& xx, const T& xy, const T& xz, const T& yy, const T& yz, const T& zz)
            : m_Data(xx, xy, xz, yy, yz, zz) {}

        SymmetricMatrix33(const matrix<3, 3, T>& m)
            : SymmetricMatrix33(m(0, 0), m(0, 1), m(0, 2), m(1, 1), m(1, 2), m(2, 2)) {}

        static SymmetricMatrix33 identity() { return SymmetricMatrix33(1.0, 0.0, 0.0, 1.0, 0.0, 1.0); }

        static SymmetricMatrix33 fromSumWithTranspose(const matrix<3, 3, T>& m)
        {
            return SymmetricMatrix33(m(0, 0) + m(0, 0), m(0, 1) + m(1, 0), m(0, 2) + m(2, 0), m(1, 1) + m(1, 1), m(1, 2) + m(2, 1), m(2, 2) + m(2, 2));
        }

        // Data access
        T operator[](size_t i) const { return m_Data[i]; }
        T& operator[](size_t i) { return m_Data[i]; }

        T operator()(size_t row, size_t column) const { return element(row, column); }
        T& operator()(size_t row, size_t column) { return element(row, column); }

        T element(size_t row, size_t column) const { return (row == 0 || column == 0) ? m_Data[row + column] : m_Data[row + column + 1]; }
        T& element(size_t row, size_t column) { return (row == 0 || column == 0) ? m_Data[row + column] : m_Data[row + column + 1]; }

        void addDiagonalMatrix(const vec<3, T>& diagonal)
        {
            this->m_Data[0] += diagonal[0];
            this->m_Data[3] += diagonal[1];
            this->m_Data[5] += diagonal[2];
        }

        matrix<3, 3, T> toMatrix() const
        {
            matrix<3, 3, T> m;
            m(0, 0) = m_Data[0];
            m(0, 1) = m(1, 0) = m_Data[1];
            m(0, 2) = m(2, 0) = m_Data[2];
            m(1, 1) = m_Data[3];
            m(1, 2) = m(2, 1) = m_Data[4];
            m(2, 2) = m_Data[5];
            return m;
        }

        T trace() const { return m_Data[0] + m_Data[3] + m_Data[5]; }

        template<typename F>
        void operator*=(F factor)
        {
            m_Data *= factor;
        }

        vec<3, T> operator*(const vec<3, T>& v) const
        {
            return vec<3, T>(
                m_Data[0] * v[0] + m_Data[1] * v[1] + m_Data[2] * v[2],
                m_Data[1] * v[0] + m_Data[3] * v[1] + m_Data[4] * v[2],
                m_Data[2] * v[0] + m_Data[4] * v[1] + m_Data[5] * v[2]);
        }

    private:
        vec<6, T> m_Data;
    };

    template<typename T>
    inline bool operator==(const SymmetricMatrix33<T>& m1, const SymmetricMatrix33<T>& m2)
    {
        for (size_t i = 0; i < 6; ++i)
            if (m1[i] != m2[i])
                return false;
        return true;
    }
    template<typename T>
    inline bool operator!=(const SymmetricMatrix33<T>& m1, const SymmetricMatrix33<T>& m2)
    {
        return !(m1 == m2);
    }

    template<typename T, typename F>
    SymmetricMatrix33<T> operator*(F factor, SymmetricMatrix33<T> m)
    {
        m *= factor;
        return m;
    }

    template<typename T>
    inline std::stringstream& operator<<(std::stringstream& stream, const SymmetricMatrix33<T>& symmetricMatrix)
    {
        stream << symmetricMatrix.toMatrix();
        return stream;
    }

    template<size_t M, size_t N, typename T>
    inline std::ostream& operator<<(std::ostream& ostream, const SymmetricMatrix33<T>& symmetricMatrix)
    {
        std::stringstream ss;
        ss << symmetricMatrix;
        return ostream << ss.rdbuf();
    }

    typedef SymmetricMatrix33<float> SymmetricMatrix33f;
    typedef SymmetricMatrix33<double> SymmetricMatrix33d;
}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
Q_DECLARE_TYPEINFO(SymmetricMatrix33f, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(SymmetricMatrix33d, Q_MOVABLE_TYPE);
Q_DECLARE_METATYPE(SymmetricMatrix33f)
Q_DECLARE_METATYPE(SymmetricMatrix33d)
#endif  // PREONMATH_QT_INTEGRATION
