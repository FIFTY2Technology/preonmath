/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "euler_fwd.h"
#include "vec.h"

#ifdef PREONMATH_QT_INTEGRATION
    #include <QDebug>
    #include <QDataStream>
#endif  // PREONMATH_QT_INTEGRATION

namespace Preon
{
namespace Math
{
    template<typename T>
    class euler
    {
    public:
        euler()
            : euler(T(0), T(0), T(0)) {}
        euler(T phi, T theta, T psi)
            : m_Data(phi, theta, psi) {}
        explicit euler(const vec<3, T>& vector)
            : m_Data(vector) {}
        template<typename Real2>
        explicit euler(const euler<Real2>& e)
            : m_Data((T)e.phi(), (T)e.theta(), (T)e.psi())
        {
        }
        template<typename Lambda>
        euler(const Lambda& f)
            : m_Data(f)
        {
        }

        T phi() const { return m_Data[0]; }
        T& phi() { return m_Data[0]; }
        T theta() const { return m_Data[1]; }
        T& theta() { return m_Data[1]; }
        T psi() const { return m_Data[2]; }
        T& psi() { return m_Data[2]; }

        void setPhi(T phi) { m_Data[0] = phi; }
        void setTheta(T theta) { m_Data[1] = theta; }
        void setPsi(T psi) { m_Data[2] = psi; }
        void set(T phi, T theta, T psi) { m_Data.set(phi, theta, psi); }
        void set(const vec<3, T>& vec) { m_Data = vec; }

        friend bool operator==(const euler& v1, const euler& v2) { return v1.m_Data == v2.m_Data; }
        friend bool operator!=(const euler& v1, const euler& v2) { return v1.m_Data != v2.m_Data; }
        friend bool operator<(const euler& v1, const euler& v2) { return v1.m_Data < v2.m_Data; }
        friend const euler operator+(const euler& v1, const euler& v2) { return euler(v1.m_Data + v2.m_Data); }

        T& operator[](PrMathSize i) { return m_Data[i]; }
        const T& operator[](PrMathSize i) const { return m_Data[i]; }
        const std::array<T, 3>& data() const { return m_Data.data(); }
        auto value() const { return euler<decltype(std::declval<T>().value())>(m_Data.value()); }

        template<typename _T>
        friend std::stringstream& operator<<(std::stringstream& stream, const euler<_T>&);

    private:
        vec<3, T> m_Data;
    };

    template<typename T>
    std::stringstream& operator<<(std::stringstream& stream, const euler<T>& val)
    {
        stream << val.m_Data;
        return stream;
    }

}  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
// Debugging streams.
template<typename T>
inline QDebug operator<<(QDebug dbg, const euler<T>& ori)
{
    dbg.nospace() << "eulerf(" << ori.phi() << ", " << ori.theta() << ", " << ori.psi() << ')';
    return dbg.space();
}
template<typename T>
inline QDataStream& operator<<(QDataStream& stream, const euler<T>& ori)
{
    stream << T(ori.phi()) << T(ori.theta()) << T(ori.psi());
    return stream;
}

template<typename T>
inline QDataStream& operator>>(QDataStream& stream, euler<T>& ori)
{
    T x, y, z;
    stream >> x;
    stream >> y;
    stream >> z;
    ori.setPhi(T(x));
    ori.setTheta(T(y));
    ori.setPsi(T(z));
    return stream;
}

template<typename T>
inline bool qFuzzyCompare(const euler<T>& e1, const euler<T>& e2)
{
    return qFuzzyCompare(e1.m_Data, e2.m_Data);
}

#endif  // PREONMATH_QT_INTEGRATION
