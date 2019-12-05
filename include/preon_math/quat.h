/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_QUAT_H
#define PREONMATH_QUAT_H

#include "compile_helper.h"

#include "vec.h"
#include "math_utils.h"

#ifdef PREONMATH_QT_INTEGRATION
#include <QDataStream>
#include <QDebug>
#endif  // PREONMATH_QT_INTEGRATION

namespace Preon
{
    namespace Math
    {
        template<typename Real>
        class quat
        {
        public:
            quat() : m_W(1), m_X(0), m_Y(0), m_Z(0) {}
            quat(Real scalar, Real xpos, Real ypos, Real zpos) : m_W(scalar), m_X(xpos), m_Y(ypos), m_Z(zpos) {}
            quat(Real scalar, const vec<3, Real>& vector) : m_W(scalar), m_X(vector.x()), m_Y(vector.y()), m_Z(vector.z()) {}
            explicit quat(const vec<4, Real>& vector) : m_W(vector.w()), m_X(vector.x()), m_Y(vector.y()), m_Z(vector.z()) {}
            template <typename Real2> explicit quat(const quat<Real2>& q) : m_W((Real)q.scalar()), m_X((Real)q.x()), m_Y((Real)q.y()), m_Z((Real)q.z()) {}

            //! Returns the zero-rotation quaternion.
            static quat<Real> identity() { return quat<Real>(1, 0, 0, 0); }

            vec<3, Real> toVector3() const { return vec<3, Real>(m_X, m_Y, m_Z); }
            void setVector(const vec<3, Real>& vector) { setVector(vector.x(), vector.y(), vector.z()); }
            inline void setVector(Real x, Real y, Real z)
            {
                m_X = x;
                m_Y = y;
                m_Z = z;
            }

            vec<4, Real> toVector4() const { return vec<4, Real>(m_X, m_Y, m_Z, m_W); }

            Real x() const { return m_X; }
            Real y() const { return m_Y; }
            Real z() const { return m_Z; }
            Real scalar() const { return m_W; }

            Real& x() { return m_X; }
            Real& y() { return m_Y; }
            Real& z() { return m_Z; }
            Real& scalar() { return m_W; }

            void setX(Real x) { m_X = x; }
            void setY(Real y) { m_Y = y; }
            void setZ(Real z) { m_Z = z; }
            void setScalar(Real scalar) { m_W = scalar; }

            Real lengthSquared() const { return m_X * m_X + m_Y * m_Y + m_Z * m_Z + m_W * m_W; }
            Real length() const { return std::sqrt(lengthSquared()); }

            void normalize()
            {
                // Need some extra precision if the length is very small.
                Real len = length();
                THROW_EXCEPTION(IsZero(len), "Cannot normalize a zero quaternion!")
                if (AreEqual(len, static_cast<Real>(1)))    // this was in the old code, do we really need / want this check?
                    return;
                *this /= len;
            }
            quat<Real> normalized() const
            {
                quat<Real> out = *this;
                out.normalize();
                return out;
            }

            inline quat<Real> conjugate() const
            {
                return quat<Real>(m_W, -m_X, -m_Y, -m_Z);
            }
            quat<Real> inverse() const
            {
                return conjugate() / length();
            }

            vec<3, Real> rotatedVector(const vec<3, Real>& vector) const
            {
                return (*this * quat<Real>(0, vector) * conjugate()).toVector3();
            }

            //! Returns the quaternion that transformed the direction from into the direction to. Both from and to must be normalized directions.
            static quat<Real> rotationBetween_nothrow(const vec<3, Real>& from, const vec<3, Real>& to, double minCrossProductLength = 0.0001)
            {
                // the only purpose of this clamping is to deal with Realing point rounding errors. We expect that from and to are normalized.
                Real angle = std::acos(MathUtils::clamp(vec<3, Real>::dotProduct(from, to), Real{-1}, Real{1}));
                vec<3, Real> axis = vec<3, Real>::crossProduct(from, to);
                Real length = axis.length();
                if (length > minCrossProductLength)
                {
                    // may occur if vectors point in the same or oposite direction
                    axis /= length;
                }
                else
                {
                    axis = from.orthogonalDirection();
                }
                return quat<Real>::fromAxisAndAngleRad(axis, angle);
            }

            //! Returns the quaternion that transformed the direction from into the direction to. Both from and to must be normalized directions.
            static quat<Real> rotationBetween(const vec<3, Real>& from, const vec<3, Real>& to, double minCrossProductLength = 0.0001)
            {
                THROW_EXCEPTION(!AreEqual(Real(1), from.lengthSquared()), "from is not normalized!")
                THROW_EXCEPTION(!AreEqual(Real(1), to.lengthSquared()), "to is not normalized!")

                return rotationBetween_nothrow(from, to, minCrossProductLength);
            }

            static quat<Real> product(const quat<Real>& q, const quat<Real>& p)
            {
                quat<Real> result;
                //w = q1p1 - q2p2 - q3p3 - q4p4
                result.setScalar(q.scalar()*p.scalar() - q.x()*p.x() - q.y()*p.y() - q.z()*p.z());
                //x = q2p1 - q1p2 - q4p3 - q3p4
                result.setX(q.x()*p.scalar() - q.scalar()*p.x() - q.z()*p.y() - q.y()*p.z());
                //y = q3p1 - q4p2 - q1p3 - q2p4
                result.setY(q.y()*p.scalar() - q.z()*p.x() - q.scalar()*p.y() - q.x()*p.z());
                //Z = q4p1 - q3p2 - q2p3 - q1p4
                result.setZ(q.z()*p.scalar() - q.y()*p.x() - q.x()*p.y() - q.scalar()*p.z());
                return result;
            }

            quat<Real>& operator+=(const quat<Real>& quat)
            {
                m_X += quat.m_X;
                m_Y += quat.m_Y;
                m_Z += quat.m_Z;
                m_W += quat.m_W;
                return *this;
            }
            quat<Real>& operator-=(const quat<Real>& quat)
            {
                m_X -= quat.m_X;
                m_Y -= quat.m_Y;
                m_Z -= quat.m_Z;
                m_W -= quat.m_W;
                return *this;
            }
            quat<Real>& operator*=(Real factor)
            {
                m_X *= factor;
                m_Y *= factor;
                m_Z *= factor;
                m_W *= factor;
                return *this;
            }
            quat<Real>& operator*=(const quat<Real>& quat)
            {
                m_X *= quat.m_X;
                m_Y *= quat.m_Y;
                m_Z *= quat.m_Z;
                m_W *= quat.m_W;
                return *this;
            }
            quat<Real>& operator/=(Real divisor)
            {
                m_X /= divisor;
                m_Y /= divisor;
                m_Z /= divisor;
                m_W /= divisor;
                return *this;
            }

            static quat<Real> fromAxisAndAngle(const vec<3, Real>& axis, Real angle)
            {
                return fromAxisAndAngleRad(axis, DegToRad(angle));
            }

            static quat<Real> fromAxisAndAngleRad(const vec<3, Real>& axis, Real angle)
            {
                // See http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q56
                // We normalize the result just in case the values are close
                // to zero, as suggested in the above FAQ.
                vec<3, Real> axisNormlized = axis.pseudoNormalized();
                Real a = angle / 2;
                Real s = std::sin(a);
                Real c = std::cos(a);
                return quat<Real>(c,
                                  axisNormlized[0] * s,
                                  axisNormlized[1] * s,
                                  axisNormlized[2] * s).normalized();
            }

            void toAxisAndAngle(vec<3, Real>* axis, Real* angle)
            {
                // Normalize just to be sure.
                normalize();
                // Calc the angle.

                *angle = Real{2} * std::acos(MathUtils::clamp(m_W, Real{-1}, Real{1}));   // Note that acos(1 + Eps) is NaN!
                Real s = std::sqrt(std::max(Real{0}, Real{1} - m_W * m_W));
                // Check for division by 0.
                if (IsZero(s))
                {
                    axis->setX(m_X);
                    axis->setY(m_Y);
                    axis->setZ(m_Z);
                }
                else
                {
                    axis->setX(m_X / s);
                    axis->setY(m_Y / s);
                    axis->setZ(m_Z / s);
                }
            }

            static quat<Real> slerp(const quat<Real>& q1, const quat<Real>& q2, Real t)
            {
                // Handle the easy cases first.
                if (t <= 0)
                    return q1;
                else if (t >= 1)
                    return q2;

                // Determine the angle between the two quaternions.
                quat<Real> q2b;
                Real dot;
                dot = q1.m_X * q2.m_X + q1.m_Y * q2.m_Y + q1.m_Z * q2.m_Z + q1.m_W * q2.m_W;
                if (dot >= 0)
                {
                    q2b = q2;
                }
                else
                {
                    q2b = -q2;
                    dot = -dot;
                }

                // Get the scale factors.  If they are too small,
                // then revert to simple linear interpolation.
                Real factor1 = 1 - t;
                Real factor2 = t;
                if ((1 - dot) > 0.0000001)
                {
                    Real angle = Real(std::acos(dot));
                    Real sinOfAngle = Real(std::sin(angle));
                    if (sinOfAngle > 0.0000001)
                    {
                        factor1 = Real(std::sin((1 - t) * angle)) / sinOfAngle;
                        factor2 = Real(std::sin(t * angle)) / sinOfAngle;
                    }
                }

                // Construct the result quaternion.
                return q1 * factor1 + q2b * factor2;
            }

            static quat<Real> nlerp(const quat<Real>& q1, const quat<Real>& q2, Real t)
            {
                // Handle the easy cases first.
                if (t <= 0)
                    return q1;
                else if (t >= 1)
                    return q2;

                // Determine the angle between the two quaternions.
                quat<Real> q2b;
                Real dot;
                dot = q1.m_X * q2.m_X + q1.m_Y * q2.m_Y + q1.m_Z * q2.m_Z + q1.m_W * q2.m_W;
                if (dot >= 0)
                    q2b = q2;
                else
                    q2b = -q2;

                // Perform the linear interpolation.
                return (q1 * (1 - t) + q2b * t).normalized();
            }

        private:
            Real m_W, m_X, m_Y, m_Z;
        };

        // Define float and double versions.
        typedef quat<float> quatf;
        typedef quat<double> quatd;

        template<typename Real>
        inline const quat<Real> operator*(const quat<Real>& q1, const quat<Real>& q2)
        {
            Real ww = (q1.z() + q1.x()) * (q2.x() + q2.y());
            Real yy = (q1.scalar() - q1.y()) * (q2.scalar() + q2.z());
            Real zz = (q1.scalar() + q1.y()) * (q2.scalar() - q2.z());
            Real xx = ww + yy + zz;
            Real qq = 0.5 * (xx + (q1.z() - q1.x()) * (q2.x() - q2.y()));

            Real w = qq - ww + (q1.z() - q1.y()) * (q2.y() - q2.z());
            Real x = qq - xx + (q1.x() + q1.scalar()) * (q2.x() + q2.scalar());
            Real y = qq - yy + (q1.scalar() - q1.x()) * (q2.y() + q2.z());
            Real z = qq - zz + (q1.z() + q1.y()) * (q2.scalar() - q2.x());

            return quat<Real>(w, x, y, z);
        }

        template<typename Real>
        inline bool operator==(const quat<Real>& q1, const quat<Real>& q2)
        {
            return q1.x() == q2.x() && q1.y() == q2.y() && q1.z() == q2.z() && q1.scalar() == q2.scalar();
        }

        template<typename Real>
        inline bool operator!=(const quat<Real>& q1, const quat<Real>& q2)
        {
            return q1.x() != q2.x() || q1.y() != q2.y() || q1.z() != q2.z() || q1.scalar() != q2.scalar();
        }

        template<typename Real>
        inline const quat<Real> operator+(const quat<Real>& q1, const quat<Real>& q2)
        {
            return quat<Real>(q1.scalar() + q2.scalar(), q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z());
        }

        template<typename Real>
        inline const quat<Real> operator-(const quat<Real>& q1, const quat<Real>& q2)
        {
            return quat<Real>(q1.scalar() - q2.scalar(), q1.x() - q2.x(), q1.y() - q2.y(), q1.z() - q2.z());
        }

        template<typename Real>
        inline const quat<Real> operator*(Real factor, const quat<Real>& q)
        {
            return quat<Real>(q.scalar() * factor, q.x() * factor, q.y() * factor, q.z() * factor);
        }

        template<typename Real>
        inline const quat<Real> operator*(const quat<Real>& q, Real factor)
        {
            return factor * q;
        }

        template<typename Real>
        inline const quat<Real> operator-(const quat<Real>& q)
        {
            return quat<Real>(-q.scalar(), -q.x(), -q.y(), -q.z());
        }

        template<typename Real>
        inline const quat<Real> operator/(const quat<Real>& q, Real divisor)
        {
            return quat<Real>(q.scalar() / divisor, q.x() / divisor, q.y() / divisor, q.z() / divisor);
        }

        template <typename Real>
        std::stringstream& operator<<(std::stringstream& stream, const quat<Real>& q)
        {
            stream << "quat(scalar:" << q.scalar()
                   << ", vector:(" << q.x() << ", " << q.y() << ", " << q.z() << "))";
            return stream;
        }

        template <typename Real>
        inline std::ostream& operator<<(std::ostream& stream,const quat<Real>& q)
        {
            std::stringstream ss;
            ss << q;
            stream << ss.str();
            return stream;
        }

        #ifdef PREONMATH_QT_INTEGRATION
        template<typename Real>
        QDebug operator<<(QDebug dbg, const quat<Real>& q)
        {
            std::stringstream ss;
            ss << q;
            dbg.nospace() << QString::fromStdString(ss.str());
            return dbg.space();
        }

        template<typename Real>
        QDataStream& operator<<(QDataStream& stream, const quat<Real>& quat)
        {
            stream << Real(quat.scalar()) << Real(quat.x())
                   << Real(quat.y()) << Real(quat.z());
            return stream;
        }

        template<typename Real>
        QDataStream& operator>>(QDataStream& stream, quat<Real>& quat)
        {
            Real scalar, x, y, z;
            stream >> scalar;
            stream >> x;
            stream >> y;
            stream >> z;
            quat.setScalar(scalar);
            quat.setX(x);
            quat.setY(y);
            quat.setZ(z);
            return stream;
        }
        #endif  // PREONMATH_QT_INTEGRATION
    }  // namespace Math
}  // namespace Preon

#ifdef PREONMATH_QT_INTEGRATION
template<typename Real>
inline bool qFuzzyCompare(const Preon::Math::quat<Real>& q1, const Preon::Math::quat<Real>& q2)
{
    return qFuzzyCompare(q1.x(), q2.x()) &&
           qFuzzyCompare(q1.y(), q2.y()) &&
           qFuzzyCompare(q1.z(), q2.z()) &&
           qFuzzyCompare(q1.scalar(), q2.scalar());
}

Q_DECLARE_TYPEINFO(quatf, Q_MOVABLE_TYPE);
Q_DECLARE_METATYPE(quatf)
Q_DECLARE_TYPEINFO(quatd, Q_MOVABLE_TYPE);
#endif  // PREONMATH_QT_INTEGRATION

#endif  // PREONMATH_QUAT_H
