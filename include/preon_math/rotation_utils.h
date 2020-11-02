/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "vec.h"
#include "vec3.h"
#include "quat.h"
#include "matrix.h"
#include "matrix33.h"
#include "euler.h"
#include "scalar_simd.h"

#include <cmath>

namespace Preon
{
namespace Math
{
    namespace RotationUtils
    {
        template<typename T>
        matrix<3, 3, T> matrixFromQuaternions(const quat<T>& q)
        {
            double qw = q.scalar();
            double qx = q.x();
            double qy = q.y();
            double qz = q.z();

            // normalization
            const double n = 1.0 / sqrt(qx * qx + qy * qy + qz * qz + qw * qw);
            qx *= n;
            qy *= n;
            qz *= n;
            qw *= n;

            return matrix<3, 3, T>(
                (T)(1 - 2 * qy * qy - 2 * qz * qz),
                (T)(2 * qx * qy - 2 * qz * qw),
                (T)(2 * qx * qz + 2 * qy * qw),
                (T)(2 * qx * qy + 2 * qz * qw),
                (T)(1 - 2 * qx * qx - 2 * qz * qz),
                (T)(2 * qy * qz - 2 * qx * qw),
                (T)(2 * qx * qz - 2 * qy * qw),
                (T)(2 * qy * qz + 2 * qx * qw),
                (T)(1 - 2 * qx * qx - 2 * qy * qy));
        }

        //! Computes the rotation matrix from euler angle e (assumes x-convention: rotations around x,y and z with angles Phi, Theta, Psi
        //! The conversion is explained here http://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotation_matrix_.E2.86.94_Euler_angles
        template<typename T>
        matrix33d matrixFromEulerAngle(const euler<T>& e)
        {
            // helper variables - Heading = First rotation axis, Attitude = Second, Bank = Last
            double phiRad = DegToRad<double>(e.phi());
            double thetaRad = DegToRad<double>(e.theta());
            double psiRad = DegToRad<double>(e.psi());

            double cosPhi(cos(phiRad));
            double sinPhi(sin(phiRad));
            double cosTheta(cos(thetaRad));
            double sinTheta(sin(thetaRad));
            double cosPsi(cos(psiRad));
            double sinPsi(sin(psiRad));

            // Computation explained here: http://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToMatrix/
            // Phi is around x-axis
            matrix33d rotPhi(1.0, 0.0, 0.0, 0.0, cosPhi, -sinPhi, 0.0, sinPhi, cosPhi);
            // Theta is around y-axis
            matrix33d rotTheta(cosTheta, 0.0, sinTheta, 0.0, 1.0, 0.0, -sinTheta, 0.0, cosTheta);
            // Psi is around z-axis
            matrix33d rotPsi(cosPsi, -sinPsi, 0.0, sinPsi, cosPsi, 0.0, 0.0, 0.0, 1.0);

            return rotPsi * rotTheta * rotPhi;
        }

        //! Returns euler Angles from matrix in ZYX convention
        template<typename T>
        euler<T> eulerAnglesFromMatrix(const matrix<3, 3, T>& m)
        {
            // Extract sin(pitch) from m23.
            T heading, pitch, bank;
            T sp = -m(2, 0);

            // Check for Gimbal lock
            if (std::abs(sp) > 9.99999)
            {
                // Looking straight up or down
                pitch = FPiHalf * sp;

                // Compute heading, slam bank to zero
                heading = atan2(-m(0, 1), m(1, 1));
                bank = 0.0f;
            }
            else
            {
                // Compute angles.  We don't have to use the "safe" asin function because we already
                // checked for range errors when checking for Gimbel lock
                heading = RadToDeg(std::atan2(m(1, 0), m(0, 0)));
                pitch = RadToDeg(std::asin(sp));
                bank = RadToDeg(std::atan2(m(2, 1), m(2, 2)));
            }

            return euler<T>(bank, pitch, heading);
        }

        template<typename T>
        euler<T> eulerAnglesFromQuaternion(const quat<T>& q)
        {
            // First convert it to a matrix and then to euler angles.
            matrix33d m = matrixFromQuaternions(static_cast<quatd>(q));

            // Somehow this has to be differently than our own eulerAnglesFromMatrix
            // method. See https://github.com/erich666/GraphicsGems
            double cy = sqrt(m(0, 0) * m(0, 0) + m(1, 0) * m(1, 0));
            if (cy > FEps)
            {
                return euler<T>((T)RadToDeg(std::atan2(m(2, 1), m(2, 2))), (T)RadToDeg(std::atan2(-m(2, 0), cy)), (T)RadToDeg(std::atan2(m(1, 0), m(0, 0))));
            }
            // Else.
            return euler<T>((T)RadToDeg(std::atan2(-m(1, 2), m(1, 1))), (T)RadToDeg(std::atan2(-m(2, 0), cy)), 0);
        }

        template<typename T>
        quat<T> quaternionsFromEulerAngles(const euler<T>& e)
        {
            quat<T> quatX, quatY, quatZ;
            quatX = quatX.fromAxisAndAngle(vec<3, T>(1, 0, 0), e.phi());
            quatY = quatY.fromAxisAndAngle(vec<3, T>(0, 1, 0), e.theta());
            quatZ = quatZ.fromAxisAndAngle(vec<3, T>(0, 0, 1), e.psi());
            return quatZ * quatY * quatX;
        }

        template<typename T>
        quat<T> quaternionsFromMatrix(const matrix<3, 3, T>& m)
        {
            return quaternionsFromEulerAngles(eulerAnglesFromMatrix(m));
        }

        inline vec3d angularVelocityViaRotationMatrix(const matrix33d& Rot1, const matrix33d& Rot2, qreal dt)
        {
            matrix33d omegaTilde = 1.0 / dt * (Rot2 - Rot1) * Rot1.transposed();
            return vec3d(omegaTilde(2, 1), omegaTilde(0, 2), omegaTilde(1, 0));
        }

        inline vec3d angularVelocityViaEuler(const eulerd& e1, const eulerd& e2, qreal dt)
        {
            // If Euler angles are the same, angular velocity must be zero
            if (e1 == e2)
                return vec3d::zero();
            return angularVelocityViaRotationMatrix(matrixFromEulerAngle(e1), matrixFromEulerAngle(e2), dt);
        }
    }  // namespace RotationUtils
}  // namespace Math
}  // namespace Preon
