/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "matrix.h"
#include "matrix33.h"
#include "vec_simd.h"
#include "quat.h"

namespace Preon
{
namespace Math
{
    namespace MatrixUtils
    {
        //! Creates a matrix that ---used as matrix-vector multiplication---
        //! expresses the a cross product.
        //! See: https://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication
        template<typename Scalar>
        PREONMATH_DEVICE matrix<3, 3, Scalar> crossProductMatrix(const vec<3, Scalar>& v)
        {
            matrix<3, 3, Scalar> r = matrix<3, 3, Scalar>::zero();
            r(0, 1) = -v(2);
            r(0, 2) = v(1);
            r(1, 0) = v(2);
            r(1, 2) = -v(0);
            r(2, 0) = -v(1);
            r(2, 1) = v(0);
            return r;
        }

        template<typename T, PrMathSize M, PrMathSize N, typename T_Mat, typename T_Vec>
        PREONMATH_DEVICE inline matrix<M, N, T> scaleColumns(const vec<M, T_Vec>& v, const matrix<M, N, T_Mat>& m)
        {
            matrix<M, N, T> out;
            for (PrMathSize c = 0; c < M; c++)
                out.setColumn(c, m.column(c) * v[c]);
            return out;
        }
        template<PrMathSize M, PrMathSize N, typename T>
        PREONMATH_DEVICE inline matrix<M, N, T> scaleColumns(const vec<M, T>& v, const matrix<M, N, T>& m)
        {
            return scaleColumns<T, M, N, T, T>(v, m);
        }

        template<typename Scalar>
        PREONMATH_DEVICE matrix<4, 4, Scalar> makeOrtho(Scalar left, Scalar right, Scalar bottom, Scalar top, Scalar nearPlane, Scalar farPlane)
        {
            // Bail out if the projection volume is zero-sized.
            if (left == right || bottom == top || nearPlane == farPlane)
                return matrix<4, 4, Scalar>::identity();

            // Construct the projection.
            Scalar width = right - left;
            Scalar invheight = top - bottom;
            Scalar clip = farPlane - nearPlane;

            matrix<4, 4, Scalar> m;  // Initialization is not important as everything is overwritten.
            m(0, 0) = 2.0 / width;
            m(0, 1) = 0.0;
            m(0, 2) = 0.0;
            m(0, 3) = -(left + right) / width;
            m(1, 0) = 0.0;
            m(1, 1) = 2.0 / invheight;
            m(1, 2) = 0.0;
            m(1, 3) = -(top + bottom) / invheight;
            m(2, 0) = 0.0;
            m(2, 1) = 0.0;
            m(2, 2) = -2.0 / clip;
            m(2, 3) = -(nearPlane + farPlane) / clip;
            m(3, 0) = 0.0;
            m(3, 1) = 0.0;
            m(3, 2) = 0.0;
            m(3, 3) = 1.0;

            // Apply the projection.
            return m;
        }

        template<typename Scalar>
        PREONMATH_DEVICE matrix<4, 4, Scalar> makeFrustum(Scalar left, Scalar right, Scalar bottom, Scalar top, Scalar nearPlane, Scalar farPlane)
        {
            // Bail out if the projection volume is zero-sized.
            if (left == right || bottom == top || nearPlane == farPlane)
                return matrix<4, 4, Scalar>::identity();

            // Construct the projection.
            matrix<4, 4, Scalar> m;  // Initialization is not important as everything is overwritten.
            Scalar width = right - left;
            Scalar invheight = top - bottom;
            Scalar clip = farPlane - nearPlane;
            m(0, 0) = 2.0 * nearPlane / width;
            m(0, 1) = 0.0;
            m(0, 2) = (left + right) / width;
            m(0, 3) = 0.0;
            m(1, 0) = 0.0;
            m(1, 1) = 2.0 * nearPlane / invheight;
            m(1, 2) = (top + bottom) / invheight;
            m(1, 3) = 0.0;
            m(2, 0) = 0.0;
            m(2, 1) = 0.0;
            m(2, 2) = -(nearPlane + farPlane) / clip;
            m(2, 3) = -2.0 * nearPlane * farPlane / clip;
            m(3, 0) = 0.0;
            m(3, 1) = 0.0;
            m(3, 2) = -1.0;
            m(3, 3) = 0.0;

            // Apply the projection.
            return m;
        }

        template<typename Scalar>
        PREONMATH_DEVICE void translate(matrix<4, 4, Scalar>* pMatrix, const vec<3, Scalar>& v)
        {
            vec<4, Scalar> v4(v.x(), v.y(), v.z(), 0.0);
            matrix<4, 4, Scalar>& m = *pMatrix;
            m(0, 3) += m.row(0) * v4;
            m(1, 3) += m.row(1) * v4;
            m(2, 3) += m.row(2) * v4;
            m(3, 3) += m.row(3) * v4;
        }

        template<typename Scalar>
        PREONMATH_DEVICE matrix<4, 4, Scalar> makeTranslation(const vec<3, Scalar>& v)
        {
            matrix<4, 4, Scalar> m = matrix<4, 4, Scalar>::identity();
            translate(&m, v);
            return m;
        }

        template<typename Scalar>
        PREONMATH_DEVICE matrix<4, 4, Scalar> makeRotation(const quat<Scalar>& quaternion)
        {
            // See http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q54
            matrix<4, 4, Scalar> m = matrix<4, 4, Scalar>::zero();
            Scalar xx = quaternion.x() * quaternion.x();
            Scalar xy = quaternion.x() * quaternion.y();
            Scalar xz = quaternion.x() * quaternion.z();
            Scalar xw = quaternion.x() * quaternion.scalar();
            Scalar yy = quaternion.y() * quaternion.y();
            Scalar yz = quaternion.y() * quaternion.z();
            Scalar yw = quaternion.y() * quaternion.scalar();
            Scalar zz = quaternion.z() * quaternion.z();
            Scalar zw = quaternion.z() * quaternion.scalar();

            m(0, 0) = 1.0 - 2 * (yy + zz);
            m(0, 1) = 2 * (xy - zw);
            m(0, 2) = 2 * (xz + yw);
            m(1, 0) = 2 * (xy + zw);
            m(1, 1) = 1.0 - 2 * (xx + zz);
            m(1, 2) = 2 * (yz - xw);
            m(2, 0) = 2 * (xz - yw);
            m(2, 1) = 2 * (yz + xw);
            m(2, 2) = 1.0 - 2 * (xx + yy);
            m(3, 3) = 1.0;

            return m;
        }

        template<typename Scalar>
        void rotate(matrix<4, 4, Scalar>* pMatrix, const quat<Scalar>& quaternion)
        {
            *pMatrix *= makeRotation(quaternion);
        }

        template<typename Scalar>
        PREONMATH_DEVICE void rotate(matrix<4, 4, Scalar>* pMatrix, Scalar angle, Scalar x, Scalar y, Scalar z)
        {
            if (angle == 0.0)
                return;

            //        quat<Scalar> q = quat<Scalar>::fromAxisAndAngle(x, y, z, angle);
            //        rotate(pMatrix, q);
            //        return;

            matrix<4, 4, Scalar> m = matrix<4, 4, Scalar>::zero();
            Scalar c, s, ic;
            if (angle == 90.0 || angle == -270.0)
            {
                s = 1.0;
                c = 0.0;
            }
            else if (angle == -90.0 || angle == 270.0)
            {
                s = -1.0;
                c = 0.0;
            }
            else if (angle == 180.0 || angle == -180.0)
            {
                s = 0.0;
                c = -1.0;
            }
            else
            {
                // Angle in degree to radians.
                Scalar a = DegToRad(angle);
                c = std::cos(a);
                s = std::sin(a);
            }
            bool quick = false;
            if (x == 0.0)
            {
                if (y == 0.0)
                {
                    if (z != 0.0)
                    {
                        // Rotate around the Z axis.
                        m.setToIdentity();
                        m(0, 0) = c;
                        m(1, 1) = c;
                        if (z < 0.0)
                        {
                            m(0, 1) = s;
                            m(1, 0) = -s;
                        }
                        else
                        {
                            m(0, 1) = -s;
                            m(1, 0) = s;
                        }
                        quick = true;
                    }
                }
                else if (z == 0.0f)
                {
                    // Rotate around the Y axis.
                    m.setToIdentity();
                    m(0, 0) = c;
                    m(2, 2) = c;
                    if (y < 0.0)
                    {
                        m(0, 2) = -s;
                        m(2, 0) = s;
                    }
                    else
                    {
                        m(0, 2) = s;
                        m(2, 0) = -s;
                    }
                    quick = true;
                }
            }
            else if (y == 0.0 && z == 0.0)
            {
                // Rotate around the X axis.
                m.setToIdentity();
                m(1, 1) = c;
                m(2, 2) = c;
                if (x < 0.0)
                {
                    m(1, 2) = s;
                    m(2, 1) = -s;
                }
                else
                {
                    m(1, 2) = -s;
                    m(2, 1) = s;
                }
                quick = true;
            }
            if (!quick)
            {
                Scalar len = x * x + y * y + z * z;
                if (!IsZero(len - 1.0) && !IsZero(len))
                {
                    len = std::sqrt(len);
                    x /= len;
                    y /= len;
                    z /= len;
                }
                ic = 1.0 - c;
                m(0, 0) = x * x * ic + c;
                m(0, 1) = x * y * ic - z * s;
                m(0, 2) = x * z * ic + y * s;
                m(0, 3) = 0.0;
                m(1, 0) = y * x * ic + z * s;
                m(1, 1) = y * y * ic + c;
                m(1, 2) = y * z * ic - x * s;
                m(1, 3) = 0.0;
                m(2, 0) = x * z * ic - y * s;
                m(2, 1) = y * z * ic + x * s;
                m(2, 2) = z * z * ic + c;
                m(2, 3) = 0.0;
                m(3, 0) = 0.0;
                m(3, 1) = 0.0;
                m(3, 2) = 0.0;
                m(3, 3) = 1.0;
            }
            *pMatrix *= m;
        }

        template<typename Scalar>
        PREONMATH_DEVICE void rotate(matrix<4, 4, Scalar>* pMatrix, Scalar angle, const vec<3, Scalar>& vector)
        {
            rotate(pMatrix, angle, vector.x(), vector.y(), vector.z());
        }

        template<typename Scalar>
        PREONMATH_DEVICE matrix<4, 4, Scalar> makeRotation(Scalar angle, const vec<3, Scalar>& vector)
        {
            matrix<4, 4, Scalar> m = matrix<4, 4, Scalar>::identity();
            rotate(&m, angle, vector.x(), vector.y(), vector.z());
            return m;
        }

        template<typename Scalar>
        PREONMATH_DEVICE void scale(matrix<4, 4, Scalar>* pMatrix, const vec<3, Scalar>& v)
        {
            (*pMatrix)(0, 0) *= v.x();
            (*pMatrix)(1, 0) *= v.x();
            (*pMatrix)(2, 0) *= v.x();
            (*pMatrix)(3, 0) *= v.x();
            (*pMatrix)(0, 1) *= v.y();
            (*pMatrix)(1, 1) *= v.y();
            (*pMatrix)(2, 1) *= v.y();
            (*pMatrix)(3, 1) *= v.y();
            (*pMatrix)(0, 2) *= v.z();
            (*pMatrix)(1, 2) *= v.z();
            (*pMatrix)(2, 2) *= v.z();
            (*pMatrix)(3, 2) *= v.z();
        }

        template<typename Scalar>
        PREONMATH_DEVICE matrix<4, 4, Scalar> makeScale(const vec<3, Scalar>& v)
        {
            matrix<4, 4, Scalar> m = matrix<4, 4, Scalar>::identity();
            scale(&m, v);
            return m;
        }

        template<typename Scalar>
        PREONMATH_DEVICE matrix<4, 4, Scalar> worldTransform(const vec<3, Scalar>& position, const quat<Scalar>& orientation, const vec<3, Scalar>& scaleVector)
        {
            // Build the world transformation.
            // First position.
            matrix<4, 4, Scalar> m = makeTranslation(position);
            // Orientation.
            rotate(&m, orientation);
            // Scale.
            scale(&m, scaleVector);
            return m;
        }

        template<typename Scalar>
        matrix<4, 4, Scalar> PREONMATH_DEVICE worldTransform(const vec<3, Scalar>& position, const quat<Scalar>& orientation, const vec<3, Scalar>& scaleVector, const vec<3, Scalar>& pivot, bool applyScaleOnPivot)
        {
            // Build the world transformation.
            matrix<4, 4, Scalar> m = makeTranslation(position);  // Position first.
            rotate(&m, orientation);  // Orientation.
            // User decided on whether to apply the scaling w.r.t. the pivot or the origin.
            if (applyScaleOnPivot)
            {
                scale(&m, scaleVector);  // Scale.
                translate(&m, pivot);  // Pivot
            }
            else
            {
                translate(&m, pivot);  // Pivot
                scale(&m, scaleVector);  // Scale.
            }
            return m;
        }

        template<typename Scalar>
        PREONMATH_DEVICE matrix<3, 3, Scalar> fromQuaternion(const quat<Scalar>& quaternion)
        {
            // See http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q54
            Scalar xx = quaternion.x() * quaternion.x();
            Scalar xy = quaternion.x() * quaternion.y();
            Scalar xz = quaternion.x() * quaternion.z();
            Scalar xw = quaternion.x() * quaternion.scalar();
            Scalar yy = quaternion.y() * quaternion.y();
            Scalar yz = quaternion.y() * quaternion.z();
            Scalar yw = quaternion.y() * quaternion.scalar();
            Scalar zz = quaternion.z() * quaternion.z();
            Scalar zw = quaternion.z() * quaternion.scalar();

            matrix<3, 3, Scalar> m;
            m(0, 0) = 1 - 2 * (yy + zz);
            m(0, 1) = 2 * (xy + zw);
            m(0, 2) = 2 * (xz - yw);
            m(1, 0) = 2 * (xy - zw);
            m(1, 1) = 1 - 2 * (xx + zz);
            m(1, 2) = 2 * (yz + xw);
            m(2, 0) = 2 * (xz + yw);
            m(2, 1) = 2 * (yz - xw);
            m(2, 2) = 1 - 2 * (xx + yy);
            // Note: This method was wrong before and produced the transpose of
            // the matrix we actually wanted. The above code is now correct.
            // However, to get the same behavior as before, we transpose the matrix
            // here.
            m.transpose();

            return m;
        }

#ifdef PREONMATH_ENABLE_SIMD
        template<typename Scalar>
        PREONMATH_DEVICE enable_if_non_simd<Scalar, bool> computeMaxEigenValue(const matrix<3, 3, Scalar>& m, Scalar& eigenMax, vec<3, Scalar>& eigenVector)
        {
            eigenVector = m.diagonal();
            eigenMax = eigenVector.length();
            Scalar prevLength = eigenMax;
            const Scalar epsilon = static_cast<Scalar>(0.001);
            for (int i = 0; i < 20; i++)
            {
                eigenVector = m * eigenVector;
                eigenMax = eigenVector.length();
                if (IsZero(eigenMax))
                    return false;
                eigenVector /= eigenMax;
                if (std::abs(prevLength - eigenMax) < epsilon)
                    break;
                prevLength = eigenMax;
            }
            return true;
        }
        PREONMATH_DEVICE inline void computeMaxEigenValue(const matrix<3, 3, float_simd>& m, float_simd& eigenMax, vec_simd<3>& eigenVector)
        {
            eigenVector = m.diagonal();
            eigenMax = eigenVector.length();
            static const float_simd zero = Simd::zero<float>();
            for (int i = 0; i < 20; i++)
            {
                eigenVector = m * eigenVector;
                eigenMax = eigenVector.length();
                eigenVector = Simd::and_mask(eigenVector / eigenMax, eigenMax > zero);
            }
        }
        PREONMATH_DEVICE inline void computeMaxEigenValue_Adaptive(const matrix<3, 3, float_simd>& m, float_simd& eigenMax, vec_simd<3>& eigenVector)
        {
            eigenVector = m.diagonal();
            eigenMax = eigenVector.length();
            float_simd prevLength = eigenMax;
            static const float_simd epsilon = Simd::make<float>(0.001f);
            static const float_simd zero = Simd::zero<float>();
            for (int i = 0; i < 20; i++)
            {
                eigenVector = m * eigenVector;
                eigenMax = eigenVector.length();
                if (Simd::testAllZero(eigenMax))
                    return;
                eigenVector = Simd::and_mask(eigenVector / eigenMax, eigenMax > zero);
                if (Simd::testAllZero(Simd::abs(prevLength - eigenMax) >= epsilon))
                    break;
                prevLength = eigenMax;
            }
        }
#endif

        template<typename Scalar>
        PREONMATH_DEVICE void jacobiRotate(matrix<3, 3, Scalar>* pMatrix, matrix<3, 3, Scalar>* R, int p, int q)
        {
            // rotates A through phi in pq-plane to set A.get(p, q) = 0
            // rotation stored in R whose columns are eigenvectors of A
            Scalar d = ((*pMatrix)(p, p) - (*pMatrix)(q, q)) / (2.0 * (*pMatrix)(q, p));
            Scalar t = 1.0 / (std::abs(d) + std::sqrt(d * d + 1.0));
            if (d < 0)
                t = -t;
            Scalar c = 1.0 / std::sqrt(t * t + 1);
            Scalar s = t * c;
            (*pMatrix)(p, p) += t * (*pMatrix)(q, p);
            (*pMatrix)(q, q) -= t * (*pMatrix)(q, p);
            (*pMatrix)(q, p) = 0;
            (*pMatrix)(p, q) = 0;
            // transform A
            for (int k = 0; k < 3; k++)
            {
                if (k != p && k != q)
                {
                    Scalar pk = (*pMatrix)(k, p) = c * (*pMatrix)(p, k) + s * (*pMatrix)(q, k);
                    Scalar qk = (*pMatrix)(k, q) = -s * (*pMatrix)(p, k) + c * (*pMatrix)(q, k);
                    (*pMatrix)(p, k) = pk;
                    (*pMatrix)(q, k) = qk;
                }
            }
            // store rotation in R
            for (int k = 0; k < 3; k++)
            {
                Scalar pk = c * (*R)(p, k) + s * (*R)(q, k);
                Scalar qk = -s * (*R)(p, k) + c * (*R)(q, k);
                (*R)(p, k) = pk;
                (*R)(q, k) = qk;
            }
        }

        template<typename Scalar>
        PREONMATH_DEVICE void eigenDecomposition(matrix<3, 3, Scalar>* pMatrix, matrix<3, 3, Scalar>* R, float jacobiEpsilon = 1e-15f, int jacobiIterations = 100)
        {
            // only for symmetric matrices!
            // A = R A' R^T, where A' is diagonal and R orthonormal
            R->setToIdentity();  // unit matrix
            int iter = 0;
            while (iter < jacobiIterations)  // 3 off diagonal elements
            {
                // find off diagonal element with maximum modulus
                int p = 0;
                int q = 1;
                Scalar max = std::abs((*pMatrix)(1, 0));
                Scalar a = std::abs((*pMatrix)(2, 0));
                if (a > max)
                {
                    p = 0;
                    q = 2;
                    max = a;
                }
                a = std::abs((*pMatrix)(2, 1));
                if (a > max)
                {
                    p = 1;
                    q = 2;
                    max = a;
                }
                // all small enough -> done
                if (max < jacobiEpsilon)
                    break;
                // rotate matrix with respect to that element
                jacobiRotate(pMatrix, R, p, q);
                ++iter;
            }
        }

        template<PrMathSize N, PrMathSize M, typename T>
        PREONMATH_DEVICE bool isZero(const matrix<N, M, T>& m)
        {
            for (PrMathSize c = 0; c < N; c++)
                for (PrMathSize r = 0; r < M; r++)
                    if (!IsZero(m(r, c)))
                        return false;
            return true;
        }

        template<typename Scalar>
        PREONMATH_DEVICE vec<3, Scalar> mapVector(const matrix<4, 4, Scalar>& m, const vec<3, Scalar>& v)
        {
            vec<4, Scalar> v4(v[0], v[1], v[2], 0.0);
            v4 = m * v4;
            return vec<3, Scalar>(v4[0], v4[1], v4[2]);
        }

        //! Inverts the MSub x MSub upper-left sub matrix of the given M x N matrix using doubles and returns the result. Returns a MSub x MSub matrix (by default MxM).
        template<PrMathSize M, PrMathSize N, typename Scalar, PrMathSize MSub = M>
        PREONMATH_DEVICE typename std::enable_if<MSub <= M && MSub <= N, matrix<MSub, MSub, Scalar>>::type invertHighPrecision(const matrix<M, N, Scalar>& m)
        {
            auto doubleInverted = matrix<MSub, MSub, double>([&](PrMathSize r, PrMathSize c) { return static_cast<double>(m(r, c)); }).inverted();
            return matrix<MSub, MSub, Scalar>([&](PrMathSize r, PrMathSize c) { return static_cast<Scalar>(doubleInverted(r, c)); });
        }

        //! Computes the 3x3 normal matrix given a 4x4 world matrix.
        template<typename Scalar>
        PREONMATH_DEVICE matrix<3, 3, Scalar> makeNormalMatrix(const matrix<4, 4, Scalar>& worldTransform)
        {
            // We need the inverted transpose to transform the normals.
            // Reference: http://gamedev.stackexchange.com/questions/44511/how-to-transform-mesh-components
            // Always perform matrix inversion using doubles.
            return invertHighPrecision<4, 4, Scalar, 3>(worldTransform).transposed();
        }

        //! Extracts the 3x3 rotation matrix given a 4x4 transform matrix.
        //! Attention: The result is only valid if no shear or stress is part of the transform, i.e., the matrix may contain translation, rotation and scale.
        template<typename Scalar>
        PREONMATH_DEVICE matrix<3, 3, Scalar> extractRotationMatrix(const matrix<4, 4, Scalar>& transform)
        {
            // Compute scale.
            vec<3, Scalar> scale(transform.column(0).length(), transform.column(1).length(), transform.column(2).length());
            matrix<3, 3, Scalar> rot;
            for (PrMathSize i = 0; i < 3; i++)
                for (PrMathSize j = 0; j < 3; j++)
                    rot(i, j) = IsZero(scale[i]) ? transform(i, j) : transform(i, j) / scale[i];
            return rot;
        }
    }  // namespace MatrixUtils

    template<typename OutV, typename InM, typename InV>
    PREONMATH_DEVICE inline vec<3, OutV> operator*(const matrix<4, 4, InM>& m, const vec<3, InV>& v)
    {
        vec<4, InV> v4(v[0], v[1], v[2], 1.0);
        v4 = operator*<InV>(m, v4);
        // THROW_EXCEPTION(!AreEqual(v4[3], static_cast<T>(1.0)), "matrix44 * vec3 results in non-normalized vec3.");
        return (1 / v4[3]) * vec<3, OutV>(v4[0], v4[1], v4[2]);
    }

    template<typename T>
    PREONMATH_DEVICE inline vec<3, T> operator*(const matrix<4, 4, T>& m, const vec<3, T>& v)
    {
        return operator*<T, T, T>(m, v);
    }

}  // namespace Math
}  // namespace Preon
