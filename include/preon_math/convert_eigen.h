/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "matrix.h"
#include "quat.h"

#include <Eigen/Geometry>

namespace Preon
{
namespace Math
{
    template<typename scalar, PrMathSize M, PrMathSize N>
    Eigen::Matrix<scalar, static_cast<int>(M), static_cast<int>(N)> toEigen(const matrix<M, N, scalar>& m)
    {
        Eigen::Matrix<scalar, M, N> em;
        for (PrMathSize row = 0; row < M; row++)
            for (PrMathSize column = 0; column < N; column++)
                em(row, column) = m(row, column);
        return em;
    }

    // Order needs to match the Eigen matrix so that the template arguments can be automatically deduced.
    template<typename scalar, int M, int N>
    matrix<static_cast<PrMathSize>(M), static_cast<PrMathSize>(N), scalar> fromEigen(const Eigen::Matrix<scalar, M, N>& em)
    {
        matrix<M, N, scalar> m;
        for (int row = 0; row < M; row++)
            for (int column = 0; column < N; column++)
                m(row, column) = em(row, column);
        return m;
    }

    template<typename scalar>
    Eigen::Quaternion<scalar> toEigen(const quat<scalar>& q)
    {
        return Eigen::Quaternion<scalar>(q.scalar(), q.x(), q.y(), q.z());
    }

    template<typename scalar>
    quat<scalar> fromEigen(const Eigen::Quaternion<scalar>& eq)
    {
        return quat<scalar>(eq.w(), eq.x(), eq.y(), eq.z());
    }
}  // namespace Math
}  // namespace Preon
