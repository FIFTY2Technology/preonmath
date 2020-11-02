/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "vec.h"
#include "matrix.h"
#include "euler.h"
#include "quat.h"

namespace Preon
{
namespace Math
{
    template<typename TOut, typename TIn, size_t D>
    vec<D, TOut> cast(const vec<D, TIn>& v)
    {
        return vec<D, TOut>([&v](size_t i) { return static_cast<TOut>(v[i]); });
    }

    template<typename TOut, typename TIn, size_t M, size_t N>
    matrix<M, N, TOut> cast(const matrix<M, N, TIn>& m)
    {
        return matrix<M, N, TOut>([&m](size_t r, size_t c) { return static_cast<TOut>(m(r, c)); });
    }

    template<typename TOut, typename TIn>
    euler<TOut> cast(const euler<TIn>& e)
    {
        return euler<TOut>(static_cast<TOut>(e.phi()), static_cast<TOut>(e.theta()), static_cast<TOut>(e.psi()));
    }

    template<typename TOut, typename TIn>
    quat<TOut> cast(const quat<TIn>& q)
    {
        return quat<TOut>(static_cast<TOut>(q.scalar()), static_cast<TOut>(q.x()), static_cast<TOut>(q.y()), static_cast<TOut>(q.z()));
    }

}  // namespace Math
}  // namespace Preon
