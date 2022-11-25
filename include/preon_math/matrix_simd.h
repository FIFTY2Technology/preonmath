/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "matrix.h"
#include "vec_simd.h"

#ifdef PREONMATH_COMPILER_MSVC
    #include <intrin.h>
#else
    #include <immintrin.h>
#endif

namespace Preon
{
namespace Math
{
    template<PrMathSize M, PrMathSize N, typename T = float>
    using matrix_simd = matrix<M, N, typename Simd::Register<T>::type>;

    namespace Simd
    {
        template<PrMathSize M, PrMathSize N, typename T>
        PREONMATH_FORCEINLINE void setValues(matrix_simd<M, N, T>* out, const matrix<M, N, T>& m0, const matrix<M, N, T>& m1, const matrix<M, N, T>& m2, const matrix<M, N, T>& m3)
        {
            for (PrMathSize c = 0; c < N; c++)
                setVecs(&out->column(c), m0.column(c), m1.column(c), m2.column(c), m3.column(c));
        }

        template<PrMathSize M, PrMathSize N, typename T, class Getter>
        PREONMATH_FORCEINLINE void setValues(matrix_simd<M, N, T>* out, const Getter& func)
        {
            setValues(out, func(0), func(1), func(2), func(3));
        }

        template<PrMathSize M, PrMathSize N, typename T, class Getter>
        PREONMATH_FORCEINLINE void setValues(matrix_simd<M, N, T>* out, const Getter& func, uint32_t numElements, const matrix<M, N, T>& fillValue)
        {
            setValues(out, func(0), (numElements > 1) ? func(1) : fillValue, (numElements > 2) ? func(2) : fillValue, (numElements > 3) ? func(3) : fillValue);
        }

        template<PrMathSize M, PrMathSize N, typename T>
        PREONMATH_FORCEINLINE matrix<M, N, T> hSum(const matrix_simd<M, N, T>& m)
        {
            matrix<M, N, T> out;
            for (PrMathSize c = 0; c < N; c++)
                for (PrMathSize r = 0; r < M; r++)
                    out(r, c) = hSum(m(r, c));
            return out;
        }
        // We have the following two methods so that we do not need to provide the
        // template parameters when using hSum for matrices since for automatic
        // template deducation, the function template arguments must be exactly
        // the same as the ones from the given matrix class. This is not the case
        // for the function above since it takes non-simd and the matrix has simd
        //  as last template argument.
        template<PrMathSize M, PrMathSize N>
        PREONMATH_FORCEINLINE matrix<M, N, float> hSum(const matrix<M, N, float_simd>& m)
        {
            return hSum<M, N, float>(m);
        }
        template<PrMathSize M, PrMathSize N>
        PREONMATH_FORCEINLINE matrix<M, N, double> hSum(const matrix<M, N, double_simd>& m)
        {
            return hSum<M, N, double>(m);
        }
    }  // namespace Simd
}  // namespace Math
}  // namespace Preon
