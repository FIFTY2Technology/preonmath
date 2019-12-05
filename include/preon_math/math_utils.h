/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_MATH_UTILS_H
#define PREONMATH_MATH_UTILS_H

#include "compile_helper.h"

#include "vec.h"
#include "scalar_simd.h"

#include <cmath>

namespace Preon
{
    namespace Math
    {
        // Todo: Move all math utilities into a class or namespace.
        namespace MathUtils
        {
            //! Rounds a floating point number to its nearest integral value.
            template<class T>
            T round(T x)
            {
                return std::round(x);
            }

            //! Returns the logarithm of the given number to the given base.
            template<class Real>
            Real logn(Real base, Real number)
            {
                return std::log(number) / std::log(base);
            }

            //! Returns the power to the given base that is closest to the given number.
            template<class Real>
            Real nearestPower(Real base, Real number)
            {
                return std::pow(base, MathUtils::round(logn(base, number)));
            }

            //! Returns the next power to the given base that is greater or equal to number.
            template<class Real>
            Real nextPower(Real base, Real number)
            {
                return std::pow(base, std::ceil(logn(base, number)));
            }

            /*!
             * Returns the next greater power of 2 for the integer n and the corresponding exponent.
             * No negative input.
             */
            template<typename IntType>
            std::tuple<IntType, IntType> nextPowerOfTwoWithExpo(IntType n)
            {
                THROW_EXCEPTION(n < 0, std::invalid_argument("negative value passed to nextPowerOfTwo!"))
                // Round up to the next power of two.
                std::tuple<IntType, IntType> out(1, 0);
                while (std::get<0>(out) < n)
                {
                    std::get<0>(out) <<= 1;
                    std::get<1>(out)++;
                }
                return out;
            }

            /*!
             * \brief Returns the next greater power of 2 of \a n.
             * \pre No negative input.
             * \param n Input variable for which the next greater power of 2 will be returned.
             * \return Returns the next greater power of 2 of \a n.
             */
            template<typename IntType>
            IntType nextPowerOfTwo(IntType n)
            {
                return std::get<0>(nextPowerOfTwoWithExpo(n));
            }

            // Returns the ceil(a / b)
            template<class Int>
            Int ceilDivision(Int a, Int b)
            {
                return (a + b - 1) / b;
            }

            //! Clamps the given value so that it is in range [min, max].
            template<class scalar>
            scalar clamp(scalar x, scalar min, scalar max)
            {
                return std::max(min, std::min(max, x));
            }

            //! Clamps the given value so that it is in range [0, 1].
            template<class scalar>
            scalar clampToZeroOne(scalar x, enable_if_non_simd<scalar, void*> = nullptr)
            {
                return clamp(x, scalar{0}, scalar{1});
            }

            //! Clamps the given value so that it is in range [0, 1].
            template<class scalar>
            scalar clampToZeroOne(scalar x, enable_if_simd<scalar, void*> = nullptr)
            {
                static const float_simd one = Simd::make<float>(1.0f);
                return Simd::and_mask(Simd::choose(x > one, one, x), x >= Simd::zero<float>());
            }

            //! Maps x in range [xMin, xMax] linearly into the range [0, 1].
            template <typename scalar>
            scalar mapToZeroOneIntervalLinear(scalar x, scalar xMin, scalar xMax)
            {
                return (x - xMin) / (xMax - xMin);
            }

            //! Maps x in range [xMin, xMax] linearly into the range [0, 1]. Clamps the result to [0, 1] if x is not in the range [xMin, xMax].
            template <typename scalar>
            scalar mapToZeroOneIntervalLinearClamped(scalar x, scalar xMin, scalar xMax)
            {
                return clampToZeroOne(mapToZeroOneIntervalLinear(x, xMin, xMax));
            }

            //! Returns the value of the gaussian probability density function with the given mean and standard deviation.
            //! See https://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c.
            template <typename scalar>
            scalar normalPdf(scalar value, scalar mean, scalar standardDeviation)
            {
                static const scalar inv_sqrt_2pi = 0.3989422804014327;
                scalar a = (value - mean) / standardDeviation;

                return (inv_sqrt_2pi / standardDeviation) * std::exp(-scalar(0.5) * a * a);
            }

            //! Takes two uniformly-distributed random numbers in [0, 1] and returns two gauss-distributed random values with the given mean and deviation.
            //! See https://en.wikipedia.org/wiki/Box-Muller_transform.
            template <typename scalar>
            vec<2, scalar> boxMullerTransform(const vec<2, scalar>& rand, scalar mean, scalar deviation)
            {
                scalar a = DTwoPi * rand[0];
                scalar b = std::sqrt(-2 * std::log(rand[1]));
                float x = std::cos(a) * b;
                float y = std::sin(a) * b;
                return vec<2, scalar>(mean + deviation * x, mean + deviation *y);
            }
        }
    }  // namespace Math
}  // namespace Preon

#endif  // PREONMATH_MATH_UTILS_H
