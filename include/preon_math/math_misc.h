/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once
#ifndef PREONMATH_MATH_MISC_H
#define PREONMATH_MATH_MISC_H

#include "compile_helper.h"

#include "vec3.h"

namespace Preon
{
    namespace Math
    {
        struct MathMisc
        {
            static float rayPlaneIntersection(const vec3f& rayOrigin, const vec3f& rayDir, const vec3f& planeNormal, const vec3f& planePos)
            {
                float denom = vec3f::dotProduct(rayDir, planeNormal);
                if (std::fabs(denom) <= 0.000001f) denom = 0.000001f;
                return vec3f::dotProduct(planeNormal, planePos - rayOrigin) / denom;
            }

            static vec3f projectPointOnLine(const vec3f& p, const vec3f& line1, const vec3f& line2, float& t)
            {
                vec3f rayDir = line2 - line1;
                t = rayPlaneIntersection(line1, rayDir, rayDir, p);
                return line1 + rayDir * t;
            }

            /// Projects a triangle on an axis, retrieves interval (tMin, tMax).
            static void projectTriangleOnAxis(const vec3f& axis, const vec3f& p1, const vec3f& p2, const vec3f& p3, float& tMin, float& tMax)
            {
                tMin = vec3f::dotProduct(axis, p1);
                tMax = tMin;
                float t = vec3f::dotProduct(axis, p2);
                if (t < tMin) tMin = t;
                else if (t > tMax) tMax = t;
                t = vec3f::dotProduct(axis, p3);
                if (t < tMin) tMin = t;
                else if (t > tMax) tMax = t;
            }

            //! Projects an AABB on an axis, retrieves interval (tMin, tMax).
            //! aabbMinMax must point to a location where the aabb min and max positions are stored (continuously and in this order).
            static void projectAABBOnAxis(const vec3f& axis, const vec3f* aabbMinMax, float& tMin, float& tMax)
            {
                tMin = 0.0f;
                tMax = 0.0f;
                for (uint d = 0; d < 3; d++)
                {
                    bool negative = axis[d] < 0.0f;
                    tMin += aabbMinMax[negative][d] * axis[d];
                    tMax += aabbMinMax[!negative][d] * axis[d];
                }
            }

            //! Projects a triangle on a unit axis (axis must be 0, 1 or 2 identifying x, y or z axis), retrieves interval (tMin, tMax).
            PREONMATH_FORCEINLINE static void projectTriangleOnUnitAxis(int axis, const vec3f& p1, const vec3f& p2, const vec3f& p3, float& tMin, float& tMax)
            {
                tMin = p1[axis];
                tMax = tMin;
                float t = p2[axis];
                if (t < tMin) tMin = t;
                else if (t > tMax) tMax = t;
                t = p3[axis];
                if (t < tMin) tMin = t;
                else if (t > tMax) tMax = t;
            }

            //! Projects two points on an axis, retrieves interval (tMin, tMax).
            PREONMATH_FORCEINLINE static void projectPointsOnAxis(const vec3f& axis, const vec3f& p1, const vec3f& p2, float& tMin, float& tMax)
            {
                tMin = vec3f::dotProduct(axis, p1);
                tMax = tMin;
                float t = vec3f::dotProduct(axis, p2);
                if (t < tMin) tMin = t;
                else if (t > tMax) tMax = t;
            }

            template<class T>
            static bool intervalDoesNotOverlap(T i1Min, T i1Max, T i2Min, T i2Max)
            {
                return (i1Min > i2Max) || (i2Min > i1Max);
            }

            static float aabbPointSquaredDistance(const vec3f& aabbMin, const vec3f& aabbMax, const vec3f& point)
            {
                float squaredDist = 0;

                // Go over all three dimensions.
                for (size_t d = 0; d < 3; ++d)
                {
                    if (point[d] < aabbMin[d])
                    {
                        float e = point[d] - aabbMin[d];
                        squaredDist += e * e;
                    }
                    else if (point[d] > aabbMax[d])
                    {
                        float e = point[d] - aabbMax[d];
                        squaredDist += e * e;
                    }
                }
                return squaredDist;
            }

            static float aabbPointMaxSquaredDistance(const vec3f& aabbMin, const vec3f& aabbMax, const vec3f& point)
            {
                float squaredDist = 0;
                squaredDist += std::max(square(point.x() - aabbMax.x()), square(point.x() - aabbMin.x()));
                squaredDist += std::max(square(point.y() - aabbMax.y()), square(point.y() - aabbMin.y()));
                squaredDist += std::max(square(point.z() - aabbMax.z()), square(point.z() - aabbMin.z()));
                return squaredDist;
            }

            static void projectPointOnAABB(const vec3f& aabbMin, const vec3f& aabbMax, vec3f& point)
            {
                if (point.x() < aabbMin.x())
                    point.x() = aabbMin.x();
                else if (point.x() > aabbMax.x())
                    point.x() = aabbMax.x();

                if (point.y() < aabbMin.y())
                    point.y() = aabbMin.y();
                else if (point.y() > aabbMax.y())
                    point.y() = aabbMax.y();

                if (point.z() < aabbMin.z())
                    point.z() = aabbMin.z();
                else if (point.z() > aabbMax.z())
                    point.z() = aabbMax.z();
            }

            template<class T>
            static T square(const T& val)
            {
                return val * val;
            }

            template<class T>
            static T linearInterpolation(const T& val1, const T& val2, float weight)
            {
                return val1 * weight + val2 * (1.0f - weight);
            }

            /// Performs trilinear interpolation.
            template<class T>
            static T  trilinearInterpolation(const T* cornerValues, float* weights)
            {
                return cornerValues[0] * (1 - weights[0]) * (1 - weights[1]) * (1 - weights[2])
                       + cornerValues[1] * (1 - weights[0]) * (1 - weights[1]) * weights[2]
                       + cornerValues[2] * (1 - weights[0]) * weights[1] * (1 - weights[2])
                       + cornerValues[3] * (1 - weights[0]) * weights[1] * weights[2]
                       + cornerValues[4] * weights[0] * (1 - weights[1]) * (1 - weights[2])
                       + cornerValues[5] * weights[0] * (1 - weights[1]) * weights[2]
                       + cornerValues[6] * weights[0] * weights[1] * (1 - weights[2])
                       + cornerValues[7] * weights[0] * weights[1] * weights[2];
            }

            template<class T>
            static T trilinearInterpolation(const T* cornerValues, const vec3f& weights)
            {
                return cornerValues[0] * (1 - weights[0]) * (1 - weights[1]) * (1 - weights[2])
                       + cornerValues[1] * (1 - weights[0]) * (1 - weights[1]) * weights[2]
                       + cornerValues[2] * (1 - weights[0]) * weights[1] * (1 - weights[2])
                       + cornerValues[3] * (1 - weights[0]) * weights[1] * weights[2]
                       + cornerValues[4] * weights[0] * (1 - weights[1]) * (1 - weights[2])
                       + cornerValues[5] * weights[0] * (1 - weights[1]) * weights[2]
                       + cornerValues[6] * weights[0] * weights[1] * (1 - weights[2])
                       + cornerValues[7] * weights[0] * weights[1] * weights[2];
            }
        };
    }  // namespace Math
}  // namespace Preon

#endif  // PREONMATH_MATH_MISC_H
