/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include "compile_helper.h"

#include "vec.h"

namespace Preon
{
namespace Math
{
#define INVALID_POSITION vec3f(F_INFINITY)
}  // namespace Math
}  // namespace Preon

// Hash Function.
namespace std
{
template<>
struct hash<Preon::Math::vec3i>
{
private:
    static const size_t PRIME1 = 73856093;
    static const size_t PRIME2 = 19349663;
    static const size_t PRIME3 = 83492791;

public:
    size_t operator()(const Preon::Math::vec3i& v) const
    {
        /*static const int size = 4;//sizeof(size_t)
            int iShift=(size*8)/Vector3i::dimension;
            size_t iMask=((size_t)1<<iShift)-1;
            return ((v.x)&(iMask)) | ((v.y&iMask)<<iShift) | ((v.z&iMask)<<(iShift<<1));*/
        return ((PRIME1 * v.x()) ^ (PRIME2 * v.y()) ^ (PRIME3 * v.z()));
    }
};

template<>
struct hash<Preon::Math::vec3f>
{
public:
    size_t operator()(const Preon::Math::vec3f& v) const
    {
        return hash<Preon::Math::vec3i>()(Preon::interpretAs<Preon::Math::vec3i>(v));
    }
};

}  // namespace std