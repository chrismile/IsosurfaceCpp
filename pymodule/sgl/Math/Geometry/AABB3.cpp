/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2017, Christoph Neuhauser
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "AABB3.hpp"

namespace sgl {

bool AABB3::intersects(const AABB3& otherAABB) const {
    if (max.x < otherAABB.min.x || min.x > otherAABB.max.x
        || max.y < otherAABB.min.y || min.y > otherAABB.max.y
        || max.z < otherAABB.min.z || min.z > otherAABB.max.z) {
        return false;
    }
    return true;
}

void AABB3::combine(const AABB3& otherAABB) {
    if (otherAABB.min.x < min.x)
        min.x = otherAABB.min.x;
    if (otherAABB.min.y < min.y)
        min.y = otherAABB.min.y;
    if (otherAABB.min.z < min.z)
        min.z = otherAABB.min.z;

    if (otherAABB.max.x > max.x)
        max.x = otherAABB.max.x;
    if (otherAABB.max.y > max.y)
        max.y = otherAABB.max.y;
    if (otherAABB.max.z > max.z)
        max.z = otherAABB.max.z;
}

void AABB3::combine(const glm::vec3& pt) {
    if (pt.x < min.x)
        min.x = pt.x;
    if (pt.y < min.y)
        min.y = pt.y;
    if (pt.z < min.z)
        min.z = pt.z;

    if (pt.x > max.x)
        max.x = pt.x;
    if (pt.y > max.y)
        max.y = pt.y;
    if (pt.z > max.z)
        max.z = pt.z;
}

bool AABB3::contains(const glm::vec3& pt) const {
    return pt.x >= min.x && pt.y >= min.y && pt.z >= min.z && pt.x <= max.x && pt.y <= max.y && pt.z <= max.z;
}

}
