/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2022, Christoph Neuhauser
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

#include <glm/glm.hpp>

#include "Util.hpp"

glm::vec3 computeNormal(
        const float* voxelGrid, int nx, int ny, int nz, float dx, float dy, float dz, const glm::ivec3& gridIndex) {
    int numCellsX = nx - 1;
    int numCellsY = ny - 1;
    int numCellsZ = nz - 1;

    float hdx = 2.0f;
    float valueXm, valueXp;
    if (gridIndex[0] > 0) {
        valueXm = voxelGrid[IDX_GRID(gridIndex[0] - 1, gridIndex[1], gridIndex[2])];
    }
    if (gridIndex[0] <= 0 || std::isnan(valueXm)) {
        valueXm = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2])];
        hdx -= 1.0f;
    }
    if (gridIndex[0] < numCellsX) {
        valueXp = voxelGrid[IDX_GRID(gridIndex[0] + 1, gridIndex[1], gridIndex[2])];
    }
    if (gridIndex[0] >= numCellsX || std::isnan(valueXp)) {
        valueXp = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2])];
        hdx -= 1.0f;
    }
    hdx *= dx;
    float normalX = (valueXp - valueXm) / hdx;

    float hdy = 2.0f;
    float valueYm, valueYp;
    if (gridIndex[1] > 0) {
        valueYm = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1] - 1, gridIndex[2])];
    }
    if (gridIndex[1] <= 0 || std::isnan(valueYm)) {
        valueYm = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2])];
        hdy -= 1.0f;
    }
    if (gridIndex[1] < numCellsY) {
        valueYp = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1] + 1, gridIndex[2])];
    }
    if (gridIndex[1] >= numCellsY || std::isnan(valueYp)) {
        valueYp = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2])];
        hdy -= 1.0f;
    }
    hdy *= dy;
    float normalY = (valueYp - valueYm) / hdy;

    float hdz = 2.0f;
    float valueZm, valueZp;
    if (gridIndex[2] > 0) {
        valueZm = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2] - 1)];
    }
    if (gridIndex[2] <= 0 || std::isnan(valueZm)) {
        valueZm = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2])];
        hdz -= 1.0f;
    }
    if (gridIndex[2] < numCellsZ) {
        valueZp = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2] + 1)];
    }
    if (gridIndex[2] >= numCellsZ || std::isnan(valueZp)) {
        valueZp = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2])];
        hdz -= 1.0f;
    }
    hdz *= dz;
    float normalZ = (valueZp - valueZm) / hdz;

    /*float normalX =
            (voxelGrid[IDX_GRID(std::min(gridIndex[0] + 1, numCellsX), gridIndex[1], gridIndex[2])]
             - voxelGrid[IDX_GRID(std::max(gridIndex[0] - 1, 0), gridIndex[1], gridIndex[2])]) /
            (-2.0f * h[0]);
    float normalY =
            (voxelGrid[IDX_GRID(gridIndex[0], std::min(gridIndex[1] + 1, numCellsY), gridIndex[2])]
             - voxelGrid[IDX_GRID(gridIndex[0], std::max(gridIndex[1] - 1, 0), gridIndex[2])]) /
            (-2.0f * h[1]);
    float normalZ =
            (voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], std::min(gridIndex[2] + 1, numCellsZ))]
             - voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], std::max(gridIndex[2] - 1, 0))]) /
            (-2.0f * h[2]);*/
    return glm::normalize(glm::vec3(normalX, normalY, normalZ));
}
