/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2021, Christoph Neuhauser, Stefan Haas
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

#ifndef ISOSURFACECPP_SNAPMC_HPP
#define ISOSURFACECPP_SNAPMC_HPP

/**
 * The following code is based on the SnapMC isosurface extraction algorithm described in the following paper.
 *
 * Sundaresan Raman and Rephael Wenger. Quality Isosurface Mesh Generation Using an Extended Marching Cubes Lookup
 * Table. Computer Graphics Forum , 27(3):791-798, 2008.
 */

#include "MarchingCubes.hpp"

struct SnapGrid {
    float* gridValues;
    bool* snapBack;
    glm::vec3* snapBackTo;
    glm::vec3* snapBackToNormals;
    float* weights;
};

void polygonizeSnapMC(
        const GridCell& gridCell, float isoLevel, const SnapGrid& snapGrid, int nx, int ny, int nz, int i, int j, int k,
        std::vector<glm::vec3>& vertexPositions, std::vector<glm::vec3>& vertexNormals);

void polygonizeSnapMC(
        const float* voxelGrid, int nx, int ny, int nz, float dx, float dy, float dz, float isoLevel, const float gamma,
        std::vector<glm::vec3>& vertexPositions, std::vector<glm::vec3>& vertexNormals);

void polygonizeSnapMC(
        const float* voxelGrid, int nx, int ny, int nz, float isoLevel, const float gamma,
        std::vector<glm::vec3>& vertexPositions, std::vector<glm::vec3>& vertexNormals);

#endif //ISOSURFACECPP_SNAPMC_HPP
