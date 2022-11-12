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

#include <algorithm>
#include <vector>
#include <cstring>
#ifdef USE_TBB
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include "Defines.hpp"
#elif defined(_OPENMP)
#include <omp.h>
#endif
#include <glm/glm.hpp>

#ifdef TRACY_ENABLE
#include <tracy/Tracy.hpp>
#endif

#include "Util.hpp"
#include "SnapMC.hpp"
#include "SnapMCTable.hpp"

/**
 * Linearly interpolates between p0 and p1 using the weight of the original grid (without the snapped values)
 * @param p0 The first point to interpolate.
 * @param p1 The second point to interpolate.
 * @param weight The weight of the interpolation in the orignal grid (without the snapped values)
 * @return The interpolated point.
 */
glm::vec3 snapBack(const glm::vec3& p0, const glm::vec3& p1, float weight) {
    glm::vec3 p;
    p.x = p0.x + weight * (p1.x - p0.x);
    p.y = p0.y + weight * (p1.y - p0.y);
    p.z = p0.z + weight * (p1.z - p0.z);

    return p;
}

/**
 * Linearly interpolates between p0 and p1 using the weight of the original grid (without the snapped values).
 * @param {vec3} p0 The first normal to interpolate.
 * @param {vec3} p1 The second normal to interpolate.
 * @param {number} weight The weight of the interpolation in the original grid (without the snapped values).
 * @return {vec3} The interpolated  normal.
 */
glm::vec3 snapBackNormal(const glm::vec3& p0, const glm::vec3& p1, float weight) {
    glm::vec3 p = snapBack(p0, p1, weight);
    return glm::normalize(p);
}

/**
 * Linearly interpolates between p0 and p1 using the values of f0, f1 and the iso-level.
 * @param isoLevel The Iso-value of the iso surface to polygonize.
 * @param p0 The first point to interpolate.
 * @param p1 The second point to interpolate.
 * @param f0 The scalar value of the first point.
 * @param f1 The scalar value of the second point.
 * @return The interpolated point.
 */
glm::vec3 vertexInterpIsoSnapMC(float isoLevel, const glm::vec3& p0, const glm::vec3& p1, float f0, float f1) {
    float weight = (isoLevel - f0) / (f1 - f0);
    glm::vec3 p;
    p.x = p0.x + weight * (p1.x - p0.x);
    p.y = p0.y + weight * (p1.y - p0.y);
    p.z = p0.z + weight * (p1.z - p0.z);

    return p;
}

/**
 * Linearly interpolates between p0 and p1 using the values of f0, f1 and the iso-level.
 * @param isoLevel The Iso-value of the iso surface to polygonize.
 * @param p0 The first normal to interpolate.
 * @param p1 The second normal to interpolate.
 * @param f0 The scalar value of the first point.
 * @param f1 The scalar value of the second point.
 * @return The interpolated normal.
 */
glm::vec3 normalInterpIsoSnapMC(float isoLevel, const glm::vec3& p0, const glm::vec3& p1, float f0, float f1) {
    glm::vec3 p = vertexInterpIsoSnapMC(isoLevel, p0, p1, f0, f1);
    return glm::normalize(p);
}


/**
 * Polygonizes a grid cell using SnapMC.
 * @param gridCell The grid cell to polygonize.
 * @param isoLevel The iso value of the iso surface to polygonize.
 * @param snapGrid The grid holding the snapping information.
 * @param i The first index ot the (in relation to the index) first vertex of the cube (z-direction).
 * @param j The second index ot the (in relation to the index) first vertex of the cube (y-direction).
 * @param k The third index ot the (in relation to the index) first vertex of the cube (x-direction).
 * @return An object containing an array of triangle points and an array of vector normals.
 */
void polygonizeSnapMC(
        const GridCell& gridCell, float isoLevel, const SnapGrid& snapGrid, int nx, int ny, int nz, int i, int j, int k,
        std::vector<glm::vec3>& vertexPositions, std::vector<glm::vec3>& vertexNormals) {
    // Check if cube intersects scalar field.
    bool cubeIntersectsScalarfield = false;
    if (gridCell.f[0] < isoLevel) {
        for (int l = 1; l < 8; l++) {
            if (gridCell.f[l] >= isoLevel) {
                cubeIntersectsScalarfield = true;
            }
        }
    }
    else {
        for (int l = 1; l < 8; l++) {
            if (gridCell.f[l] < isoLevel) {
                cubeIntersectsScalarfield = true;
            }
        }
    }
    if (!cubeIntersectsScalarfield) {
        return;
    }
    if (std::any_of(gridCell.f, gridCell.f + 8, [](float val) { return std::isnan(val); })) {
        return;
    }

    // Calculate the table index to find the right array in isoTable.
    int tableIndex = 0;
    for (int i = 0; i < 8; i++) {
        if (gridCell.f[i] > isoLevel) {
            tableIndex += positive[i];
        }
        else if (gridCell.f[i] == isoLevel) {
            tableIndex += equals[i];
        }
    }

    // Create the triangle vertices.
    glm::vec3 isoPoint;
    glm::vec3 normal;
    for (int l = 0; l < isoTableSize[tableIndex]; l++) {
        if (isoTable[tableIndex][l] < 8) {
            // If the isosurface vertex lies directly on a grid vertex, calculate the global indices of the current cube
            // vertex.
            int i_curr = i;
            int j_curr = j;
            int k_curr = k;
            int j0 = isoTable[tableIndex][l];
            for (int d = 0; d < 3; d++) {
                if ((j0 % 2) == 1) {
                    switch (d) {
                        case 0:
                            k_curr++;
                            break;
                        case 1:
                            j_curr++;
                            break;
                        case 2:
                            i_curr++;
                            break;
                    }
                }
                j0 = j0 / 2;
            }

            int idx1d = IDX_GRID(i_curr, j_curr, k_curr);
            if (snapGrid.snapBack[idx1d]) {
                // Uf the current vertex was snapped to.
                isoPoint = snapBack(
                        gridCell.v[isoTable[tableIndex][l]], snapGrid.snapBackTo[idx1d], snapGrid.weights[idx1d]);
                vertexPositions.push_back(isoPoint);
                normal = snapBackNormal(
                        gridCell.n[isoTable[tableIndex][l]], snapGrid.snapBackToNormals[idx1d], snapGrid.weights[idx1d]);
                vertexNormals.push_back(normal);
            }
            else {
                // If the current vertex has the value of the isoLevel without snapping.
                isoPoint = gridCell.v[isoTable[tableIndex][l]];
                vertexPositions.push_back(isoPoint);
                normal = gridCell.n[isoTable[tableIndex][l]];
                vertexNormals.push_back(normal);
            }
        }
        else {
            // If the isosurface vertex lies on a grid edge.
            switch (isoTable[tableIndex][l]) {
                case 8:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[0], gridCell.v[1], gridCell.f[0], gridCell.f[1]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[0], gridCell.n[1], gridCell.f[0], gridCell.f[1]));
                    break;
                case 9:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[0], gridCell.v[2], gridCell.f[0], gridCell.f[2]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[0], gridCell.n[2], gridCell.f[0], gridCell.f[2]));
                    break;
                case 10:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[1], gridCell.v[3], gridCell.f[1], gridCell.f[3]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[1], gridCell.n[3], gridCell.f[1], gridCell.f[3]));
                    break;
                case 11:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[2], gridCell.v[3], gridCell.f[2], gridCell.f[3]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[2], gridCell.n[3], gridCell.f[2], gridCell.f[3]));
                    break;
                case 12:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[0], gridCell.v[4], gridCell.f[0], gridCell.f[4]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[0], gridCell.n[4], gridCell.f[0], gridCell.f[4]));
                    break;
                case 13:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[1], gridCell.v[5], gridCell.f[1], gridCell.f[5]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[1], gridCell.n[5], gridCell.f[1], gridCell.f[5]));
                    break;
                case 14:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[2], gridCell.v[6], gridCell.f[2], gridCell.f[6]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[2], gridCell.n[6], gridCell.f[2], gridCell.f[6]));
                    break;
                case 15:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[3], gridCell.v[7], gridCell.f[3], gridCell.f[7]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[3], gridCell.n[7], gridCell.f[3], gridCell.f[7]));
                    break;
                case 16:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[4], gridCell.v[5], gridCell.f[4], gridCell.f[5]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[4], gridCell.n[5], gridCell.f[4], gridCell.f[5]));
                    break;
                case 17:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[4], gridCell.v[6], gridCell.f[4], gridCell.f[6]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[4], gridCell.n[6], gridCell.f[4], gridCell.f[6]));
                    break;
                case 18:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[5], gridCell.v[7], gridCell.f[5], gridCell.f[7]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[5], gridCell.n[7], gridCell.f[5], gridCell.f[7]));
                    break;
                case 19:
                    isoPoint = vertexInterpIsoSnapMC(
                            isoLevel, gridCell.v[6], gridCell.v[7], gridCell.f[6], gridCell.f[7]);
                    vertexPositions.push_back(isoPoint);
                    vertexNormals.push_back(normalInterpIsoSnapMC(
                            isoLevel, gridCell.n[6], gridCell.n[7], gridCell.f[6], gridCell.f[7]));
                    break;
            }
        }
    }
}

/**
 * Snaps the isoLevel value to the first given vertex of the edge, if the isoLevel has not been already
 * snaped to that vertex. If the isoSurface value on this edge lies nearer than the already snapped value of
 * another edge, the snapBackTo vertex is set to the second vertex given as well as the snapBackToNormal is
 * set to the normal of the second vertex. If the isoSurface value lies exactlyat the same distance as at the
 * already snapped value, the snapBackTo vertex is set to the vertex with the lowest overall index.
 * @param snapGrid The grid where the snaps are happening. An object storing "gridValues"(a 3D array of 1D grid values),
 * "snapBackTo" (a 3D array of 3D grid indices) and "weights" (a 3D array of 1D weights to snap vertices back).
 * @param gridPoints A three-dimensional array of three-dimensional grid points.
 * @param gridNormals A three-dimensional array of normals at the grid points.
 * @param isoLevel The Iso-value of the iso surface to construct.
 * @param weight The ratio to which the original point is between the first vertex and the second vertex.
 * @param i0 The first index of the the first vertex of the edge.
 * @param j0 The second index of the the first vertex of the edge.
 * @param k0 The third index of the the first vertex of the edge.
 * @param i1 The first index of the the second vertex of the edge.
 * @param j1 The second index of the the second vertex of the edge.
 * @param k1 The third index of the the second vertex of the edge.
 */
void snapToVertex(
        SnapGrid& snapGrid, const glm::vec3* gridPoints, const glm::vec3* gridNormals, float isoLevel, float weight,
        int nx, int ny, int nz, int i0, int j0, int k0, int i1, int j1, int k1) {
    int idx0 = IDX_GRID(i0, j0, k0);
    int idx1 = IDX_GRID(i1, j1, k1);
    if (snapGrid.gridValues[idx0] != isoLevel) {
        // If no snap happend to the first vertex already.
        snapGrid.gridValues[idx0] = isoLevel;
        snapGrid.snapBack[idx0] = true;
        snapGrid.snapBackTo[idx0] = gridPoints[idx1];
        snapGrid.snapBackToNormals[idx0] = gridNormals[idx1];
        snapGrid.weights[idx0] = weight;
    } else if (snapGrid.snapBack[idx0] && snapGrid.weights[idx0] > weight) {
        // If the isoSurface value on the already snapped edge is further away than on the current one.
        snapGrid.snapBackTo[idx0] = gridPoints[idx1];
        snapGrid.snapBackToNormals[idx0] = gridNormals[idx1];
        snapGrid.weights[idx0] = weight;
    }
}

/**
 * Snaps the isoLevel value at an edge to one of the vertices if the isosurface is near enough to one of them.
 * @param cartesianGrid The original grid. An object storing "gridPoints" (a 3D array of 3D grid points) and
 * "gridValues" (a 3D array of 1D grid values).
 * @param snapGrid The grid where the snaps are happening. An object storing "gridValues"(a 3D array of 1D grid values),
 * "snapBackTo" (a 3D array of 3D grid indices) and "weights" (a 3D array of 1D weights to snap vertices back).
 * @param isoLevel The Iso-value of the iso surface to construct.
 * @param i0 The first index of the the first vertex of the edge.
 * @param j0 The second index of the the first vertex of the edge.
 * @param k0 The third index of the the first vertex of the edge.
 * @param i1 The first index of the the second vertex of the edge.
 * @param j1 The second index of the the second vertex of the edge.
 * @param k1 The third index of the the second vertex of the edge.
 */
void snapAtEdge(
        const float* cartesianGrid, const glm::vec3* gridPoints, const glm::vec3* gridNormals,
        SnapGrid& snapGrid, float isoLevel, const float gamma,
        int nx, int ny, int nz, int i0, int j0, int k0, int i1, int j1, int k1) {
    float epsilon = 0.00001f;

    // Weight or distance.
    float weight0 = gamma;
    float weight1 = gamma;

    int idx0 = IDX_GRID(i0, j0, k0);
    int idx1 = IDX_GRID(i1, j1, k1);
    if ((cartesianGrid[idx0] < isoLevel && cartesianGrid[idx1] > isoLevel)
            || (cartesianGrid[idx0] > isoLevel && cartesianGrid[idx1] < isoLevel)) {
        // If this is a +/- edge.
        float value_difference = cartesianGrid[idx1] - cartesianGrid[idx0];

        if (value_difference > epsilon || value_difference < -epsilon) {
            weight0 = (cartesianGrid[idx1] - isoLevel) / value_difference;
            weight1 = (isoLevel - cartesianGrid[idx0]) / value_difference;
        }
        else {
            weight0 = 0.5f;
            weight1 = 0.5f;
        }
    }
    if (weight1 < gamma) {
        // Snap to vertex i0, j0, k0.
        snapToVertex(snapGrid, gridPoints, gridNormals, isoLevel, weight1, nx, ny, nz, i0, j0, k0, i1, j1, k1);

    }
    else if (weight0 < gamma) {
        // Snap to vertex i1, j1, k1.
        snapToVertex(snapGrid, gridPoints, gridNormals, isoLevel, weight0, nx, ny, nz, i1, j1, k1, i0, j0, k0);
    }
}

/**
 * Constructs a grid with snapped values from an cartesian grid.
 * @param cartesianGrid The original grid. An object storing "gridPoints" (a 3D array of 3D grid points) and
 * "gridValues" (a 3D array of 1D grid values).
 * @param isoLevel The Iso-value of the iso surface to construct.
 * @return An object storing "gridValues"(a 3D array of 1D grid values),
 * "snapBackTo" (a 3D array of 3D grid indices) and "weights" (a 3D array of 1D weights to snap vertices back).
 */
SnapGrid constructCartesianSnapGridScalarField(
        SnapGrid& snapGrid, const float* cartesianGrid, const glm::vec3* gridPoints, const glm::vec3* gridNormals,
        float isoLevel, const float gamma, int nx, int ny, int nz) {
    memcpy(snapGrid.gridValues, cartesianGrid, sizeof(float) * nx * ny * nz);

    // Go over all vertices of the grid (just not the last ones in the respective direction since we want to got over
    // the edges).
#ifdef USE_TBB
    tbb::parallel_for(tbb::blocked_range<int>(0, nz - 1), [&](auto const& r) {
        for (int k = r.begin(); k != r.end(); k++) {
#else
#ifdef _MSC_VER
    #pragma omp parallel for shared(snapGrid, cartesianGrid, gridPoints, gridNormals, isoLevel, nx, ny, nz) \
    firstprivate(gamma)
#else
    #pragma omp parallel for shared(snapGrid, cartesianGrid, gridPoints, gridNormals, isoLevel, nx, ny, nz) \
    firstprivate(gamma), default(none)
#endif
    for (int k = 0; k < nz - 1; k++) {
#endif
        for (int j = 0; j < ny - 1; j++) {
            for (int i = 0; i < nx - 1; i++) {
                snapAtEdge(
                        cartesianGrid, gridPoints, gridNormals, snapGrid,
                        isoLevel, gamma, nx, ny, nz, i, j, k, i + 1, j, k);
                snapAtEdge(
                        cartesianGrid, gridPoints, gridNormals, snapGrid,
                        isoLevel, gamma, nx, ny, nz, i, j, k, i, j + 1, k);
                snapAtEdge(
                        cartesianGrid, gridPoints, gridNormals, snapGrid,
                        isoLevel, gamma, nx, ny, nz, i, j, k, i, j, k + 1);
            }
        }
    }
#ifdef USE_TBB
    });
#endif

    return snapGrid;
}

void polygonizeSnapMC(
        const float* voxelGrid, int nx, int ny, int nz, float dx, float dy, float dz, float isoLevel, const float gamma,
        std::vector<glm::vec3>& vertexPositions, std::vector<glm::vec3>& vertexNormals) {
#ifdef TRACY_ENABLE
    ZoneScoped;
#endif

    int numCellsX = nx - 1;
    int numCellsY = ny - 1;
    int numCellsZ = nz - 1;

    auto* gridPoints = new glm::vec3[nx * ny * nz];
    auto* gridNormals = new glm::vec3[nx * ny * nz];
#ifdef USE_TBB
    tbb::parallel_for(tbb::blocked_range<int>(0, nz), [&](auto const& r) {
        for (int z = r.begin(); z != r.end(); z++) {
#else
#ifdef _MSC_VER
    #pragma omp parallel for shared(voxelGrid, gridPoints, gridNormals, nx, ny, nz, dx, dy, dz)
#else
    #pragma omp parallel for default(none) shared(voxelGrid, gridPoints, gridNormals, nx, ny, nz, dx, dy, dz)
#endif
    for (int z = 0; z < nz; z++) {
#endif
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                // Compute the normal vector.
                glm::vec3 n = computeNormal(voxelGrid, nx, ny, nz, dx, dy, dz, glm::ivec3(x, y, z));

                int idx = IDX_GRID(x, y, z);
                gridPoints[idx] = glm::vec3(float(x) * dx, float(y) * dy, float(z) * dz);
                gridNormals[idx] = glm::vec3(n[0], n[1], n[2]);
            }
        }
    }
#ifdef USE_TBB
    });
#endif

    SnapGrid snapGrid{};
    snapGrid.gridValues = new float[nx * ny * nz];
    snapGrid.snapBack = new bool[nx * ny * nz];
    snapGrid.snapBackTo = new glm::vec3[nx * ny * nz];
    snapGrid.snapBackToNormals = new glm::vec3[nx * ny * nz];
    snapGrid.weights = new float[nx * ny * nz];
    memset(snapGrid.snapBack, 0, sizeof(bool) * nx * ny * nz);

    constructCartesianSnapGridScalarField(snapGrid, voxelGrid, gridPoints, gridNormals, isoLevel, gamma, nx, ny, nz);

#ifdef USE_TBB
    VertexNormalArrayBlock geometryBlock = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, numCellsZ), VertexNormalArrayBlock(),
            [&](tbb::blocked_range<int> const& r, VertexNormalArrayBlock init) -> VertexNormalArrayBlock {
                std::vector<glm::vec3>& vertexPositionsLocal = init.vertexPositionsLocal;
                std::vector<glm::vec3>& vertexNormalsLocal = init.vertexNormalsLocal;
                for (int z = r.begin(); z != r.end(); z++) {
#else
#ifdef _MSC_VER
    #pragma omp parallel shared(numCellsX, numCellsY, numCellsZ, nx, ny, nz, dx, dy, dz, isoLevel) \
    shared(voxelGrid, snapGrid, vertexPositions, vertexNormals)
#else
    #pragma omp parallel default(none) shared(numCellsX, numCellsY, numCellsZ, nx, ny, nz, dx, dy, dz, isoLevel) \
    shared(voxelGrid, snapGrid, vertexPositions, vertexNormals)
#endif
    {
        std::vector<glm::vec3> vertexPositionsLocal;
        std::vector<glm::vec3> vertexNormalsLocal;

        #pragma omp for
        for (int z = 0; z < numCellsZ; z++) {
#endif
            for (int y = 0; y < numCellsY; y++) {
                for (int x = 0; x < numCellsX; x++) {
                    GridCell gridCell;

                    for (int l = 0; l < 8; l++) {
                        glm::ivec3 gridIndex(x, y, z);
                        if (l == 1 || l == 3 || l == 5 || l == 7) {
                            gridIndex[0] += 1;
                        }
                        if (l == 2 || l == 3 || l == 6 || l == 7) {
                            gridIndex[1] += 1;
                        }
                        if (l == 4 || l == 5 || l == 6 || l == 7) {
                            gridIndex[2] += 1;
                        }

                        // Compute the normal vector.
                        glm::vec3 n = computeNormal(voxelGrid, nx, ny, nz, dx, dy, dz, gridIndex);

                        gridCell.v[l] = glm::vec3{
                            float(gridIndex[0]) * dx, float(gridIndex[1]) * dy, float(gridIndex[2]) * dz};
                        gridCell.n[l] = glm::vec3{n[0], n[1], n[2]};
                        gridCell.f[l] = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2])];
                    }

                    polygonizeSnapMC(
                            gridCell, isoLevel, snapGrid, nx, ny, nz, x, y, z,
                            vertexPositionsLocal, vertexNormalsLocal);
                }
            }
        }
#ifdef USE_TBB
                return init;
            },
            [&](VertexNormalArrayBlock lhs, VertexNormalArrayBlock rhs) -> VertexNormalArrayBlock {
                VertexNormalArrayBlock blockOut = std::move(lhs);
                blockOut.vertexPositionsLocal.insert(
                        blockOut.vertexPositionsLocal.end(),
                        std::make_move_iterator(rhs.vertexPositionsLocal.begin()),
                        std::make_move_iterator(rhs.vertexPositionsLocal.end()));
                blockOut.vertexNormalsLocal.insert(
                        blockOut.vertexNormalsLocal.end(),
                        std::make_move_iterator(rhs.vertexNormalsLocal.begin()),
                        std::make_move_iterator(rhs.vertexNormalsLocal.end()));
                return blockOut;
            });
    vertexPositions = std::move(geometryBlock.vertexPositionsLocal);
    vertexNormals = std::move(geometryBlock.vertexNormalsLocal);
#else

#ifdef _OPENMP
        #pragma omp for ordered schedule(static, 1)
        for (int threadIdx = 0; threadIdx < omp_get_num_threads(); ++threadIdx) {
#else
        for (int threadIdx = 0; threadIdx < 1; ++threadIdx) {
#endif
            #pragma omp ordered
            {
                vertexPositions.reserve(vertexPositions.size() + vertexPositionsLocal.size());
                for (auto& vertexPosition : vertexPositionsLocal) {
                    vertexPositions.push_back(vertexPosition);
                }

                vertexNormals.reserve(vertexNormals.size() + vertexNormalsLocal.size());
                for (auto& vertexNormal : vertexNormalsLocal) {
                    vertexNormals.push_back(vertexNormal);
                }
            }
        }
    }
#endif

    delete[] snapGrid.gridValues;
    delete[] snapGrid.snapBack;
    delete[] snapGrid.snapBackTo;
    delete[] snapGrid.snapBackToNormals;
    delete[] snapGrid.weights;

    delete[] gridPoints;
    delete[] gridNormals;
}

void polygonizeSnapMC(
        const float* voxelGrid, int nx, int ny, int nz, float isoLevel, const float gamma,
        std::vector<glm::vec3>& vertexPositions, std::vector<glm::vec3>& vertexNormals) {
    polygonizeSnapMC(voxelGrid, nx, ny, nz, 1.0f, 1.0f, 1.0f, isoLevel, gamma, vertexPositions, vertexNormals);
}
