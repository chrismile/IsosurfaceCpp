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

#include <vector>
#include <cstring>
#include <glm/glm.hpp>

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

// For indexing the 3D arrays.
#define IDX_GRID(x, y, z) int((x) + ((y) + (z) * ny) * nx)

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
        const GridCell& gridCell, float isoLevel, SnapGrid& snapGrid, int nx, int ny, int nz, int i, int j, int k,
        std::vector<glm::vec3>& vertexPositions, std::vector<glm::vec3>& vertexNormals) {
    // Check if cube intersects scalar field.
    bool cubeIntersectsScalarfield = false;
    if (gridCell.f[0] < isoLevel) {
        for (int j = 1; j < 8; j++) {
            if (gridCell.f[j] >= isoLevel) {
                cubeIntersectsScalarfield = true;
            }
        }
    }
    else {
        for (int j = 1; j < 8; j++) {
            if (gridCell.f[j] < isoLevel) {
                cubeIntersectsScalarfield = true;
            }
        }
    }
    if (!cubeIntersectsScalarfield) {
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
                isoPoint = snapBack(gridCell.v[isoTable[tableIndex][l]], snapGrid.snapBackTo[idx1d], snapGrid.weights[idx1d]);
                vertexPositions.push_back(isoPoint);
                normal = snapBackNormal(gridCell.n[isoTable[tableIndex][l]], snapGrid.snapBackToNormals[idx1d], snapGrid.weights[idx1d]);
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
    if (snapGrid.gridValues[IDX_GRID(i0, j0, k0)] != isoLevel) {
        // If no snap happend to the first vertex already.
        snapGrid.gridValues[IDX_GRID(i0, j0, k0)] = isoLevel;
        snapGrid.snapBack[IDX_GRID(i0, j0, k0)] = true;
        snapGrid.snapBackTo[IDX_GRID(i0, j0, k0)] = gridPoints[IDX_GRID(i1, j1, k1)];
        snapGrid.snapBackToNormals[IDX_GRID(i0, j0, k0)] = gridNormals[IDX_GRID(i1, j1, k1)];
        snapGrid.weights[IDX_GRID(i0, j0, k0)] = weight;
    } else if (snapGrid.snapBack[IDX_GRID(i0, j0, k0)] && snapGrid.weights[IDX_GRID(i0, j0, k0)] > weight) {
        // If the isoSurface value on the already snapped edge is further away than on the current one.
        snapGrid.snapBackTo[IDX_GRID(i0, j0, k0)] = gridPoints[IDX_GRID(i1, j1, k1)];
        snapGrid.snapBackToNormals[IDX_GRID(i0, j0, k0)] = gridNormals[IDX_GRID(i1, j1, k1)];
        snapGrid.weights[IDX_GRID(i0, j0, k0)] = weight;
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
    float epsilon = 0.00001;

    // Weight or distance.
    float weight0 = gamma;
    float weight1 = gamma;

    if ((cartesianGrid[IDX_GRID(i0, j0, k0)] < isoLevel && cartesianGrid[IDX_GRID(i1, j1, k1)] > isoLevel)
            || (cartesianGrid[IDX_GRID(i0, j0, k0)] > isoLevel && cartesianGrid[IDX_GRID(i1, j1, k1)] < isoLevel)) {
        // If this is a +/- edge.
        float value_difference = cartesianGrid[IDX_GRID(i1, j1, k1)] - cartesianGrid[IDX_GRID(i0, j0, k0)];

        if (value_difference > epsilon || value_difference < -epsilon) {
            weight0 = (cartesianGrid[IDX_GRID(i1, j1, k1)] - isoLevel) / value_difference;
            weight1 = (isoLevel - cartesianGrid[IDX_GRID(i0, j0, k0)]) / value_difference;
        }
        else {
            weight0 = 0.5;
            weight1 = 0.5;
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
    snapGrid.gridValues = new float[nx * ny * ny];
    memcpy(snapGrid.gridValues, cartesianGrid, sizeof(float) * nx * ny * ny);

    // Go over all vertices of the grid (just not the last ones in the respective direction since we want to got over
    // the edges).
    for (int i = 0; i < nz - 1; i++) {
        for (int j = 0; j < ny - 1; j++) {
            for (int k = 0; k < nx - 1; k++) {
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
    return snapGrid;
}

void polygonizeSnapMC(
        const float* voxelGrid, int nx, int ny, int nz, float isoLevel, const float gamma,
        std::vector<glm::vec3>& vertexPositions, std::vector<glm::vec3>& vertexNormals) {
    int numCellsX = nx - 1;
    int numCellsY = ny - 1;
    int numCellsZ = nz - 1;
    int numCells = numCellsX * numCellsY * numCellsZ;
    GridCell* gridCells = new GridCell[numCells];

    glm::vec3* gridPoints = new glm::vec3[nx * ny * nz];
    glm::vec3* gridNormals = new glm::vec3[nx * ny * nz];
    #pragma omp parallel for default(none) shared(voxelGrid, gridPoints, gridNormals, nx, ny, nz)
    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                // Compute the normal vector.
                glm::vec3 h(1.0f / float(nx), 1.0f / float(ny), 1.0f / float(nz));
                float normalX = (voxelGrid[IDX_GRID(std::min(x+1, nx-1), y, z)]
                                 - voxelGrid[IDX_GRID(std::max(x-1, 0), y, z)]) / (-2.0f*h[0]);
                float normalY = (voxelGrid[IDX_GRID(x, std::min(y+1, ny-1), z)]
                                 - voxelGrid[IDX_GRID(x, std::max(y-1, 0), z)]) / (-2.0f*h[1]);
                float normalZ = (voxelGrid[IDX_GRID(x, y, std::min(z+1, nz-1))]
                                 - voxelGrid[IDX_GRID(x, y, std::max(z-1, 0))]) / (-2.0f*h[2]);
                glm::vec3 n = glm::normalize(glm::vec3(normalX, normalY, normalZ));

                int idx = IDX_GRID(x, y, z);
                gridPoints[idx] = glm::vec3{float(x), float(y), float(z)};
                gridNormals[idx] = glm::vec3{n[0], n[1], n[2]};
            }
        }
    }

    #pragma omp parallel for default(none) shared(voxelGrid, gridCells, numCellsX, numCellsY, numCellsZ, nx, ny, nz)
    for (int z = 0; z < numCellsZ; z++) {
        for (int y = 0; y < numCellsY; y++) {
            for (int x = 0; x < numCellsX; x++) {
                GridCell& gridCell = gridCells[x + (y + z * numCellsY) * numCellsX];

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
                    glm::vec3 h(1.0f / float(nx), 1.0f / float(ny), 1.0f / float(nz));
                    float normalX = (voxelGrid[IDX_GRID(std::min(gridIndex[0]+1, numCellsX), gridIndex[1], gridIndex[2])]
                                     - voxelGrid[IDX_GRID(std::max(gridIndex[0]-1, 0), gridIndex[1], gridIndex[2])]) / (-2.0f*h[0]);
                    float normalY = (voxelGrid[IDX_GRID(gridIndex[0], std::min(gridIndex[1]+1, numCellsY), gridIndex[2])]
                                     - voxelGrid[IDX_GRID(gridIndex[0], std::max(gridIndex[1]-1, 0), gridIndex[2])]) / (-2.0f*h[1]);
                    float normalZ = (voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], std::min(gridIndex[2]+1, numCellsZ))]
                                     - voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], std::max(gridIndex[2]-1, 0))]) / (-2.0f*h[2]);
                    glm::vec3 n = glm::normalize(glm::vec3(normalX, normalY, normalZ));

                    gridCell.v[l] = glm::vec3{float(gridIndex[0]), float(gridIndex[1]), float(gridIndex[2])};
                    gridCell.n[l] = glm::vec3{n[0], n[1], n[2]};
                    gridCell.f[l] = voxelGrid[IDX_GRID(gridIndex[0], gridIndex[1], gridIndex[2])];
                }
            }
        }
    }

    SnapGrid snapGrid;
    snapGrid.gridValues = new float[nx * ny * nz];
    snapGrid.snapBack = new bool[nx * ny * nz];
    snapGrid.snapBackTo = new glm::vec3[nx * ny * nz];
    snapGrid.snapBackToNormals = new glm::vec3[nx * ny * nz];
    snapGrid.weights = new float[nx * ny * nz];
    memset(snapGrid.snapBack, 0, sizeof(bool) * nx * ny * nz);

    constructCartesianSnapGridScalarField(snapGrid, voxelGrid, gridPoints, gridNormals, isoLevel, gamma, nx, ny, nz);
    for (int z = 0; z < numCellsZ; z++) {
        for (int y = 0; y < numCellsY; y++) {
            for (int x = 0; x < numCellsX; x++) {
                GridCell& gridCell = gridCells[x + (y + z * numCellsY) * numCellsX];
                polygonizeSnapMC(
                        gridCell, isoLevel, snapGrid, nx, ny, nz, x, y, z, vertexPositions, vertexNormals);
            }
        }
    }

    delete[] snapGrid.gridValues;
    delete[] snapGrid.snapBack;
    delete[] snapGrid.snapBackTo;
    delete[] snapGrid.snapBackToNormals;
    delete[] snapGrid.weights;

    delete[] gridPoints;
    delete[] gridNormals;
    delete[] gridCells;
}
