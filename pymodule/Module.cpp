/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2025, Christoph Neuhauser
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
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <glm/glm.hpp>
#include <Utils/Mesh/IndexMesh.hpp>
#include <MarchingCubes.hpp>
#include <SnapMC.hpp>

pybind11::tuple pyPolygonizeMC(
        pybind11::array_t<float, pybind11::array::c_style | pybind11::array::forcecast> gridData,
        float dx, float dy, float dz, float isoLevel) {
    if (gridData.ndim() != 3) {
        throw std::runtime_error("Number of dimensions must be 3.");
    }
    auto* gridDataPtr = gridData.data();
    auto nx = int(gridData.shape(2));
    auto ny = int(gridData.shape(1));
    auto nz = int(gridData.shape(0));
    std::vector<glm::vec3> isosurfaceVertexPositions;
    std::vector<glm::vec3> isosurfaceVertexNormals;
    polygonizeMarchingCubes(
            gridDataPtr, nx, ny, nz, dx, dy, dz, isoLevel,
            isosurfaceVertexPositions, isosurfaceVertexNormals);
    std::vector<uint32_t> triangleIndices;
    std::vector<glm::vec3> vertexPositions;
    std::vector<glm::vec3> vertexNormals;

    float step = std::min(dx, std::min(dy, dz));
    sgl::computeSharedIndexRepresentation(
            isosurfaceVertexPositions, isosurfaceVertexNormals,
            triangleIndices, vertexPositions, vertexNormals,
            1e-5f * step);

    size_t numTriangles = triangleIndices.size() / 3;
    pybind11::array_t<uint32_t> pyTriangleIndices(std::vector<size_t>{numTriangles, 3});
    pybind11::array_t<float> pyVertexPositions(std::vector<size_t>{vertexPositions.size(), 3});
    pybind11::array_t<float> pyVertexNormals(std::vector<size_t>{vertexNormals.size(), 3});
    memcpy(pyTriangleIndices.mutable_data(), triangleIndices.data(), triangleIndices.size() * sizeof(uint32_t));
    memcpy(pyVertexPositions.mutable_data(), vertexPositions.data(), vertexPositions.size() * sizeof(glm::vec3));
    memcpy(pyVertexNormals.mutable_data(), vertexNormals.data(), vertexNormals.size() * sizeof(glm::vec3));
    return pybind11::make_tuple(pyTriangleIndices, pyVertexPositions, pyVertexNormals);
}

pybind11::tuple pyPolygonizeSnapMC(
        pybind11::array_t<float, pybind11::array::c_style | pybind11::array::forcecast> gridData,
        float dx, float dy, float dz, float isoLevel, float gamma) {
    if (gridData.ndim() != 3) {
        throw std::runtime_error("Number of dimensions must be 3.");
    }
    auto* gridDataPtr = gridData.data();
    auto nx = int(gridData.shape(2));
    auto ny = int(gridData.shape(1));
    auto nz = int(gridData.shape(0));
    std::vector<glm::vec3> isosurfaceVertexPositions;
    std::vector<glm::vec3> isosurfaceVertexNormals;
    polygonizeSnapMC(
            gridDataPtr, nx, ny, nz, dx, dy, dz, isoLevel, gamma,
            isosurfaceVertexPositions, isosurfaceVertexNormals);
    std::vector<uint32_t> triangleIndices;
    std::vector<glm::vec3> vertexPositions;
    std::vector<glm::vec3> vertexNormals;

    float step = std::min(dx, std::min(dy, dz));
    sgl::computeSharedIndexRepresentation(
            isosurfaceVertexPositions, isosurfaceVertexNormals,
            triangleIndices, vertexPositions, vertexNormals,
            1e-5f * step);

    size_t numTriangles = triangleIndices.size() / 3;
    pybind11::array_t<uint32_t> pyTriangleIndices(std::vector<size_t>{numTriangles, 3});
    pybind11::array_t<float> pyVertexPositions(std::vector<size_t>{vertexPositions.size(), 3});
    pybind11::array_t<float> pyVertexNormals(std::vector<size_t>{vertexNormals.size(), 3});
    memcpy(pyTriangleIndices.mutable_data(), triangleIndices.data(), triangleIndices.size() * sizeof(uint32_t));
    memcpy(pyVertexPositions.mutable_data(), vertexPositions.data(), vertexPositions.size() * sizeof(glm::vec3));
    memcpy(pyVertexNormals.mutable_data(), vertexNormals.data(), vertexNormals.size() * sizeof(glm::vec3));
    return pybind11::make_tuple(pyTriangleIndices, pyVertexPositions, pyVertexNormals);
}

PYBIND11_MODULE(isosurfacecpp, m) {
    m.def("polygonize_mc", pyPolygonizeMC,
          "Polygonizes the passed grid using Marching Cubes.",
          pybind11::arg("grid_data"), pybind11::arg("dx"), pybind11::arg("dy"), pybind11::arg("dz"),
          pybind11::arg("iso_level"));
    m.def("polygonize_mc", pyPolygonizeSnapMC,
          "Polygonizes the passed grid using SnapMC.",
          pybind11::arg("grid_data"), pybind11::arg("dx"), pybind11::arg("dy"), pybind11::arg("dz"),
          pybind11::arg("iso_level"), pybind11::arg("gamma"));
}
