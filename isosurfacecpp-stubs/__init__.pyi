# BSD 2-Clause License
#
# Copyright (c) 2025, Christoph Neuhauser
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import annotations
import numpy as np
import isosurfacecpp
import typing

__all__ = [
    "polygonize_mc",
    "polygonize_snapmc",
]


def polygonize_mc(
        grid_data: np.ndarray, dx: float, dy: float, dz: float, iso_level: float
    ) -> np.ndarray: -> typing.Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Polygonizes the passed grid using Marching Cubes.
    """

def polygonize_snapmc(
        grid_data: np.ndarray, dx: float, dy: float, dz: float, iso_level: float, gamma: float
    ) -> np.ndarray: -> typing.Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Polygonizes the passed grid using SnapMC.
    """
