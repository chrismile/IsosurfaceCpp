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

import os
import setuptools
from setuptools import setup
from setuptools.command.egg_info import egg_info
from pybind11.setup_helpers import Pybind11Extension, build_ext


extra_compile_args = []
extra_link_args = []
if os.name == 'nt':
    extra_compile_args.append('/Zc:__cplusplus')
    extra_compile_args.append('/openmp')
    extra_link_args.append('/openmp')
else:
    extra_compile_args.append('-fopenmp')
    extra_link_args.append('-fopenmp')


class EggInfoInstallLicense(egg_info):
    def run(self):
        if not self.distribution.have_run.get('install', True):
            self.mkpath(self.egg_info)
            self.copy_file('LICENSE', self.egg_info)
        egg_info.run(self)


def find_all_sources_in_dir(root_dir):
    source_files = []
    for root, subdirs, files in os.walk(root_dir):
        for filename in files:
            if filename.endswith('.cpp') or filename.endswith('.cc') or filename.endswith('.c'):
                source_files.append(root + "/" + filename)
    return source_files


include_dirs = [
    'src',
    'isosurfacecpp',
    'isosurfacecpp/sgl',
]
source_files = []
source_files += find_all_sources_in_dir('src')
source_files += find_all_sources_in_dir('isosurfacecpp')

libraries = []
defines = [
    ('USE_GLM',),
    ('DLL_OBJECT', ''),
]

for define in defines:
    if os.name == 'nt':
        if len(define) == 1:
            extra_compile_args.append('/D')
            extra_compile_args.append(f'{define[0]}')
        else:
            extra_compile_args.append('/D')
            extra_compile_args.append(f'{define[0]}={define[1]}')
    else:
        if len(define) == 1:
            extra_compile_args.append(f'-D{define[0]}')
        else:
            extra_compile_args.append(f'-D{define[0]}={define[1]}')

setup(
    name='IsosurfaceCpp',
    author='Christoph Neuhauser',
    ext_modules=[
        Pybind11Extension(
            'isosurfacecpp',
            source_files,
            cxx_std=17,
            libraries=libraries,
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args
        )
    ],
    packages=['isosurfacecpp', 'isosurfacecpp-stubs'],
    package_data={ 'isosurfacecpp-stubs': ['*.pyi'] },
    cmdclass={
        'build_ext': build_ext,
        'egg_info': EggInfoInstallLicense
    },
    license_files=('LICENSE',),
    include_dirs=include_dirs
)
