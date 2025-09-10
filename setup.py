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

import sys
import os
import glob
import shutil
from pathlib import Path
from urllib.request import urlopen
import setuptools
from setuptools import setup, Extension, find_packages
from setuptools.command.egg_info import egg_info
from setuptools.dist import Distribution
from setuptools.command import bdist_egg
from setuptools.command.build_ext import build_ext


extra_compile_args = []
extra_link_args = []
if os.name == 'nt':
    extra_compile_args.append('/std:c++17')
    extra_compile_args.append('/Zc:__cplusplus')
    extra_compile_args.append('/openmp')
    extra_link_args.append('/openmp')
else:
    extra_compile_args.append('-std=c++17')
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
    'pymodule',
    'pymodule/sgl',
]
source_files = []
source_files += find_all_sources_in_dir('src')
source_files += find_all_sources_in_dir('pymodule')


data_files_all = []
data_files = ['pymodule/isosurfacecpp.pyi']
libraries = []
defines = [
    ('USE_GLM',),
    ('DLL_OBJECT', ''),
]

data_files_all.append(('.', data_files))


def update_data_files_recursive(data_files_all, directory):
    files_in_directory = []
    for filename in os.listdir(directory):
        abs_file = directory + "/" + filename
        if os.path.isdir(abs_file):
            update_data_files_recursive(data_files_all, abs_file)
        else:
            files_in_directory.append(abs_file)
    if len(files_in_directory) > 0:
        data_files_all.append((directory, files_in_directory))


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

uses_pip = \
    ('_' in os.environ and (os.environ['_'].endswith('pip') or os.environ['_'].endswith('pip3'))) \
    or 'PIP_BUILD_TRACKER' in os.environ
if os.path.exists('isosurfacecpp'):
    shutil.rmtree('isosurfacecpp')
if uses_pip:
    Path('isosurfacecpp/Data').mkdir(parents=True, exist_ok=True)
    shutil.copy('pymodule/isosurfacecpp.pyi', 'isosurfacecpp/__init__.pyi')
    shutil.copy('LICENSE', 'isosurfacecpp/LICENSE')
    pkg_data = ['**/LICENSE']
    ext_modules = [
        Extension(
            'isosurfacecpp.isosurfacecpp',
            source_files,
            libraries=libraries,
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args
        )
    ]
    dist = Distribution(attrs={'name': 'isosurfacecpp', 'version': '0.0.0', 'ext_modules': ext_modules})
    bdist_egg_cmd = dist.get_command_obj('bdist_egg')
    build_cmd = bdist_egg_cmd.get_finalized_command('build_ext')
    isosurfacecpp_so_file = ''
    for ext in build_cmd.extensions:
        fullname = build_cmd.get_ext_fullname(ext.name)
        filename = build_cmd.get_ext_filename(fullname)
        isosurfacecpp_so_file = os.path.basename(filename)
    with open('isosurfacecpp/__init__.py', 'w') as file:
        file.write('import numpy\n\n')
        file.write('def __bootstrap__():\n')
        file.write('    global __bootstrap__, __loader__, __file__\n')
        file.write('    import sys, pkg_resources, importlib.util\n')
        file.write(f'    __file__ = pkg_resources.resource_filename(__name__, \'{isosurfacecpp_so_file}\')\n')
        file.write('    __loader__ = None; del __bootstrap__, __loader__\n')
        file.write('    spec = importlib.util.spec_from_file_location(__name__,__file__)\n')
        file.write('    mod = importlib.util.module_from_spec(spec)\n')
        file.write('    spec.loader.exec_module(mod)\n')
        file.write('__bootstrap__()\n')
    setup(
        name='IsosurfaceCpp',
        author='Christoph Neuhauser',
        ext_modules=ext_modules,
        packages=find_packages(include=['isosurfacecpp', 'isosurfacecpp.*']),
        package_data={'isosurfacecpp': ['**/*.py', '**/*.pyi', '**/*.md', '**/*.txt', '**/*.xml', '**/*.glsl'] + pkg_data},
        #include_package_data=True,
        cmdclass={
            'build_ext': build_ext,
            'egg_info': EggInfoInstallLicense
        },
        license_files=('LICENSE',),
        include_dirs=include_dirs
    )
else:
    setup(
        name='IsosurfaceCpp',
        author='Christoph Neuhauser',
        ext_modules=[
            Extension(
                'isosurfacecpp',
                source_files,
                libraries=libraries,
                extra_compile_args=extra_compile_args,
                extra_link_args=extra_link_args
            )
        ],
        data_files=data_files_all,
        cmdclass={
            'build_ext': build_ext,
            'egg_info': EggInfoInstallLicense
        },
        license_files=('LICENSE',),
        include_dirs=include_dirs
    )
