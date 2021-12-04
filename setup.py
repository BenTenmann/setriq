from __future__ import print_function
import sys
if sys.version_info < (3,):
    print('Python 2 has reached end-of-life and is not supported by setriq.')
    sys.exit(-1)
if sys.platform == 'win32' and sys.maxsize.bit_length() == 31:
    print('32-bit Windows Python runtime is not supported. Please switch to 64-bit Python.')
    sys.exit(-1)

import logging
import pathlib
import platform
import re
import traceback
from glob import glob

from pybind11.setup_helpers import ParallelCompile, Pybind11Extension
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext as base_build_ext

python_min_version = (3, 7, 2)
python_min_version_str = '.'.join(map(str, python_min_version))
if sys.version_info < python_min_version:
    print('Minimum Python version required is {}. Current version is {}'.format(python_min_version_str,
                                                                                platform.python_version()))

logging.basicConfig(
    level=logging.DEBUG,
    format='%(name)s - %(asctime)s.%(msecs)03d - %(levelname)s - %(module)s.%(funcName)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# find the sem-rel version
dirname = pathlib.Path(__file__).parent
long_description = (dirname / 'README.md').read_text()

try:
    changelog = (dirname / 'CHANGELOG.md').read_text()
    __version__, *_ = re.findall(r"\[([0-9.]+)]", changelog)
except (FileNotFoundError, ValueError) as ex:
    __version__ = '0.1.0'
    logging.error(ex)
    logging.error(traceback.print_exc())
    logging.warning(f'Unable to get semantic release version. Setting version to {__version__}.')

SOURCE_DIR = 'src'
PROJECT_NAME = 'setriq'

ParallelCompile('NPY_NUM_BUILD_JOBS').install()


def get_compile_args():
    if platform.system() == 'Darwin':
        return ['-Xpreprocessor', '-fopenmp']

    return ['-fopenmp']


def get_link_args():
    if platform.system() == 'Darwin':
        return ['-lomp']

    return ['-fopenmp']


extensions = [
    Pybind11Extension(
        f'{PROJECT_NAME}._C',
        sources=sorted(glob(f'{SOURCE_DIR}/{PROJECT_NAME}/_C/**/*.cpp', recursive=True)),
        cxx_std=14,
        define_macros=[('VERSION_INFO', __version__)],
        include_dirs=['include/setriq'],
        extra_compile_args=get_compile_args(),
        extra_link_args=get_link_args()
    ),
]

setup(
    name=PROJECT_NAME,
    version=__version__,
    description='Python package written in C++ for pairwise distance computation for sequences.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='benjamin-tenmann',
    author_email='b.tenmann@me.com',
    url='https://github.com/BenTenmann/setriq',
    ext_modules=extensions,
    license='MIT',
    python_requires='>=3.7,<3.10',
    install_requires=[
        'glom>=20.0.0,<21.0.0',
        'numpy>=1.0.0,<2.0.0',
        'srsly>=2.0.0,<3.0.0',
    ],
    package_dir={f'{PROJECT_NAME}': f'{SOURCE_DIR}/{PROJECT_NAME}'},
    packages=find_packages(where=f'{SOURCE_DIR}', exclude=['tests', 'scripts']),
    package_data={f'{PROJECT_NAME}': ['data/*.json']},
    include_package_data=True
)
