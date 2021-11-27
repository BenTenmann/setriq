import logging
import os
import pathlib
import re
import traceback
from glob import glob

from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup, find_packages

logging.basicConfig(
    level=logging.DEBUG,
    format='%(name)s - %(asctime)s.%(msecs)03d - %(levelname)s - %(module)s.%(funcName)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

url = os.environ.get('URL', '')

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

extensions = [
    Pybind11Extension(
        'setriq',
        sources=sorted(glob('src/**/*.cpp', recursive=True)),
        cxx_std=14,
        define_macros=[('VERSION_INFO', __version__)],
        include_dirs=['include/setriq'],
        extra_compile_args=['-O3']
    ),
]

setup(
    name='setriq',
    version=__version__,
    description='Python package written in C++ for pairwise distance computation for sequences.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='benjamin-tenmann',
    author_email='b.tenmann@me.com',
    url=url,
    ext_modules=extensions,
    license='MIT',
    requires=[],
    python_requires='>=3.7,<3.10',
    packages=find_packages(exclude=['tests']),
    data_files=[('data', ['data/blosum-62.json'])]
)
