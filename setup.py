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
import shutil
import subprocess
import traceback
from glob import glob

from pybind11.setup_helpers import ParallelCompile, Pybind11Extension
from setuptools import setup, find_packages

python_min_version = (3, 7, 2)
python_min_version_str = '.'.join(map(str, python_min_version))
if sys.version_info < python_min_version:
    logging.error('Minimum Python version required is {}. Current version is {}'.format(python_min_version_str,
                                                                                        platform.python_version()))
    sys.exit(-1)

logging.basicConfig(
    level=logging.DEBUG,
    format='%(name)s - %(asctime)s.%(msecs)03d - %(levelname)s - %(module)s.%(funcName)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# find the sem-rel version
DIR = pathlib.Path(__file__).parent

SOURCE_DIR = 'src'
PROJECT_NAME = 'setriq'
EXTENSION_NAME = '_C'

ParallelCompile('NPY_NUM_BUILD_JOBS').install()


class BuildFlags:
    _tools_formulae = {
        'Darwin': ('brew', 'libomp'),
    }
    _args = {
        'Darwin': {
            'compiler': ['-Xpreprocessor', '-fopenmp'],
            'linker': ['-lomp']
        }
    }
    _default_args = {
        'compiler': ['-fopenmp'],
        'linker': ['-fopenmp']
    }
    compiler: list
    linker: list

    def __init__(self):
        self._system = platform.system()
        tool, formula = self._tools_formulae.get(self._system, ('apt', 'libomp-dev'))
        args = self._args.get(self._system, self._default_args)

        not_found = self._libomp_check(tool, formula)
        if not_found is not None:
            logging.warning(f'{repr(not_found)} not found -- cannot compile parallelized code')
            logging.info('for information on how to enable CPU parallelization, '
                         'please see https://github.com/BenTenmann/setriq#requirements')
            for key in args:
                args[key] = []

        for key, val in args.items():
            self.__setattr__(key, val)

    @staticmethod
    def _libomp_check(tool, formula):
        if shutil.which(tool) is None:
            return tool

        formulae = subprocess.check_output([tool, 'list']).decode()
        if formula not in formulae:
            return formula

        return None


def main():
    long_description = (DIR / 'README.md').read_text()

    try:
        changelog = (DIR / 'CHANGELOG.md').read_text()
        __version__, *_ = re.findall(r"\[([0-9.]+)]", changelog)
    except (FileNotFoundError, ValueError) as ex:
        __version__ = '0.1.0'
        logging.error(ex)
        logging.error(traceback.print_exc())
        logging.warning(f'Unable to get semantic release version. Setting version to {__version__}.')

    flags = BuildFlags()
    extensions = [
        Pybind11Extension(
            f'{PROJECT_NAME}.{EXTENSION_NAME}',
            sources=sorted(glob(f'{SOURCE_DIR}/{PROJECT_NAME}/{EXTENSION_NAME}/**/*.cpp', recursive=True)),
            cxx_std=14,
            define_macros=[
                ('VERSION_INFO', __version__),
                ('EXTENSION_NAME', EXTENSION_NAME)
            ],
            include_dirs=['include/setriq'],
            extra_compile_args=flags.compiler,
            extra_link_args=flags.linker
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
            'pandas>=1.0.0,<2.0.0',
            'srsly>=2.0.0,<3.0.0',
        ],
        package_dir={f'{PROJECT_NAME}': f'{SOURCE_DIR}/{PROJECT_NAME}'},
        packages=find_packages(where=f'{SOURCE_DIR}', exclude=['tests', 'scripts']),
        package_data={f'{PROJECT_NAME}': ['data/*.json']},
        include_package_data=True,
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            'Topic :: Software Development :: Build Tools',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
        ]
    )


if __name__ == '__main__':
    main()
