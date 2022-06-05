"""
modules
=======

This sub-package contains all of the `setriq` Python API.

Examples
--------
>>> from setriq import modules
>>> sub_mat = modules.BLOSUM62
>>> sub_mat('A', 'L')
... -1

"""

from .substitution import (
    SubstitutionMatrix,
    BLOSUM45,
    BLOSUM62,
    BLOSUM90
)
from .distances import (
    CdrDist,
    Levenshtein,
    TcrDist,
    Hamming,
    Jaro,
    JaroWinkler
)

from . import single_dispatch

__all__ = [
    'SubstitutionMatrix',
    'BLOSUM45',
    'BLOSUM62',
    'BLOSUM90',
    'CdrDist',
    'Levenshtein',
    'TcrDist',
    'Hamming',
    'Jaro',
    'JaroWinkler',
    'single_dispatch'
]
