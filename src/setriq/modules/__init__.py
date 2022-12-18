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

from . import single_dispatch
from .distances import (
    CdrDist,
    Hamming,
    Jaro,
    JaroWinkler,
    Levenshtein,
    LongestCommonSubstring,
    OptimalStringAlignment,
    TcrDist,
)
from .substitution import BLOSUM45, BLOSUM62, BLOSUM90, SubstitutionMatrix

__all__ = [
    "SubstitutionMatrix",
    "BLOSUM45",
    "BLOSUM62",
    "BLOSUM90",
    "CdrDist",
    "Levenshtein",
    "TcrDist",
    "Hamming",
    "Jaro",
    "JaroWinkler",
    "LongestCommonSubstring",
    "OptimalStringAlignment",
    "single_dispatch",
]
