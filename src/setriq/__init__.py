"""
setriq
======

Provides:
    1. fast computation of pairwise (immunoglobulin) sequence distances
    2. convenience methods for handling immunoglobulin sequences
    3. parsed BLOSUM substitution matrices

How to use the documentation
----------------------------

"""

from .modules import (
    CdrDist, TcrDist,
    SubstitutionMatrix, BLOSUM45, BLOSUM62, BLOSUM90
)
