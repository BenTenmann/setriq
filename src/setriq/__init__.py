"""
setriq
======

Provides:
    1. fast computation of pairwise (immunoglobulin) sequence distances
    2. parsed BLOSUM substitution matrices

About
-----
This package is a one-trick-pony which does its trick pretty well: efficient computation of pairwise sequence distances
on CPU.

Examples
--------
>>> import setriq
>>> metric = setriq.CdrDist()
>>>
>>> sequences = [
...     'CASSLKPNTEAFF',
...     'CASSAHIANYGYTF',
...     'CASRGATETQYF'
... ]
>>> distances = metric(sequences)

References
----------
[1] Dash, P., Fiore-Gartland, A.J., Hertz, T., Wang, G.C., Sharma, S., Souquette, A., Crawford, J.C., Clemens, E.B.,
    Nguyen, T.H., Kedzierska, K. and La Gruta, N.L., 2017. Quantifiable predictive features define epitope-specific
    T cell receptor repertoires. Nature, 547(7661), pp.89-93. (https://doi.org/10.1038/nature22383)

[2] Levenshtein, V.I., 1966, February. Binary codes capable of correcting deletions, insertions, and reversals. In
    Soviet physics doklady (Vol. 10, No. 8, pp. 707-710).

[3] python-Levenshtein (https://github.com/ztane/python-Levenshtein)

[4] Thakkar, N. and Bailey-Kellogg, C., 2019. Balancing sensitivity and specificity in distinguishing TCR groups by CDR
    sequence similarity. BMC bioinformatics, 20(1), pp.1-14. (https://doi.org/10.1186/s12859-019-2864-8)
"""

from .modules import (
    CdrDist, Levenshtein, TcrDist, Hamming, Jaro, JaroWinkler,
    SubstitutionMatrix, BLOSUM45, BLOSUM62, BLOSUM90,
    single_dispatch
)
