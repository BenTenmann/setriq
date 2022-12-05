r"""
single_dispatch
===============

Functional implementations of the normal Metric classes, with the difference that the comparison happens only between
two strings. This can be useful for integration with other tools such as PySpark, where we want to access the fast
distance function implementations, without the forced pairwise comparisons.

Examples
--------

An example for computing the pairwise sequence distances using PySpark and setriq:

>>> from pyspark.sql import SparkSession
>>> from pyspark.sql.functions import udf
>>> from pyspark.sql.types import DoubleType
>>>
>>> spark = SparkSession \
...    .builder \
...    .master("local[2]") \
...    .appName("setriq-spark") \
...    .getOrCreate()
>>>
>>> df = spark.createDataFrame([('CASSLKPNTEAFF',), ('CASSAHIANYGYTF',), ('CASRGATETQYF',)], ['sequence'])
>>> df = df.withColumnRenamed('sequence', 'a').crossJoin(df.withColumnRenamed('sequence', 'b'))
>>>
>>> lev_udf = udf(levenshtein, returnType=DoubleType())  # single dispatch levenshtein distance
>>> df = df.withColumn('distance', lev_udf('a', 'b'))
>>> df.show()
+--------------+--------------+--------+
|             a|             b|distance|
+--------------+--------------+--------+
| CASSLKPNTEAFF| CASSLKPNTEAFF|     0.0|
| CASSLKPNTEAFF|CASSAHIANYGYTF|     8.0|
| CASSLKPNTEAFF|  CASRGATETQYF|     8.0|
|CASSAHIANYGYTF| CASSLKPNTEAFF|     8.0|
|CASSAHIANYGYTF|CASSAHIANYGYTF|     0.0|
|CASSAHIANYGYTF|  CASRGATETQYF|     9.0|
|  CASRGATETQYF| CASSLKPNTEAFF|     8.0|
|  CASRGATETQYF|CASSAHIANYGYTF|     9.0|
|  CASRGATETQYF|  CASRGATETQYF|     0.0|
+--------------+--------------+--------+

"""

from typing import List, Optional

import setriq._C as C

from .substitution import BLOSUM45, SubstitutionMatrix
from .utils import (
    check_jaro_weights,
    check_jaro_winkler_params,
    ensure_equal_sequence_length_sd,
    single_dispatch,
    tcr_dist_sd_component_check,
)

__all__ = [
    "cdr_dist",
    "levenshtein",
    "tcr_dist",
    "hamming",
    "jaro",
    "jaro_winkler",
    "longest_common_substring",
    "optimal_string_alignment",
]


@single_dispatch
def cdr_dist(
    a: str,
    b: str,
    substitution_matrix: SubstitutionMatrix = BLOSUM45,
    gap_opening_penalty: float = 10.0,
    gap_extension_penalty: float = 1.0,
) -> float:
    """
    Compute the CDRdist [1]_ metric between two sequences.

    {params}
    substitution_matrix: SubstitutionMatrix
        A substitution matrix object to inform the alignment scoring.
    gap_opening_penalty: float
        The penalty given to an alignment based on a gap opening. Values other than 0 give affine gap scoring.
        (default=10.0)
    gap_extension_penalty: float
        The penalty used to score the extension of the gap after opening. (default=1.0)

    {returns}

    Examples
    --------
    >>> seq = ('AASQ', 'PASQ')
    >>> cdr_dist(*seq)  # default params
    >>> cdr_dist(*seq, substitution_matrix=BLOSUM45,  # custom params
    ...          gap_opening_penalty=5.0, gap_extension_penalty=2.0)

    References
    ----------
    .. [1] Thakkar, N. and Bailey-Kellogg, C., 2019. Balancing sensitivity and specificity in distinguishing TCR groups
       by CDR sequence similarity. BMC bioinformatics, 20(1), pp.1-14. (https://doi.org/10.1186/s12859-019-2864-8)

    """
    distance = C.cdr_dist_sd(
        a,
        b,
        **substitution_matrix,  # type: ignore[arg-type]
        gap_opening_penalty=gap_opening_penalty,
        gap_extension_penalty=gap_extension_penalty
    )
    return distance


@single_dispatch
def levenshtein(a: str, b: str, extra_cost: float = 0.0) -> float:
    """
    Compute the Levenshtein distance [1]_ between two sequences. Based on the implementation in [2]_.

    {params}
    extra_cost: float
        Additional cost assigned by Levenshtein algorithm.

    {returns}

    Examples
    --------
    >>> levenshtein('AASQ', 'PASQ')

    References
    ----------
    .. [1] Levenshtein, V.I., 1966, February. Binary codes capable of correcting deletions, insertions, and reversals.
       In Soviet physics doklady (Vol. 10, No. 8, pp. 707-710).
    .. [2] python-Levenshtein (https://github.com/ztane/python-Levenshtein)

    """
    distance = C.levenshtein_sd(a, b, extra_cost=extra_cost)
    return distance


@single_dispatch
@ensure_equal_sequence_length_sd
def tcr_dist_component(
    a: str,
    b: str,
    substitution_matrix: SubstitutionMatrix,
    gap_penalty: float,
    gap_symbol: str = "-",
    weight: float = 1.0,
) -> float:
    distance = C.tcr_dist_component_sd(
        a,
        b,
        **substitution_matrix,  # type: ignore[arg-type]
        gap_penalty=gap_penalty,
        gap_symbol=gap_symbol,
        weight=weight
    )
    return distance


@tcr_dist_sd_component_check
def tcr_dist(a: dict, b: dict, **component_def) -> float:
    """
    Compute the TCRdist [1]_ metric between two sequences.

    Parameters
    ----------
    a: dict

    b: dict
    component_def

    Returns
    -------
    distance: float
        The distance between the two sequences.

    References
    ----------
    .. [1] Dash, P., Fiore-Gartland, A.J., Hertz, T., Wang, G.C., Sharma, S., Souquette, A., Crawford, J.C., Clemens,
       E.B., Nguyen, T.H., Kedzierska, K. and La Gruta, N.L., 2017. Quantifiable predictive features define
       epitope-specific T cell receptor repertoires. Nature, 547(7661), pp.89-93. (https://doi.org/10.1038/nature22383)

    """
    distance = 0.0
    for name, component in component_def.items():
        distance += tcr_dist_component(a[name], b[name], **component_def[name])
    return distance


@single_dispatch
@ensure_equal_sequence_length_sd
def hamming(a: str, b: str, mismatch_score: float = 1.0) -> float:
    """
    Compute the Hamming [1]_ distance between two sequences.

    {params}
    mismatch_score: float
        The weight given to a mismatch between the two sequence positions.

    {returns}

    Examples
    --------
    >>> hamming('AASQ', 'PASQ')
    >>> hamming('AASQ', 'PAS')  # error! different length sequences

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Hamming_distance

    """
    distance = C.hamming_sd(a, b, mismatch_score=mismatch_score)
    return distance


@single_dispatch
def jaro(a: str, b: str, jaro_weights: Optional[List[float]] = None) -> float:
    """
    Compute the Jaro [1] distance between two sequences. Adapted from [2].

    {params}
    jaro_weights: List[float]

    {returns}

    Examples
    --------
    >>> jaro('AASQ', 'PASQ')

    References
    ----------
    .. [1] Jaro, M.A., 1989. Advances in record-linkage methodology as applied to matching the 1985 census of Tampa,
       Florida. Journal of the American Statistical Association, 84(406), pp.414-420.
    .. [2] Van der Loo, M.P., 2014. The stringdist package for approximate string matching. R J., 6(1), p.111.

    """
    jaro_weights = check_jaro_weights(jaro_weights)
    distance = C.jaro_sd(a, b, jaro_weights=jaro_weights)
    return distance


@single_dispatch
@check_jaro_winkler_params
def jaro_winkler(
    a: str, b: str, p: float, max_l: int = 4, jaro_weights: Optional[List[float]] = None
) -> float:
    """
    Compute the Jaro-Winkler [1]_ distance between two sequences.

    {params}
    p: float
        The scaling factor applied to the common prefix re-weighting. The value needs to be in the range [0.0, 0.25]. If
        set to 0.0, Jaro-Winkler reduces down to Jaro.
    max_l: int
        The maximum length common prefix. (``default=4``)
    jaro_weights: List[float]

    {returns}

    Examples
    --------
    >>> jaro_winkler('AASQ', 'PASQ', p=0.10)

    References
    ----------
    .. [1] Winkler, W.E., 1990. String comparator metrics and enhanced decision rules in the Fellegi-Sunter model of
       record linkage.

    """
    jaro_weights = check_jaro_weights(jaro_weights)
    distance = C.jaro_winkler_sd(a, b, p=p, max_l=max_l, jaro_weights=jaro_weights)
    return distance


@single_dispatch
def longest_common_substring(a: str, b: str) -> float:
    """
    Compute the LCS distance [1]_ between two sequences.

    {params}

    {returns}

    Examples
    --------
    >>> longest_common_substring('AASQ', 'PASQ')

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Longest_common_substring_problem

    """
    distance = C.longest_common_substring_sd(a, b)
    return distance


@single_dispatch
def optimal_string_alignment(a: str, b: str) -> float:
    """
    Compute the OSA [1]_ between two sequences.

    {params}

    {returns}

    Examples
    --------
    >>> optimal_string_alignment('AASQ', 'PASQ')

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance

    """
    distance = C.optimal_string_alignment_sd(a, b)
    return distance
