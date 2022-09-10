import pytest

from setriq import BLOSUM62
from setriq import single_dispatch


class Cases:
    SEQUENCES = [
        ("CASSLKPNTEAFF", "CASSAHIANYGYTF"),
        ("CASSLKPNTEAFF", "CASRGATETQYF"),
        ("CASSAHIANYGYTF", "CASRGATETQYF"),
    ]
    EQUAL_SEQUENCE_LENGTH = [
        ("AASQ", "PASQ"),
        ("GTA", "HLA"),
        ("SEQVENCES", "SEQVENCES"),
    ]


class Results:
    CDR_DIST = [0.7121679380884349, 0.6498905737037513, 0.75209911355881]
    LEVENSHTEIN = [8.0, 8.0, 9.0]
    TCR_DIST = [4.0, 8.0, 0.0]
    HAMMING = [1.0, 2.0, 0.0]
    JARO = [0.3335622710622711, 0.3641636141636142, 0.3373015873015873]
    JARO_WINKLER = [0.20013736263736265, 0.25491452991452995, 0.2361111111111111]
    LONGEST_COMMON_SUBSTRING = [
        13.0,
        13.0,
        14.0,
    ]
    OPTIMAL_STRING_ALIGNMENT = [
        8.0,
        8.0,
        9.0,
    ]


@pytest.mark.parametrize(
    ["sequences", "distance"], zip(Cases.SEQUENCES, Results.CDR_DIST)
)
def test_cdr_dist(sequences, distance):
    result = single_dispatch.cdr_dist(*sequences)
    assert result == distance


@pytest.mark.parametrize(
    ["sequences", "distance"], zip(Cases.SEQUENCES, Results.LEVENSHTEIN)
)
def test_levenshtein(sequences, distance):
    result = single_dispatch.levenshtein(*sequences)
    assert result == distance


@pytest.mark.parametrize(
    ["sequences", "distance"], zip(Cases.EQUAL_SEQUENCE_LENGTH, Results.TCR_DIST)
)
def test_tcr_dist(sequences, distance):
    a, b = sequences
    component = {"substitution_matrix": BLOSUM62, "gap_penalty": 4.0}
    result = single_dispatch.tcr_dist({"cmp_1": a}, {"cmp_1": b}, cmp_1=component)
    assert result == distance


@pytest.mark.parametrize(
    ["sequences", "distance"], zip(Cases.EQUAL_SEQUENCE_LENGTH, Results.HAMMING)
)
def test_hamming(sequences, distance):
    result = single_dispatch.hamming(*sequences)
    assert result == distance


@pytest.mark.parametrize(["sequences", "distance"], zip(Cases.SEQUENCES, Results.JARO))
def test_jaro(sequences, distance):
    result = single_dispatch.jaro(*sequences)
    assert result == distance


@pytest.mark.parametrize(
    ["sequences", "distance"], zip(Cases.SEQUENCES, Results.JARO_WINKLER)
)
def test_jaro_winkler(sequences, distance):
    result = single_dispatch.jaro_winkler(*sequences, p=0.10)
    assert result == distance


@pytest.mark.parametrize(
    ["sequences", "distance"], zip(Cases.SEQUENCES, Results.LONGEST_COMMON_SUBSTRING)
)
def test_longest_common_substring(sequences, distance):
    result = single_dispatch.longest_common_substring(*sequences)
    assert result == distance


@pytest.mark.parametrize(
    ["sequences", "distance"], zip(Cases.SEQUENCES, Results.OPTIMAL_STRING_ALIGNMENT)
)
def test_optimal_string_alignment(sequences, distance):
    result = single_dispatch.optimal_string_alignment(*sequences)
    assert result == distance
