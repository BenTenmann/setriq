import decimal as dc
import itertools
import warnings

import numpy as np
import pytest
from sklearn import preprocessing

import setriq

ROUNDING = dc.Decimal("0.0001")

# ------ Test Examples ----------------------------------------------------------------------------------------------- #
test_cases = [["AASQ", "PASQ"], ["GTA", "HLA", "KKR"], ["SEQVENCES", "SEQVENCES"]]

cdr_dist_results = [
    [dc.Decimal("0.3153")],
    [dc.Decimal("0.7288"), dc.Decimal("1.0"), dc.Decimal("1.0")],
    [dc.Decimal("0.0")],
]

levensthein_test_results = [
    [dc.Decimal("1.0")],
    [dc.Decimal("2.0"), dc.Decimal("3.0"), dc.Decimal("3.0")],
    [dc.Decimal("0.0")],
]

tcr_dist_component_results = [
    [dc.Decimal("4.0")],
    [dc.Decimal("8.0"), dc.Decimal("12.0"), dc.Decimal("12.0")],
    [dc.Decimal("0.0")],
]

tcr_dist_results = [
    [dc.Decimal("24.0")],
    [dc.Decimal("48.0"), dc.Decimal("72.0"), dc.Decimal("72.0")],
    [dc.Decimal("0.0")],
]

hamming_results = [
    [dc.Decimal("1.0")],
    [dc.Decimal("2.0"), dc.Decimal("3.0"), dc.Decimal("3.0")],
    [dc.Decimal("0.0")],
]

jaro_results = [
    [dc.Decimal("0.1667")],
    [dc.Decimal("0.4444"), dc.Decimal("1.0"), dc.Decimal("1.0")],
    [dc.Decimal("0.0")],
]

jaro_winkler_results = jaro_results  # change this in the future

longest_common_substring_results = [
    [dc.Decimal("2.0")],
    [dc.Decimal("4.0"), dc.Decimal("6.0"), dc.Decimal("6.0")],
    [dc.Decimal("0.0")],
]

optimal_string_alignment_results = [
    [dc.Decimal("1.0")],
    [dc.Decimal("2.0"), dc.Decimal("3.0"), dc.Decimal("3.0")],
    [dc.Decimal("0.0")],
]


# ------ Fixtures ---------------------------------------------------------------------------------------------------- #
@pytest.fixture()
def cdr_dist():
    def _method():
        return setriq.CdrDist()

    return _method


@pytest.fixture()
def levenshtein():
    def _method():
        return setriq.Levenshtein()

    return _method


@pytest.fixture()
def tcr_dist_component():
    def _method():
        return setriq.modules.distances.TcrDistComponent(
            substitution_matrix=setriq.BLOSUM62,
            gap_penalty=4.0,
            gap_symbol="-",
            weight=1.0,
        )

    return _method


@pytest.fixture()
def tcr_dist_keys():
    return ["cdr_1", "cdr_2", "cdr_2_5", "cdr_3"]


@pytest.fixture()
def tcr_dist_base():
    def _method():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            metric = setriq.TcrDist()

        return metric

    return _method


@pytest.fixture()
def tcr_dist_custom():
    def _method():
        cdr_1 = setriq.modules.distances.TcrDistComponent(
            substitution_matrix=setriq.BLOSUM62, gap_penalty=4.0
        )
        cdr_2 = setriq.modules.distances.TcrDistComponent(
            substitution_matrix=setriq.BLOSUM62, gap_penalty=4.0
        )
        cdr_2_5 = setriq.modules.distances.TcrDistComponent(
            substitution_matrix=setriq.BLOSUM62, gap_penalty=4.0
        )
        cdr_3 = setriq.modules.distances.TcrDistComponent(
            substitution_matrix=setriq.BLOSUM62, gap_penalty=8.0, weight=3.0
        )
        return setriq.TcrDist(cdr_1=cdr_1, cdr_2=cdr_2, cdr_2_5=cdr_2_5, cdr_3=cdr_3)

    return _method


@pytest.fixture()
def mock_abc(monkeypatch):
    monkeypatch.setattr("setriq.modules.distances.Metric.__abstractmethods__", set())


# ------ Helper Functions -------------------------------------------------------------------------------------------- #
def response_to_decimal(response):
    res = [
        dc.Decimal(r).quantize(ROUNDING, rounding=dc.ROUND_HALF_UP) for r in response
    ]

    return res


def convert_to_tcr_dist_format(tc, rs):
    out = zip(
        (
            [{"cdr_1": seq, "cdr_2": seq, "cdr_2_5": seq, "cdr_3": seq} for seq in case]
            for case in tc
        ),
        rs,
    )

    return out


# ------ Tests ------------------------------------------------------------------------------------------------------- #
def test_metric(mock_abc):
    # test abstract class passes
    metric = setriq.modules.distances.Metric()
    metric.forward([])


@pytest.mark.parametrize(["sequences", "distances"], zip(test_cases, cdr_dist_results))
def test_cdr_dist(cdr_dist, sequences, distances):
    metric = cdr_dist()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(
    ["sequences", "distances"], zip(test_cases, levensthein_test_results)
)
def test_levenshtein(levenshtein, sequences, distances):
    metric = levenshtein()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(
    ["sequences", "distances"], zip(test_cases, tcr_dist_component_results)
)
def test_tcr_dist_component(tcr_dist_component, sequences, distances):
    metric = tcr_dist_component()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(
    ["sequences", "distances"], convert_to_tcr_dist_format(test_cases, tcr_dist_results)
)
def test_tcr_dist(tcr_dist_base, tcr_dist_keys, sequences, distances):
    metric = tcr_dist_base()
    assert metric.required_input_keys == tcr_dist_keys
    assert all(
        k == key for (k, _), key in zip(metric.default_definition, tcr_dist_keys)
    )

    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(["sequences", "distances"], zip(test_cases, tcr_dist_results))
def test_tcr_dist_error(tcr_dist_base, sequences, distances):
    metric = tcr_dist_base()

    with pytest.raises(ValueError):
        metric(sequences)


@pytest.mark.parametrize(
    ["sequences", "distances"], convert_to_tcr_dist_format(test_cases, tcr_dist_results)
)
def test_tcr_dist_custom(tcr_dist_custom, tcr_dist_keys, sequences, distances):
    metric = tcr_dist_custom()
    assert metric.required_input_keys == tcr_dist_keys

    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


def test_tcr_dist_custom_error():
    rainbows = setriq.modules.distances.TcrDistComponent(setriq.BLOSUM62, 4.0)
    butterflies = {"wings": "beat"}
    with pytest.raises(TypeError):
        setriq.TcrDist(rainbows=rainbows, butterflies=butterflies)


@pytest.mark.parametrize(["sequences", "distances"], zip(test_cases, hamming_results))
def test_hamming(sequences, distances):
    metric = setriq.Hamming()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(["sequences", "distances"], zip(test_cases, jaro_results))
def test_jaro(sequences, distances):
    metric = setriq.Jaro()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(
    ["sequences", "distances"], zip(test_cases, jaro_winkler_results)
)
def test_jaro_winkler(sequences, distances):
    metric = setriq.JaroWinkler(p=0.10)
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(
    ["sequences", "distances"], zip(test_cases, longest_common_substring_results)
)
def test_longest_common_substring(sequences, distances):
    metric = setriq.LongestCommonSubstring()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(
    ["sequences", "distances"], zip(test_cases, optimal_string_alignment_results)
)
def test_optimal_string_alignment(sequences, distances):
    metric = setriq.OptimalStringAlignment()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(
    ["metric", "case"],
    itertools.product(
        [mt(return_squareform=True) for mt in (setriq.CdrDist,)], test_cases
    ),
)
def test_squareform(metric, case):
    res = metric(case)
    assert isinstance(res, np.ndarray)
    assert np.allclose(res, res.T)
    assert np.all(res >= 0.0)
    assert np.all(np.diag(res) == 0.0)


@pytest.mark.parametrize(
    ["metric", "case"],
    itertools.product([mt() for mt in (setriq.CdrDist,)], test_cases),
)
def test_to_sklearn(metric, case):
    estimator = metric.to_sklearn()
    assert isinstance(estimator, preprocessing.FunctionTransformer)
    assert isinstance(estimator.transform(case), np.ndarray)
