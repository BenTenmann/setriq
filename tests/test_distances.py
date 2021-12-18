import decimal as dc
import warnings

import pytest

import setriq

ROUNDING = dc.Decimal('0.0001')

# ------ Test Examples ----------------------------------------------------------------------------------------------- #
test_cases = [
    ['AASQ', 'PASQ'],
    ['GTA', 'HLA', 'KKR'],
    ['SEQVENCES', 'SEQVENCES']
]

cdr_dist_results = [
    [dc.Decimal('0.3153')],
    [dc.Decimal('0.7288'), dc.Decimal('1.0'), dc.Decimal('1.0')],
    [dc.Decimal('0.0')],
]

levensthein_test_results = [
    [dc.Decimal('1.0')],
    [dc.Decimal('2.0'), dc.Decimal('3.0'), dc.Decimal('3.0')],
    [dc.Decimal('0.0')]
]

tcr_dist_component_results = [
    [dc.Decimal('4.0')],
    [dc.Decimal('8.0'), dc.Decimal('12.0'), dc.Decimal('12.0')],
    [dc.Decimal('0.0')]
]

tcr_dist_results = [
    [dc.Decimal('24.0')],
    [dc.Decimal('48.0'), dc.Decimal('72.0'), dc.Decimal('72.0')],
    [dc.Decimal('0.0')],
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
            gap_penalty=4.,
            gap_symbol='-',
            weight=1.
        )

    return _method


@pytest.fixture()
def tcr_dist_keys():
    return ['cdr_1', 'cdr_2', 'cdr_2_5', 'cdr_3']


@pytest.fixture()
def tcr_dist_base():
    def _method():
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            metric = setriq.TcrDist()

        return metric

    return _method


@pytest.fixture()
def tcr_dist_custom():
    def _method():
        cdr_1 = setriq.modules.distances.TcrDistComponent(
            substitution_matrix=setriq.BLOSUM62,
            gap_penalty=4.
        )
        cdr_2 = setriq.modules.distances.TcrDistComponent(
            substitution_matrix=setriq.BLOSUM62,
            gap_penalty=4.
        )
        cdr_2_5 = setriq.modules.distances.TcrDistComponent(
            substitution_matrix=setriq.BLOSUM62,
            gap_penalty=4.
        )
        cdr_3 = setriq.modules.distances.TcrDistComponent(
            substitution_matrix=setriq.BLOSUM62,
            gap_penalty=8.,
            weight=3.
        )
        return setriq.TcrDist(cdr_1=cdr_1, cdr_2=cdr_2, cdr_2_5=cdr_2_5, cdr_3=cdr_3)

    return _method


# ------ Helper Functions -------------------------------------------------------------------------------------------- #
def response_to_decimal(response):
    res = [dc.Decimal(r).quantize(ROUNDING, rounding=dc.ROUND_HALF_UP) for r in response]

    return res


def convert_to_tcr_dist_format(tc, rs):
    out = zip(([{'cdr_1': seq, 'cdr_2': seq, 'cdr_2_5': seq, 'cdr_3': seq} for seq in case]
               for case in tc), rs)

    return out


# ------ Tests ------------------------------------------------------------------------------------------------------- #
@pytest.mark.parametrize(['sequences', 'distances'], zip(test_cases, cdr_dist_results))
def test_cdr_dist(cdr_dist, sequences, distances):
    metric = cdr_dist()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(['sequences', 'distances'], zip(test_cases, levensthein_test_results))
def test_levenshtein(levenshtein, sequences, distances):
    metric = levenshtein()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(['sequences', 'distances'], zip(test_cases, tcr_dist_component_results))
def test_tcr_dist_component(tcr_dist_component, sequences, distances):
    metric = tcr_dist_component()
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(['sequences', 'distances'], convert_to_tcr_dist_format(test_cases, tcr_dist_results))
def test_tcr_dist(tcr_dist_base, tcr_dist_keys, sequences, distances):
    metric = tcr_dist_base()
    assert metric.required_input_keys == tcr_dist_keys
    assert all(k == key for (k, _), key in zip(metric.default_definition, tcr_dist_keys))

    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize(['sequences', 'distances'], zip(test_cases, tcr_dist_results))
def test_tcr_dist_error(tcr_dist_base, sequences, distances):
    metric = tcr_dist_base()

    with pytest.raises(ValueError):
        metric(sequences)


@pytest.mark.parametrize(['sequences', 'distances'], convert_to_tcr_dist_format(test_cases, tcr_dist_results))
def test_tcr_dist_custom(tcr_dist_custom, tcr_dist_keys, sequences, distances):
    metric = tcr_dist_custom()
    assert metric.required_input_keys == tcr_dist_keys

    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = response_to_decimal(response)
    assert all(r == tgt for r, tgt in zip(res, distances))


def test_tcr_dist_custom_error():
    rainbows = setriq.modules.distances.TcrDistComponent(setriq.BLOSUM62, 4.)
    butterflies = {'wings': 'beat'}
    with pytest.raises(TypeError):
        setriq.TcrDist(rainbows=rainbows, butterflies=butterflies)
