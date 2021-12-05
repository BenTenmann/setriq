import decimal
import warnings

import pytest

import setriq

ROUNDING = decimal.Decimal('0.0001')

test_cases = [
    ['AASQ', 'PASQ'],
    ['GTA', 'HLA', 'KKR'],
    ['SEQVENCES', 'SEQVENCES']
]

cdr_dist_results = [
    [decimal.Decimal('0.2950')],
    [decimal.Decimal('0.7418'), decimal.Decimal('1.0'), decimal.Decimal('1.0')],
    [decimal.Decimal('0.0')],
]

tcr_dist_results = [
    [decimal.Decimal('24.0')],
    [decimal.Decimal('48.0'), decimal.Decimal('72.0'), decimal.Decimal('72.0')],
    [decimal.Decimal('0.0')],
]


@pytest.fixture()
def cdr_dist():
    def _method():
        return setriq.CdrDist()

    return _method


@pytest.fixture()
def tcr_dist_base():
    def _method():
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            metric = setriq.TcrDist()

        return metric

    return _method


def convert_to_tcr_dist_format(tc, rs):
    out = zip(([{'cdr_1': seq, 'cdr_2': seq, 'cdr_2_5': seq, 'cdr_3': seq} for seq in case]
               for case in tc), rs)

    return out


@pytest.mark.parametrize('test_case', zip(test_cases, cdr_dist_results))
def test_cdr_dist(cdr_dist, test_case):
    metric = cdr_dist()

    sequences, distances = test_case
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = [decimal.Decimal(r).quantize(ROUNDING, rounding=decimal.ROUND_HALF_UP) for r in response]
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize('test_case', convert_to_tcr_dist_format(test_cases, tcr_dist_results))
def test_tcr_dist(tcr_dist_base, test_case):
    metric = tcr_dist_base()

    sequences, distances = test_case
    response = metric(sequences)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)

    res = [decimal.Decimal(r).quantize(ROUNDING, rounding=decimal.ROUND_HALF_UP) for r in response]
    assert all(r == tgt for r, tgt in zip(res, distances))


@pytest.mark.parametrize('test_case', zip(test_cases, tcr_dist_results))
def test_tcr_dist_error(tcr_dist_base, test_case):
    metric = tcr_dist_base()

    sequences, distances = test_case
    with pytest.raises(ValueError):
        metric(sequences)
