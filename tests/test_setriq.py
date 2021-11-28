import decimal

import pytest
import srsly

import setriq

decimal.getcontext().prec = 4
ROUNDING = decimal.Decimal('0.0001')

test_cases = [
    (['AASQ', 'PASQ'], [decimal.Decimal('0.3564')]),
    (['GTA', 'HLA', 'KKR'], [decimal.Decimal(''),
                             decimal.Decimal(''),
                             decimal.Decimal('')]),
    (['SEQVENCES', 'SEQVENCES'], [decimal.Decimal('0.0')])
]


@pytest.fixture()
def blosum_62():
    def _method():
        out = srsly.read_json('data/blosum-62.json')
        return out
    
    return _method


@pytest.mark.parametrize('test_case', test_cases)
def test_cdr_dist(test_case, blosum_62):
    blosum = blosum_62()

    sequences, distances = test_case
    response = setriq.cdr_dist(sequences, **blosum)

    n = len(sequences)
    assert len(response) == (n * (n - 1) / 2)
    assert all(decimal.Decimal(r).quantize(ROUNDING, rounding=decimal.ROUND_HALF_UP) == tgt
               for r, tgt in zip(response, distances))
