from itertools import product

import pytest
import setriq

# ------ Test Examples ----------------------------------------------------------------------------------------------- #
hello_world = [
    (('hello', 'hello'), 1.),
    (('hello', 'world'), 0.),
    (('world', 'hello'), 0.),
    (('world', 'world'), 1.)
]

file = {
    'complete_0': {'index': {'hello': 0, 'world': 1},
                   'substitution_matrix': [[1., 0.],
                                           [0., 1.]]},
    'complete_1': {'index': {'hello': 0, 'world': 1},
                   'substitution_matrix': [[1., 0.],
                                           [0., 1.]],
                   'some_var': 0xffff},
    'incomplete_0': {'index': {'hello': 0, 'world': 1}},
    'incomplete_1': {'substitution_matrix': [[1., 0.],
                                             [0., 1.]]}
}
complete_files = [
    'complete_0',
    'complete_1'
]

incomplete_files = [
    'incomplete_0',
    'incomplete_1'
]


# ------ Fixtures ---------------------------------------------------------------------------------------------------- #
@pytest.fixture()
def substitution_matrix_parts():
    def _method():
        idx = {'hello': 0, 'world': 1}
        scoring = [[1., 0.],
                   [0., 1.]]

        return idx, scoring

    return _method


@pytest.fixture()
def mock_json_reader(monkeypatch):
    monkeypatch.setattr(
        'setriq.modules.substitution.srsly.read_json',
        lambda name, *args, **kwargs: file[str(name)]
    )


# ------ Tests ------------------------------------------------------------------------------------------------------- #
@pytest.mark.parametrize(['sub', 'tgt'], hello_world)
def test_substitution_matrix_init(substitution_matrix_parts, sub, tgt):
    idx, scoring = substitution_matrix_parts()

    sm = setriq.SubstitutionMatrix(index=idx, substitution_matrix=scoring)
    assert hasattr(sm, 'index')
    assert hasattr(sm, 'substitution_matrix')
    assert len(sm) == len(scoring)

    assert sm(*sub) == tgt


@pytest.mark.parametrize(['file_name', 'test_case'], product(complete_files, hello_world))
def test_substitution_matrix_from_json(mock_json_reader, file_name, test_case):
    sm = setriq.SubstitutionMatrix.from_json(file_name)

    sub, tgt = test_case
    assert sm(*sub) == tgt


@pytest.mark.parametrize('file_name', incomplete_files)
def test_substitution_matrix_from_json_error(mock_json_reader, file_name):
    with pytest.raises(ValueError):
        setriq.SubstitutionMatrix.from_json(file_name)


def test_substitution_matrix_add_token():
    pass
