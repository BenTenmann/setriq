from itertools import product

import pytest

import setriq

# ------ Test Examples ----------------------------------------------------------------------------------------------- #
hello_world = [
    (("hello", "hello"), 1.0),
    (("hello", "world"), 0.0),
    (("world", "hello"), 0.0),
    (("world", "world"), 1.0),
]

file = {
    "complete_0": {
        "index": {"hello": 0, "world": 1},
        "substitution_matrix": [[1.0, 0.0], [0.0, 1.0]],
    },
    "complete_1": {
        "index": {"hello": 0, "world": 1},
        "substitution_matrix": [[1.0, 0.0], [0.0, 1.0]],
        "some_var": 0xFFFF,
    },
    "incomplete_0": {"index": {"hello": 0, "world": 1}},
    "incomplete_1": {"substitution_matrix": [[1.0, 0.0], [0.0, 1.0]]},
}
complete_files = ["complete_0", "complete_1"]

incomplete_files = ["incomplete_0", "incomplete_1"]

add_tokens_single_value = [("-", 1.0), ("+", 1.0), ("*", 1.0)]

add_tokens_multiple_values = [
    ("-", [0.0, 0.0, 1.0]),
    ("+", [0.0, 0.0, 1.0]),
    ("*", [0.0, 0.0, 1.0]),
]

add_tokens_existing = [("hello", 1.0), ("world", 1.0)]

add_tokens_multiple_values_wrong_dim = [
    ("-", [0.0, 1.0]),
    ("+", [0.0, 0.0, 1.0, 0.0]),
    ("*", [0.0]),
]


# ------ Fixtures ---------------------------------------------------------------------------------------------------- #
@pytest.fixture()
def substitution_matrix_parts():
    def _method():
        idx = {"hello": 0, "world": 1}
        scoring = [[1.0, 0.0], [0.0, 1.0]]

        return idx, scoring

    return _method


@pytest.fixture()
def mock_json_reader(monkeypatch):
    monkeypatch.setattr(
        "setriq.modules.substitution.srsly.read_json",
        lambda name, *args, **kwargs: file[str(name)],
    )


# ------ Tests ------------------------------------------------------------------------------------------------------- #
@pytest.mark.parametrize(["sub", "tgt"], hello_world)
def test_substitution_matrix_init(substitution_matrix_parts, sub, tgt):
    idx, scoring = substitution_matrix_parts()

    sm = setriq.SubstitutionMatrix(index=idx, substitution_matrix=scoring)
    assert hasattr(sm, "index")
    assert hasattr(sm, "substitution_matrix")
    assert len(sm) == len(scoring)

    assert sm(*sub) == tgt


@pytest.mark.parametrize(
    ["file_name", "test_case"], product(complete_files, hello_world)
)
def test_substitution_matrix_from_json(mock_json_reader, file_name, test_case):
    sm = setriq.SubstitutionMatrix.from_json(file_name)

    sub, tgt = test_case
    assert sm(*sub) == tgt


@pytest.mark.parametrize("file_name", incomplete_files)
def test_substitution_matrix_from_json_error(mock_json_reader, file_name):
    with pytest.raises(ValueError):
        setriq.SubstitutionMatrix.from_json(file_name)


@pytest.mark.parametrize(["token", "value"], add_tokens_single_value)
def test_substitution_matrix_add_token_single_value(
    substitution_matrix_parts, token, value
):
    idx, scoring = substitution_matrix_parts()

    sm = setriq.SubstitutionMatrix(index=idx, substitution_matrix=scoring)
    n = len(sm)

    lm = sm.add_token(token, value)
    assert len(lm) == (n + 1)
    assert lm(token, token) == value

    sm.add_token(token, value, inplace=True)
    assert len(sm) == (n + 1)
    assert sm(token, token) == value

    with pytest.raises(ValueError):
        sm.add_token(token, value)


@pytest.mark.parametrize(["token", "values"], add_tokens_multiple_values)
def test_substitution_matrix_add_token_multiple_values(
    substitution_matrix_parts, token, values
):
    idx, scoring = substitution_matrix_parts()

    sm = setriq.SubstitutionMatrix(index=idx, substitution_matrix=scoring)
    n = len(sm)

    lm = sm.add_token(token, values)
    assert len(lm) == (n + 1)
    assert lm(token, token) == values[-1]

    sm.add_token(token, values, inplace=True)
    assert len(sm) == (n + 1)
    assert sm(token, token) == values[-1]


@pytest.mark.parametrize(["token", "value"], add_tokens_existing)
def test_substitution_matrix_add_token_existing(
    substitution_matrix_parts, token, value
):
    idx, scoring = substitution_matrix_parts()

    sm = setriq.SubstitutionMatrix(index=idx, substitution_matrix=scoring)
    with pytest.raises(ValueError):
        sm.add_token(token, value)


@pytest.mark.parametrize(["token", "values"], add_tokens_multiple_values_wrong_dim)
def test_substitution_matrix_add_wrong_dim(substitution_matrix_parts, token, values):
    idx, scoring = substitution_matrix_parts()

    sm = setriq.SubstitutionMatrix(index=idx, substitution_matrix=scoring)
    with pytest.raises(ValueError):
        sm.add_token(token, values)
