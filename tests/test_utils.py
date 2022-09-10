import pytest

import numpy as np

from setriq.modules import utils


class Cases:
    ENFORCE_LIST = [1, [], [1, 2, 3], "hello world", np.array([1, 2, 3])]
    ENSURE_SEQ_LEN = [["AASQ", "PSQ"], ["GTA", "LA", "KKR"], ["SQVENCES", "SEQVENCES"]]
    SINGLE_DISPATCH = [(1, 2), ("hello", 3), ([], [])]
    TCR_DIST_COMPONENT = [
        ({}, ValueError),  # checks the default error
        ({"comp_1": {"gap_penalty": 1.0}}, ValueError),  # checks missing keys
        ({"comp_1": {"substitution_matrix": [], "gap_penalty": 1.0}}, TypeError),
    ]
    JARO_WEIGHTS = [([1 / 3] * 3,), (None,), ([0.5, 0.3, 0.2],)]
    JARO_WINKLER_PARAMS = [
        ({"p": 0.50, "max_l": 4}, r"`p` must be in range \[0\.0, 0\.25\]"),
        ({"p": 0.10, "max_l": -1}, r"`max_l` must be a non-negative integer"),
    ]


@pytest.mark.parametrize("test_case", Cases.ENFORCE_LIST)
def test_enforce_list(test_case):
    @utils.enforce_list()
    def f(x):
        return x

    assert isinstance(f(test_case), list)


@pytest.mark.parametrize("test_case", Cases.ENSURE_SEQ_LEN)
def test_ensure_equal_sequence_length(test_case):
    @utils.ensure_equal_sequence_length(argnum=0)
    def f(x):
        return x

    with pytest.raises(ValueError):
        f(test_case)


@pytest.mark.parametrize("test_case", Cases.SINGLE_DISPATCH)
def test_single_dispatch(test_case):
    @utils.single_dispatch
    def f(a, b):
        return a, b

    with pytest.raises(TypeError):
        f(*test_case)


@pytest.mark.parametrize(["arguments", "exception"], Cases.TCR_DIST_COMPONENT)
def test_tcr_dist_component_check_sd(arguments, exception):
    @utils.tcr_dist_sd_component_check
    def f(a, b, **component_df):
        return a, b, component_df

    with pytest.raises(exception):
        f({"comp_1": ""}, {"comp_1": ""}, **arguments)


@pytest.mark.parametrize("test_case", Cases.ENSURE_SEQ_LEN)
def test_ensure_equal_sequence_length_sd(test_case):
    @utils.ensure_equal_sequence_length_sd
    def f(a, b):
        return a, b

    with pytest.raises(ValueError):
        f(*test_case[:2])


@pytest.mark.parametrize("test_case", Cases.JARO_WEIGHTS)
def test_check_jaro_weights(test_case):
    @utils.check_jaro_weights
    def f(a, b, jaro_weights=None):
        return jaro_weights

    (arg,) = test_case
    assert f("", "", *test_case) == (arg or [1 / 3] * 3)
    assert f("", "", jaro_weights=arg) == (arg or [1 / 3] * 3)

    with pytest.raises(ValueError, match="`jaro_weights` has to be of length 3"):
        f("", "", (arg or [1 / 3] * 3)[:2])

    with pytest.raises(ValueError, match="`jaro_weights` has to sum to 1.0"):
        f("", "", [1] * 3 if not arg else [elem + 1 for elem in arg])


@pytest.mark.parametrize(["arguments", "exception_message"], Cases.JARO_WINKLER_PARAMS)
def test_check_jaro_winkler_params(arguments, exception_message):
    @utils.check_jaro_winkler_params
    def f(p, max_l):
        return p, max_l

    with pytest.raises(ValueError, match=exception_message):
        f(**arguments)
