"""
Package utilities. Not meant for outside use.
"""

import enum
import inspect
from functools import WRAPPER_ASSIGNMENTS, wraps
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union

from .substitution import BLOSUM62, SubstitutionMatrix

__all__ = [
    "enforce_list",
    "ensure_equal_sequence_length",
    "single_dispatch",
    "tcr_dist_sd_component_check",
    "ensure_equal_sequence_length_sd",
    "check_jaro_weights",
    "check_jaro_winkler_params",
    "TCR_DIST_DEFAULT",
    "TcrDistDef",
]

TcrDistComponentDef = Dict[str, Union[float, SubstitutionMatrix]]
NamedTCRDD = Tuple[str, TcrDistComponentDef]
TcrDistDef = List[NamedTCRDD]

WRAPPER_ASSIGNMENTS = (*WRAPPER_ASSIGNMENTS, "__signature__")  # type: ignore[assignment]
TCR_DIST_DEFAULT: TcrDistDef = [
    ("cdr_1", {"substitution_matrix": BLOSUM62, "gap_penalty": 4.0, "weight": 1.0}),
    ("cdr_2", {"substitution_matrix": BLOSUM62, "gap_penalty": 4.0, "weight": 1.0}),
    ("cdr_2_5", {"substitution_matrix": BLOSUM62, "gap_penalty": 4.0, "weight": 1.0}),
    ("cdr_3", {"substitution_matrix": BLOSUM62, "gap_penalty": 8.0, "weight": 3.0}),
]


class Argument(enum.Enum):
    POSITIONAL = 0
    KEYWORD = 1
    DEFAULT = 2


def _get_argument_index_from_argnum(n_params: int, argnum: int):
    # updates the argument index provided -- fixes issues with `_get_argument` when `argnum` < 0 (negative indexing)
    if argnum < 0:
        return argnum + n_params
    return argnum


def _get_argument(params, argname: str, argnum: int, args: tuple, kwargs: dict):
    # convenience function for getting a specific argument from all the args and kwargs passed to an arbitrary function
    # in Python, any argument can be passed as either positional or keyword (except if explicitly encoded otherwise)
    # and so we need to account for this, as well as the default case, when the argument is not passed
    if argname in kwargs:
        # we also return where we found the argument, i.e. the enum describes the argument type. This helps us
        # downstream, if we want to put an augmented argument back into args, kwargs
        return kwargs[argname], Argument.KEYWORD
    if len(args) > argnum:
        return args[argnum], Argument.POSITIONAL
    return params[argname].default, Argument.DEFAULT


def _put_argument(
    params,
    argname: str,
    argnum: int,
    arg_type: enum.Enum,
    argval,
    args: tuple,
    kwargs: dict,
):
    # this function is the complement to `_get_argument`
    if arg_type == Argument.POSITIONAL:
        args = tuple(arg if argnum != idx else argval for idx, arg in enumerate(args))
    elif arg_type == Argument.KEYWORD:
        kwargs[argname] = argval
    elif arg_type == Argument.DEFAULT:
        if params[argname].default != argval:
            kwargs[argname] = argval
    return args, kwargs


def _get_func_argument_info(fn: Callable, argnum: int):
    signature = inspect.signature(fn)
    params = signature.parameters
    argname = list(params)[argnum]

    return signature, params, argname


def _get_func_argument_info_from_name(fn: Callable, argname: str):
    signature = inspect.signature(fn)
    params = signature.parameters
    argnum = list(params).index(argname)

    return signature, params, argnum


def _add_func_signature(fn: Callable, signature: inspect.Signature) -> Callable:
    setattr(fn, "__signature__", getattr(fn, "__signature__", None) or signature)
    return fn


def enforce_list(argnum: int = 0, convert_iterable: bool = True):
    """
    Enforce that a specified argument is always passed as a list to a given function. This is a decorator factory.

    Parameters
    ----------
    argnum: int
        The positional (integer) index of the argument to be forced into list format. Works like regular positional
        indexing. (default=0)
    convert_iterable: bool
        Boolean defining whether to force convert any iterable (except ``str``) into a list.

    Returns
    -------
    decorator: Callable
        A new decorator function, which can be used to decorate an arbitrary function / method.

    Examples
    --------

    A basic example:

    >>> @enforce_list()
    ... def f(x):
    ...     return x * 3
    >>>
    >>> f([3])
    ... [3, 3, 3]
    >>> f(3)
    ... [3, 3, 3]

    Notice, that we ommit the ``argnum`` parameter here. This is because ``argnum`` is 0 by default, i.e. it looks at
    the first argument passed to ``f``. Note, that ``enforce_list`` needs to be called before decorating a function.

    Enforce list can also be composed arbitrarily, to enforce multiple arguments to be lists:

    >>> @enforce_list(argnum=0)
    ... @enforce_list(argnum=1)
    ... def f(x, y):
    ...     return x + y

    Finally, ``enforce_list`` can also force convert other iterables (excluding str) into lists:

    >>> @enforce_list(argnum=0, convert_iterable=True)
    ... def f(x):
    ...     return x * 3
    >>>
    >>> x = np.array([1.])
    >>> f(x)
    ... [1., 1., 1.]

    """
    if callable(argnum):
        raise TypeError(
            f"Make sure to call {repr(enforce_list.__name__)} before decorating."
        )

    def decorator(fn):
        signature, params, argname = _get_func_argument_info(fn, argnum)
        argidx = _get_argument_index_from_argnum(len(params), argnum)
        fn = _add_func_signature(fn, signature)

        @wraps(fn, assigned=WRAPPER_ASSIGNMENTS)
        def _fn(*args, **kwargs):
            argument, arg_type = _get_argument(params, argname, argidx, args, kwargs)
            if isinstance(argument, Iterable) and not isinstance(argument, str):
                if not isinstance(argument, list) and convert_iterable:
                    argument = list(argument)
            else:
                argument = [argname]
            args, kwargs = _put_argument(
                params, argname, argidx, arg_type, argument, args, kwargs
            )
            out = fn(*args, **kwargs)
            return out

        return _fn

    return decorator


def ensure_equal_sequence_length(argnum: int):
    """
    Ensure that all input sequences are of equal length.

    Parameters
    ----------
    argnum: int
        The positional (integer) index of the argument to be forced into list format. Works like regular positional
        indexing. (``default=0``)

    Returns
    -------
    decorator: Callable
        The decorator used to wrap the function where all sequences ought to be of equal length.

    Examples
    --------
    >>> @ensure_equal_sequence_length(argnum=0)
    ... def f(sequences):
    ...     return sequences
    >>>
    >>> a = ['AASQ', 'PWSQ']  # sequences of equal length
    >>> b = ['GAT', 'AAFFD']  # sequences with varying length
    >>>
    >>> f(a)  # no error
    >>> f(b)  # error!

    """
    import pandas as pd

    def decorator(fn):
        signature, params, argname = _get_func_argument_info(fn, argnum)
        argidx = _get_argument_index_from_argnum(len(params), argnum)
        fn = _add_func_signature(fn, signature)

        @wraps(fn, assigned=WRAPPER_ASSIGNMENTS)
        def _fn(*args, **kwargs):
            argument, _ = _get_argument(params, argname, argidx, args, kwargs)
            if (
                argument
                and not (len(argument[0]) == pd.Series(argument).str.len()).all()
            ):
                raise ValueError("Sequences must be of equal length")
            out = fn(*args, **kwargs)
            return out

        return _fn

    return decorator


def ensure_equal_sequence_length_sd(fn: Callable) -> Callable:
    signature = inspect.signature(fn)
    fn = _add_func_signature(fn, signature)

    add_doc = """Note
    ----
    `a` and `b` must be of equal length.
    """

    fn.__doc__ = (fn.__doc__ or "") + add_doc

    @wraps(fn, assigned=WRAPPER_ASSIGNMENTS)
    def _fn(a, b, *args, **kwargs):
        if len(a) != len(b):
            raise ValueError("Sequences must be of equal length")
        out = fn(a, b, *args, **kwargs)
        return out

    return _fn


def single_dispatch(fn: Callable) -> Callable:
    signature = inspect.signature(fn)
    fn = _add_func_signature(fn, signature)

    empty_doc = f"""
    Compute the `{fn.__name__}` metric between two sequences.

    {{params}}

    {{returns}}
    """

    param_doc = """
    Parameters
    ----------
    a: str
        A sequence to be compared.
    b: str
        A sequence to be compared."""

    return_doc = """Returns
    -------
    distance: float
        The computed distance between the two sequences."""

    fn.__doc__ = (fn.__doc__ or empty_doc).format(params=param_doc, returns=return_doc)

    @wraps(fn, assigned=WRAPPER_ASSIGNMENTS)
    def _fn(a, b, *args, **kwargs):
        if any(not isinstance(sequence, str) for sequence in (a, b)):
            raise TypeError("`a` and `b` must be of type str")
        out = fn(a, b, *args, **kwargs)
        return out

    return _fn


def tcr_dist_sd_component_check(fn):
    # wrapper specifically used for the single dispatch tcr_dist function to check inputs (components)
    import os

    signature = inspect.signature(fn)
    fn = _add_func_signature(fn, signature)

    def _check_component(name, component: dict):
        essential_keys = ["substitution_matrix", "gap_penalty"]
        optional_keys = ["gap_symbol", "weight"]
        missing_keys = set(essential_keys).difference(component)
        if missing_keys:
            msg = ", ".join(map(repr, missing_keys))
            raise ValueError(f"missing keys in component def {repr(name)}: {msg}")

        init_types: List[Union[type, Tuple[type, ...]]] = [
            SubstitutionMatrix,
            (float, int),
            str,
            (float, int),
        ]
        for key, _type in zip(essential_keys + optional_keys, init_types):
            elem = component.get(key)
            if elem is not None and not isinstance(elem, _type):
                given_type = type(elem)
                raise TypeError(
                    f"{repr(key)} needs to be of type {repr(_type)}, not {repr(given_type)}"
                )

    @wraps(fn, assigned=WRAPPER_ASSIGNMENTS)
    def _fn(a, b, **component_def):
        if not os.environ.get("SKIP_TCR_DIST_COMPONENT_CHECK"):
            for name, component in component_def.items():
                _check_component(name, component)
            if not component_def:
                component_def = dict(TCR_DIST_DEFAULT)
            if set(component_def).difference(set(a).union(b)):
                raise ValueError(
                    "key mismatch between payloads (`a` and `b`) and defined components."
                )
        out = fn(a, b, **component_def)
        return out

    return _fn


def check_jaro_weights(weights: Optional[List[float]]) -> List[float]:
    if weights is None:
        weights = [1 / 3] * 3
    if len(weights) != 3:
        raise ValueError("`jaro_weights` has to be of length 3")
    if sum(weights) != 1.0:
        raise ValueError("`jaro_weights` has to sum to 1.0")
    return weights


def check_jaro_winkler_params(fn: Callable):
    # checks that Jaro-Winkler parameters are sensibly defined
    argname_p = "p"
    argname_l = "max_l"
    signature, params, argnum_p = _get_func_argument_info_from_name(fn, argname_p)
    _, _, argnum_l = _get_func_argument_info_from_name(fn, argname_l)
    fn = _add_func_signature(fn, signature)

    @wraps(fn, assigned=WRAPPER_ASSIGNMENTS)
    def _fn(*args, **kwargs):
        p, _ = _get_argument(params, argname_p, argnum_p, args, kwargs)
        max_l, _ = _get_argument(params, argname_l, argnum_l, args, kwargs)

        if not (0.0 <= p <= 0.25):
            raise ValueError("`p` must be in range [0.0, 0.25]")
        if max_l < 0 or not isinstance(max_l, int):
            raise ValueError("`max_l` must be a non-negative integer")
        out = fn(*args, **kwargs)
        return out

    return _fn
