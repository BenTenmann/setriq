"""
Python API for sequence distance functions.

"""

import abc
import warnings
from typing import (
    Any,
    Dict,
    Generic,
    List,
    Optional,
    Sequence,
    Set,
    TypeVar,
    Union,
)

import numpy as np
import numpy.typing as npt
import pandas as pd
from glom import glom
from scipy import spatial
from sklearn import preprocessing

import setriq._C as C

from .substitution import BLOSUM45, SubstitutionMatrix
from .utils import (
    TCR_DIST_DEFAULT,
    TcrDistDef,
    check_jaro_weights,
    check_jaro_winkler_params,
    enforce_list,
    ensure_equal_sequence_length,
)

__all__ = [
    "CdrDist",
    "Levenshtein",
    "TcrDist",
    "TcrDistComponent",
    "Hamming",
    "Jaro",
    "JaroWinkler",
    "LongestCommonSubstring",
    "OptimalStringAlignment",
]

FloatArray = npt.NDArray[np.float64]
SeqRecord = TypeVar("SeqRecord", bound=Union[str, Dict[str, str]])


class Metric(abc.ABC, Generic[SeqRecord]):
    """
    The Metric abstract base class. Users familiar with the torch paradigm will recognize the overall structure of
    the metric subclasses.

    Methods
    -------
    forward(self, *args, **kwargs):
        an abstract method which needs to be implemented in every subclass. It is accessed via the ``__call__`` method
        of the base class.

    """

    call_args: Dict[str, Any]

    def __init__(self, return_squareform: bool = False):
        self.return_squareform = return_squareform

    @abc.abstractmethod
    def forward(self, sequences: Sequence[SeqRecord]) -> List[float]:
        pass

    @enforce_list(argnum=1, convert_iterable=True)
    def __call__(self, sequences: Sequence[SeqRecord]) -> FloatArray:
        out = np.array(self.forward(sequences))
        if self.return_squareform:
            out = spatial.distance.squareform(out)

        return out

    def to_sklearn(self) -> preprocessing.FunctionTransformer:
        """Creates a FunctionTransformer from a given Metric instance.

        This is useful for integrating the distance functions into scikit-learn pipelines.

        Returns
        -------
        The FunctionTransformer which can be used with other scikit-learn tooling.

        """
        return preprocessing.FunctionTransformer(
            func=self,
            check_inverse=False,
        )


class CdrDist(Metric[str]):
    """
    The CdrDist [1]_ class. Inherits from Metric.

    Examples
    --------
    >>> sequences = ['CASSLKPNTEAFF', 'CASSAHIANYGYTF', 'CASRGATETQYF']
    >>>
    >>> metric = CdrDist()
    >>> distances = metric(sequences)

    References
    ----------
    .. [1] Thakkar, N. and Bailey-Kellogg, C., 2019. Balancing sensitivity and specificity in distinguishing TCR groups
       by CDR sequence similarity. BMC bioinformatics, 20(1), pp.1-14. (https://doi.org/10.1186/s12859-019-2864-8)

    """

    def __init__(
        self,
        substitution_matrix: SubstitutionMatrix = BLOSUM45,
        gap_opening_penalty: float = 10.0,
        gap_extension_penalty: float = 1.0,
        return_squareform: bool = False,
    ):
        super(CdrDist, self).__init__(return_squareform)
        self.call_args = {
            **substitution_matrix,
            "gap_opening_penalty": gap_opening_penalty,
            "gap_extension_penalty": gap_extension_penalty,
        }
        self.fn = C.cdr_dist

    def forward(self, sequences: Sequence[str]) -> List[float]:
        out = self.fn(sequences, **self.call_args)

        return out


class Levenshtein(Metric[str]):
    """
    The Levenshtein [1]_ class. Inherits from Metric. It uses a refactor of the ``python-Levenshtein`` [2]_
    implementation in the backend.

    Examples
    --------
    >>> sequences = ['CASSLKPNTEAFF', 'CASSAHIANYGYTF', 'CASRGATETQYF']
    >>>
    >>> metric = Levenshtein()
    >>> distances = metric(sequences)

    References
    ----------
    .. [1] Levenshtein, V.I., 1966, February. Binary codes capable of correcting deletions, insertions, and reversals.
       In Soviet physics doklady (Vol. 10, No. 8, pp. 707-710).
    .. [2] python-Levenshtein (https://github.com/ztane/python-Levenshtein)

    """

    def __init__(self, extra_cost: float = 0.0, return_squareform: bool = False):
        super(Levenshtein, self).__init__(return_squareform)
        self.call_args = {"extra_cost": extra_cost}
        self.fn = C.levenshtein

    def forward(self, sequences: Sequence[str]) -> List[float]:
        out = self.fn(sequences, **self.call_args)

        return out


class TcrDistComponent(Metric[str]):
    """
    The TcrDistComponent class. Inherits from Metric.

    Examples
    --------
    >>> from setriq import BLOSUM62
    >>> sequences = ['CASSLKPNTE', 'CASS-HIANY', 'CASRGAT--Q']
    >>>
    >>> metric = TcrDistComponent(substitution_matrix=BLOSUM62, gap_penalty=4., gap_symbol='-', weight=1.)
    >>> distances = metric(sequences)

    """

    def __init__(
        self,
        substitution_matrix: SubstitutionMatrix,
        gap_penalty: float,
        gap_symbol: str = "-",
        weight: float = 1.0,
        return_squareform: bool = False,
    ):
        """
        Initialize a TcrDistComponent object.

        Parameters
        ----------
        substitution_matrix : SubstitutionMatrix
            a SubstitutionMatrix object
        gap_penalty : float
            the gap penalty, i.e. the score given for a <not-gap-symbol> -> <gap-symbol> substitution
        gap_symbol : str
            the gap symbol (default = '-')
        weight : float
            the weighting of the component weight

        """
        super(TcrDistComponent, self).__init__(return_squareform)
        self.call_args = {
            **substitution_matrix,
            "gap_penalty": gap_penalty,
            "gap_symbol": gap_symbol,
            "weight": weight,
        }
        self.fn = C.tcr_dist_component

    @ensure_equal_sequence_length(argnum=1)
    def forward(self, sequences: Sequence[str]) -> List[float]:
        out = self.fn(sequences, **self.call_args)

        return out


class TcrDist(Metric[Dict[str, str]]):
    """
    TcrDist [1]_ class. Inherits from Metric. It is a container class for individual TcrDistComponent instances.
    Components are executed sequentially and their results aggregated at the end (summation).

    Attributes
    ----------
    components : List[str]
        holds the names of the components to be executed

    Examples
    --------
    >>> sequences = [
    ...     {'cdr_1': 'TSG------FNG', 'cdr_2': 'VVL----DGL', 'cdr_2_5': 'SRSN-GY', 'cdr_3': 'CAVR-----'},
    ...     {'cdr_1': 'TSG------FYG', 'cdr_2': 'NGL----DGL', 'cdr_2_5': 'SRSD-SY', 'cdr_3': 'CA-------'},
    ...     {'cdr_1': 'NSA------FQY', 'cdr_2': 'TYS----SGN', 'cdr_2_5': 'DKSSKY-', 'cdr_3': 'CAMS-----'}
    ... ]
    >>> metric = TcrDist()  # will produce a warning stating default configuration (Dash et al)
    >>> distances = metric(sequences)

    References
    ----------
    .. [1] Dash, P., Fiore-Gartland, A.J., Hertz, T., Wang, G.C., Sharma, S., Souquette, A., Crawford, J.C., Clemens,
       E.B., Nguyen, T.H., Kedzierska, K. and La Gruta, N.L., 2017. Quantifiable predictive features define
       epitope-specific T cell receptor repertoires. Nature, 547(7661), pp.89-93. (https://doi.org/10.1038/nature22383)

    """

    _default: TcrDistDef = TCR_DIST_DEFAULT
    _default_msg: str = (
        "TcrDist has been initialized using the default configuration. "
        "Please ensure that the input is a list of dictionaries, each with keys: {}"
    ).format(", ".join(repr(key) for key, _ in _default))

    def __init__(self, return_squareform: bool = False, **components: TcrDistComponent):
        """
        Initialize a TcrDist object. Initialization can happen in two ways:
            1. no arguments are passed, instantiating the default configuration of TcrDist -- i.e. the configuration
               described in Dash et al.
            2. a set of keyword arguments is passed, where each value is a TcrDistComponent instance and the key gives
               the component its name. An arbitrary number of components can be set and their naming is also arbitrary.

        Parameters
        ----------
        components : keyword arguments
            either a set of keyword arguments, where each value is a TcrDistComponent instance which will be stored as
            an attribute with the key as its name OR `None` -- in which case the default configuration is loaded (Dash
            et al).

        Examples
        --------

        Initialize with default parameters

        >>> metric = TcrDist()  # will produce a warning

        Initialize with custom components

        >>> component_1 = TcrDistComponent(substitution_matrix=BLOSUM45, gap_penalty=3, weight=1)
        >>> component_2 = TcrDistComponent(substitution_matrix=BLOSUM45, gap_penalty=5, weight=10)
        >>>
        >>> metric = TcrDist(cmp_1=component_1, cmp_2=component_2)

        Keep in mind that the keys will be used to assiociate the components to the relevant input fields, i.e. in this
        case the input should take the shape of:

        >>> [{'cmp_1': '...', 'cmp_2': '...'}, ...]

        additional keys will have no effect.

        """
        super(TcrDist, self).__init__(return_squareform)
        parts: List[str] = []

        # user-defined configuration
        if components:
            for name, component in components.items():
                # some type checking
                if not isinstance(component, TcrDistComponent):
                    raise TypeError(
                        f"{name!r} is not of type {TcrDistComponent.__class__.__name__}"
                    )

                setattr(self, name, component)
                parts.append(name)

        # default configuration
        else:
            for name, definition in self._default:
                setattr(self, name, TcrDistComponent(**definition))  # type: ignore[arg-type]
                parts.append(name)

            # warn user that default has been initialised and inform required input format
            warnings.warn(self._default_msg, UserWarning)

        self.components = parts

    def _check_input_format(self, ipt: pd.Index) -> None:
        pts: Set[str] = set(self.components)

        diff: Set[str] = pts.difference(ipt)
        if diff:
            raise ValueError("Missing key(s): {}".format(", ".join(map(repr, diff))))

    @property
    def required_input_keys(self) -> List[str]:
        """
        Get the keys (=fields) required in the input to TcrDist instance.

        Returns
        -------
        required_input_keys : List[str]
            returns a list of strings signifying the keys required in the input

        """
        return self.components

    @property
    def default_definition(self) -> TcrDistDef:
        """
        Get the default TcrDistComponent schema as defined by Dash et al.

        Returns
        -------
        default_schema : List[tuple]
            returns the schema for the TcrDistComponent instances held in the default configuration

        """
        return self._default

    def forward(self, sequences: Sequence[Dict[str, str]]) -> List[float]:
        if not sequences:
            return []

        # check the input keys provided -- assumes consistency
        self._check_input_format(pd.DataFrame(sequences).columns)

        # iterate through components and collect component output
        out: List[FloatArray] = []
        for part in self.components:
            # gather sequences of the associated field into a list
            sqs: List[str] = glom(sequences, [part])
            component: TcrDistComponent = getattr(self, part)

            # execute component on list of associated sequences
            result: FloatArray = component(sqs)
            out.append(result)

        # aggregate the component outputs
        res: FloatArray = np.stack(out).sum(axis=0)
        return res.tolist()


class Hamming(Metric[str]):
    """
    Hamming distance class. Inherits from Metric. Sequences must be of equal length.

    Examples
    --------
    >>> metric = Hamming(mismatch_score=2.0)
    >>> sequences = ['CASSLKPNTEAFF', 'CASSAHIANYGYTF', 'CASRGATETQYF']
    >>> distances = metric(sequences)

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Hamming_distance
    """

    def __init__(self, mismatch_score: float = 1.0, return_squareform: bool = False):
        super(Hamming, self).__init__(return_squareform)
        self.call_args = {"mismatch_score": mismatch_score}
        self.fn = C.hamming

    @ensure_equal_sequence_length(argnum=1)
    def forward(self, sequences: Sequence[str]) -> List[float]:
        out = self.fn(sequences, **self.call_args)
        return out


class Jaro(Metric[str]):
    """
    Jaro distance class. Inherits from Metric. Adapted from [2].

    Examples
    --------
    >>> metric = Jaro()
    >>> sequences = ['CASSLKPNTEAFF', 'CASSAHIANYGYTF', 'CASRGATETQYF']
    >>> distances = metric(sequences)

    References
    ----------
    [1] Jaro, M.A., 1989. Advances in record-linkage methodology as applied to matching the 1985 census of Tampa,
        Florida. Journal of the American Statistical Association, 84(406), pp.414-420.
    [2] Van der Loo, M.P., 2014. The stringdist package for approximate string matching. R J., 6(1), p.111.
    """

    def __init__(
        self,
        jaro_weights: Optional[List[float]] = None,
        return_squareform: bool = False,
    ):
        super(Jaro, self).__init__(return_squareform)
        jaro_weights = check_jaro_weights(jaro_weights)
        self.call_args = {"jaro_weights": jaro_weights}
        self.fn = C.jaro

    def forward(self, sequences: Sequence[str]) -> List[float]:
        out = self.fn(sequences, **self.call_args)
        return out


class JaroWinkler(Jaro):
    """
    Jaro-Winkler distance class. Inherits from Jaro.

    Examples
    --------
    >>> metric = JaroWinkler(p=0.10)
    >>> sequences = ['CASSLKPNTEAFF', 'CASSAHIANYGYTF', 'CASRGATETQYF']
    >>> distances = metric(sequences)

    References
    ----------
    .. [1] Winkler, W.E., 1990. String comparator metrics and enhanced decision rules in the Fellegi-Sunter model of
       record linkage.

    """

    @check_jaro_winkler_params
    def __init__(
        self,
        p: float = 0.10,
        max_l: int = 4,
        jaro_weights: Optional[List[float]] = None,
        return_squareform: bool = False,
    ):
        super(JaroWinkler, self).__init__(
            jaro_weights=jaro_weights, return_squareform=return_squareform
        )
        self.call_args["p"] = p
        self.call_args["max_l"] = max_l
        self.fn = C.jaro_winkler  # type: ignore[assignment]


class LongestCommonSubstring(Metric[str]):
    """
    Longest common substring distance class. Inherits from Metric.

    Examples
    --------
    >>> metric = LongestCommonSubstring()
    >>> sequences = ['CASSLKPNTEAFF', 'CASSAHIANYGYTF', 'CASRGATETQYF']
    >>> distances = metric(sequences)

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Longest_common_substring_problem

    """

    def __init__(self, return_squareform: bool = False):
        super(LongestCommonSubstring, self).__init__(return_squareform)
        self.fn = C.longest_common_substring

    def forward(self, sequences: Sequence[str]) -> List[float]:
        out = self.fn(sequences)
        return out


class OptimalStringAlignment(Metric[str]):
    """
    Optimal string alignment distance class. Inherits from Metric.

    Examples
    --------
    >>> metric = OptimalStringAlignment()
    >>> sequences = ['CASSLKPNTEAFF', 'CASSAHIANYGYTF', 'CASRGATETQYF']
    >>> distances = metric(sequences)

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance

    """

    def __init__(self, return_squareform: bool = False):
        super(OptimalStringAlignment, self).__init__(return_squareform)
        self.fn = C.optimal_string_alignment

    def forward(self, sequences: Sequence[str]) -> List[float]:
        out = self.fn(sequences)
        return out
