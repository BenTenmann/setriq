"""
Python API for sequence distance functions.

"""

import abc
import warnings
from typing import Dict, List

import numpy as np
from glom import glom

import setriq._C as C
from ._substitution import *

__all__ = [
    'CdrDist',
    'TcrDist',
    'TcrDistComponent',
]


class Metric(abc.ABC):
    """
    The `Metric` abstract base class. Users familiar with the `torch` paradigm will recognize the overall structure of
    the metric subclasses.

    Methods
    -------
    forward(self, *args, **kwargs):
        an abstract method which needs to be implemented in every subclass. It is accessed via the `__call__` method of
        the base class

    """

    @abc.abstractmethod
    def forward(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        out = self.forward(*args, **kwargs)

        return out


class CdrDist(Metric):
    """
    The `CdrDist` class. Inherits from `Metric`.

    Examples
    --------
    >>> sequences = ['CASSLKPNTEAFF', 'CASSAHIANYGYTF', 'CASRGATETQYF']
    >>>
    >>> metric = CdrDist()
    >>> distances = metric(sequences)

    References
    ----------
    [1] Thakkar, N. and Bailey-Kellogg, C., 2019. Balancing sensitivity and specificity in distinguishing TCR groups by
        CDR sequence similarity. BMC bioinformatics, 20(1), pp.1-14. (https://doi.org/10.1186/s12859-019-2864-8)
    """
    def __init__(self, substitution_matrix: SubstitutionMatrix = BLOSUM62, gap_penalty: float = 1.):
        self.call_args = {
            **substitution_matrix,
            # implement gap-penalty option
        }
        self.fn = C.cdr_dist

    def forward(self, sequences: List[str]) -> List[float]:
        out = self.fn(sequences, **self.call_args)

        return out


class TcrDistComponent(Metric):
    """
    The TcrDistComponent class.

    Examples
    --------
    >>> sequences = ['CASSLKPNTE', 'CASS-HIANY', 'CASRGAT--Q']
    >>>
    >>> metric = TcrDistComponent(substitution_matrix=BLOSUM62, gap_penalty=4., gap_symbol='-', weight=1.)
    >>> distances = metric(sequences)
    """

    def __init__(self, 
                 substitution_matrix: SubstitutionMatrix, 
                 gap_penalty: float,
                 gap_symbol: str = '-',
                 weight: float = 1.):
        self.call_args = {
            **substitution_matrix,
            'gap_penalty': gap_penalty,
            'gap_symbol': gap_symbol,
            'weight': weight
        }
        self.fn = C.tcr_dist_component

    def forward(self, sequences: List[str]) -> List[float]:
        out = self.fn(sequences, **self.call_args)

        return out


class TcrDist(Metric):
    """
    TcrDist class. This is a container class for individual TcrDistComponent instances. Components are executed
    sequentially and their results aggregated at the end (summation).

    Attributes
    ----------
    components : List[str]
        holds the names of the components to be executed
    """
    _default = [
        ('cdr_1', {'substitution_matrix': BLOSUM62, 'gap_penalty': 4., 'weight': 1.}),
        ('cdr_2', {'substitution_matrix': BLOSUM62, 'gap_penalty': 4., 'weight': 1.}),
        ('cdr_3', {'substitution_matrix': BLOSUM62, 'gap_penalty': 8., 'weight': 3.})
    ]
    _default_msg = (
        'TcrDist has been initialized using the default configuration. '
        'Please ensure that the input is a list of dictionaries, each with keys: {}'
    ).format(', '.join(repr(key) for key, _ in _default))
    
    def __init__(self, **components):
        parts = []
        if components:
            for name, component in components.items():
                if not isinstance(component, TcrDistComponent):
                    raise TypeError(f'{repr(name)} is not of type {TcrDistComponent}')

                self.__setattr__(name, component)
                parts.append(name)
        else:
            for name, definition in self._default:
                self.__setattr__(name, TcrDistComponent(**definition))
                parts.append(name)

            warnings.warn(self._default_msg, UserWarning)
        
        self.components = parts

    def _check_input_format(self, ipt):
        pts = set(self.components)

        diff = pts.difference(ipt)
        if diff:
            raise ValueError('please inspect inputs - missing key(s): {}'.format(', '.join(map(repr, diff))))

    @property
    def required_input_keys(self) -> List[str]:
        return self.components

    @property
    def default_definition(self) -> List[tuple]:
        return self._default

    def forward(self, sequences: List[Dict[str, str]]) -> List[float]:
        self._check_input_format(sequences[0])

        out = []
        for part in self.components:
            sqs = glom(sequences, [part])
            component = self.__getattribute__(part)

            result = component(sqs)
            out.append(result)

        out = np.array(out).sum(axis=0)
        return out.tolist()
