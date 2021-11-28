import abc
import warnings
from typing import Dict, List

import numpy as np
from glom import glom

import setriq._C as C
from ._substitution import *

__all__ = [
    'CdrDist',
]


class Metric(abc.ABC):

    @abc.abstractmethod
    def forward(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        out = self.forward(*args, **kwargs)

        return out


class CdrDist(Metric):
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
        
        self._parts = parts

    def _check_input_format(self, ipt):
        pts = set(self._parts)

        diff = pts.difference(ipt)
        if diff:
            raise ValueError('please inspect inputs - missing key(s): {}'.format(', '.join(map(repr, diff))))

    @property
    def required_input_keys(self):
        return self._parts

    @property
    def default_definition(self):
        return self._default

    def forward(self, sequences: List[Dict[str, str]]) -> List[float]:
        self._check_input_format(sequences[0])

        out = []
        for part in self._parts:
            sqs = glom(sequences, [part])
            component = self.__getattribute__(part)

            result = component(sqs)
            out.append(result)

        out = np.array(out).sum(axis=0)
        return out.tolist()
