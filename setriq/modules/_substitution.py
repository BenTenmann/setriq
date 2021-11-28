import abc
import pathlib
import pkg_resources as pkg
from typing import (
    Dict,
    List,
    Union
)

import srsly

__all__ = [
    'SubstitutionMatrix',
    'BLOSUM45',
    'BLOSUM62',
    'BLOSUM90',
]


class SubstitutionMatrix(abc.ABC):
    _required_keys = (
        'index',
        'substitution_matrix',
    )
    index: Dict[str, int]
    substitution_matrix: List[List[float]]

    def __init__(self, index: Dict[str, int], substitution_matrix: List[List[float]], *args, **kwargs):
        self.index = index
        self.substitution_matrix = substitution_matrix

    @classmethod
    def from_json(cls, file_path: Union[str, pathlib.Path]) -> "SubstitutionMatrix":
        if isinstance(file_path, str):
            file_path = pathlib.Path(file_path)

        values = srsly.read_json(file_path)
        for key in cls._required_keys:
            if key not in values:
                raise ValueError(f'required key {repr(key)} not in provided file: {file_path}')

        model = cls(**values)
        return model

    def __getitem__(self, key: str):
        out = self.__getattribute__(key)
        return out

    def keys(self):
        return self._required_keys

    def __call__(self, a: str, b: str) -> float:
        i = self.index[a]
        j = self.index[b]

        out = self.substitution_matrix[i][j]
        return out


PKG_NAME = __name__.split('.')[0]
DATA_DIR = pathlib.Path(pkg.resource_filename(PKG_NAME, 'data/'))

BLOSUM45 = SubstitutionMatrix.from_json(DATA_DIR / 'blosum-45.json')
BLOSUM62 = SubstitutionMatrix.from_json(DATA_DIR / 'blosum-62.json')
BLOSUM90 = SubstitutionMatrix.from_json(DATA_DIR / 'blosum-90.json')
