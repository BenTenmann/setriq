"""
Substitution matrix convenience interface.

"""

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
    r"""
    The SubstitutionMatrix abstract base class. It holds convenience methods for loading and using the classic
    biological sequence substitution matrices (e.g. `BLOSUM`).

    Attributes
    ----------
    index : Dict[str, int]
        a mapping of strings (amino acids) to integers (matrix index)
    substitution_matrix : List[List[float]]
        the substitution scoring ($N \times N$) matrix

    Methods
    -------
    from_json(self, file_path: Union[str, Path]) -> SubstitutionMatrix:
        load a substitution matrix from a json file

    Examples
    --------
    suppose we have a token index `idx` and a substitution matrix `scores`
    >>> idx = {'hello': 1, 'world': 2}
    >>> scores = [[1., -1.],
    ...           [-1., 1.]]
    >>> matrix = SubstitutionMatrix(index=idx, substitution_matrix=scores)

    here we can see that we can provide any arbitrary substitution matrix, but in general it is advised to use the
    pre-loaded BLOSUM matrices
    >>> [BLOSUM45, BLOSUM62, BLOSUM90]  # choose one of the following

    these are just instances of `SubstitutionMatrix`, initialised through `from_json`
    """
    _required_keys = (
        'index',
        'substitution_matrix',
    )
    index: Dict[str, int]
    substitution_matrix: List[List[float]]

    def __init__(self, index: Dict[str, int], substitution_matrix: List[List[float]], *args, **kwargs):
        """
        Construct an instance of a SubstitutionMatrix object.

        Parameters
        ----------
        index : Dict[str, int]
            a token index which maps the string tokens to their index in the matrix
        substitution_matrix : List[List[float]]
            the substitution scoring matrix
        args
        kwargs
        """
        self.index = index
        self.substitution_matrix = substitution_matrix

    @classmethod
    def from_json(cls, file_path: Union[str, pathlib.Path]) -> "SubstitutionMatrix":
        """
        Load a SubstitutionMatrix from a json file.

        Parameters
        ----------
        file_path : Union[str, Path]
            a path to a json file holding at least the token index and the substitution scoring matrix

        Returns
        -------
        substitution_matrix : SubstitutionMatrix
            returns an instance of the SubstitutionMatrix class, holding the values found at `file_path`

        Examples
        --------
        >>> sub_mat = SubstitutionMatrix.from_json('/path/to/file.json')

        """
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
        """
        Get the score for a given substitution.

        Parameters
        ----------
        a : str
            a token in the index which substitutes
        b : str
            another token in the index which is substituted

        Returns
        -------
        substitution_score : float
            returns the score of a given substitution

        Examples
        --------
        >>> sub_mat = BLOSUM62
        >>> sub_mat('A', 'R')
        ... -1

        """
        i = self.index[a]
        j = self.index[b]

        out = self.substitution_matrix[i][j]
        return out


# below we load the matrices which come with the package -- this exposes them to the user
# they will be used for default settings in a number of metrics
PKG_NAME = __name__.split('.')[0]
DATA_DIR = pathlib.Path(pkg.resource_filename(PKG_NAME, 'data/'))

BLOSUM45 = SubstitutionMatrix.from_json(DATA_DIR / 'blosum-45.json')
BLOSUM62 = SubstitutionMatrix.from_json(DATA_DIR / 'blosum-62.json')
BLOSUM90 = SubstitutionMatrix.from_json(DATA_DIR / 'blosum-90.json')
