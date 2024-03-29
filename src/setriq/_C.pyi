from typing import Dict, List, Sequence

def cdr_dist(
    sequences: Sequence[str],
    substitution_matrix: List[List[float]],
    index: Dict[str, int],
    gap_opening_penalty: float,
    gap_extension_penalty: float,
) -> List[float]: ...
def levenshtein(sequences: Sequence[str], extra_cost: float) -> List[float]: ...
def tcr_dist_component(
    sequences: Sequence[str],
    substitution_matrix: List[List[float]],
    index: Dict[str, int],
    gap_penalty: float,
    gap_symbol: str,
    weight: float,
) -> List[float]: ...
def hamming(sequences: Sequence[str], mismatch_score: float) -> List[float]: ...
def jaro(sequences: Sequence[str], jaro_weights: List[float]) -> List[float]: ...
def jaro_winkler(
    sequences: Sequence[str], p: float, max_l: int, jaro_weights: List[float]
) -> List[float]: ...
def longest_common_substring(sequences: Sequence[str]) -> List[float]: ...
def optimal_string_alignment(sequences: Sequence[str]) -> List[float]: ...
def cdr_dist_sd(
    a: str,
    b: str,
    substitution_matrix: List[List[float]],
    index: Dict[str, int],
    gap_opening_penalty: float,
    gap_extension_penalty: float,
) -> float: ...
def levenshtein_sd(a: str, b: str, extra_cost: float) -> float: ...
def tcr_dist_component_sd(
    a: str,
    b: str,
    substitution_matrix: List[List[float]],
    index: Dict[str, int],
    gap_penalty: float,
    gap_symbol: str,
    weight: float,
) -> float: ...
def hamming_sd(a: str, b: str, mismatch_score: float) -> float: ...
def jaro_sd(a: str, b: str, jaro_weights: List[float]) -> float: ...
def jaro_winkler_sd(
    a: str, b: str, p: float, max_l: int, jaro_weights: List[float]
) -> float: ...
def longest_common_substring_sd(a: str, b: str) -> float: ...
def optimal_string_alignment_sd(a: str, b: str) -> float: ...
