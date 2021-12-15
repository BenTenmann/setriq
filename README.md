# setriq: pairwise sequence distances
[![CircleCI](https://circleci.com/gh/BenTenmann/setriq/tree/main.svg?style=shield&circle-token=11d21cf82d1b29647f02543f6bfee9703a8f7bfe)](https://circleci.com/gh/BenTenmann/setriq/tree/main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![logo](fig/logo.png)

A `Python` package written in `C++` for computing pairwise distances between (immunoglobulin) sequences. 

## Install
This package is available on PyPI
```bash
pip install setriq
```

## Quickstart

`setriq` inherits from the `torch` philosophy of callable objects. Each `Metric` subclass is a callable upon 
initialisation, taking a list of objects (usually `str`) and returning a list of `float` values.

```python
import setriq
metric = setriq.CdrDist()

sequences = [
    'CASSLKPNTEAFF',
    'CASSAHIANYGYTF',
    'CASRGATETQYF'
]
distances = metric(sequences)
```

The returned list is flat and contains `N * (N - 1) / 2` elements, i.e. the lower (or upper) triangle of the distance 
matrix. To get the square form of the matrix, use `scipy.spatial.distance.squareform` on the returned distances.

## About

As the header suggests, `setriq` is a no-frills Python package for fast computation of pairwise sequence distances, with
a focus on immunoglobulins. It is a declarative framework and borrows many concepts from the popular `torch` library. It 
has been optimized for parallel compute on CPU architectures.

It can **only** perform pairwise, all-v-all distance computations. This decision was made to maximize consistency. 
