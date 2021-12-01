# setriq: pairwise sequence distances
[![CircleCI](https://circleci.com/gh/BenTenmann/setriq/tree/main.svg?style=shield&circle-token=11d21cf82d1b29647f02543f6bfee9703a8f7bfe)](https://circleci.com/gh/BenTenmann/setriq/tree/main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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
