# Quantum arrival time
A Python package for modelling and analyzing arrival time distributions of superpositions of Gaussian wave packets. Currently, it supports the <i>quantum clock</i> (Page-Wootters) and the <i>quantum flux</i> arrival time proposals. Kijowski's and Leavens' distributions are in the pipeline.

[![arXiv](https://img.shields.io/badge/arXiv-2205.02219-b31b1b.svg)](https://doi.org/10.48550/arXiv.2205.02219)

Contributors: Simone Roncallo [@simoneroncallo](https://github.com/simoneroncallo) <br>
Reference: Simone Roncallo, Krzysztof Sacha and Lorenzo Maccone <i>“When does a particle arrive?”</i> [Quantum 7, 968 (2023)](https://doi.org/10.22331/q-2023-03-30-968)

## Installation
The package can be downloaded and set up in a [conda](https://docs.conda.io/) environment with
```bash
conda create --name arrival-env
conda activate arrival-env
conda install --file requirements.txt
```

## Example
Consider a single Gaussian wave packet with parameters 
```python
X = [-10] # Initial position
P = [7] # Initial momentum
S = [1] # Initial width (standard deviation)
M = [1] # Mass

```
The arrival time distribution can be computed as
```python
from toa import GaussianTrain, ArrivalTime
train = GaussianTrain(X, P, S, M)
toa = ArrivalTime(train)
toa.visualize(numPoints, tLim, xLim, 0)
```
with `numPoints` the number of samples and `tLim, xLim` the temporal and spatial domain, respectively.

## Structure
This repository has the following structure
```bash
toa/
├── distributions.py # Arrival time distributions
└── gaussian.py # Gaussian wave packets

example.ipynb # Example notebook
requirements.txt # Dependencies
```

