<h1 align="center">
    <img src="https://github.com/cribbslab/tresor-dev/blob/main/img/Tresor-logo.png?raw=true" width="276" height="114">
    <br>
</h1>

![](https://img.shields.io/pypi/v/tresor?logo=PyPI)
![](https://img.shields.io/badge/last_released-Jul14._2023-green.svg)
![](https://img.shields.io/badge/tresor-executable-519dd9.svg)


<!-- ![Build](https://github.com/2003100127/tresor/actions/workflows/build.yml/badge.svg) -->

###### tags: `computational biology`, `sequencing read simulation`

## Overview

```angular2html

```

Tresor is a Python toolkit for simulating sequencing reads at the single-locus, bulk RNA-seq, and single-cell levels. It is implemented based on phylogenetic tree-based methods, which allows for ultra-fast simulation read generation. Tresor allows both short-read and long-read sequencing read simulation, and substitution and indel (insertion and deletion) errors added to reads. Tresor implements a very flexible read generation framework, which allows users to design their simulated reads in any forms and structures. Tresor can vastly help both computational and experimental researchers to swiftly test their sequencing method ideas.

## Documentation

Please check how to use the full functionalities of Tresor in the documentation https://cribbslab.github.io/tresor.

## Installation

### Using pip (recommended)

```sh
# create a conda environment
conda create --name tresor python=3.11

# activate the conda environment
conda activate tresor

# the latest version
pip install tresor --upgrade
```

## Citation

Please cite our work if you use Tresor in your research.
```angular2html
@article{tresor,
    title = {Tresor: high-performance simulation of sequencing reads},
    author = {Jianfeng Sun and Adam P. Cribbs},
    url = {https://github.com/cribbslab/tresor},
    year = {2023},
}
```

## Contact

Please report any questions on [issue](https://github.com/2003100127/tresor/issues) pages.
