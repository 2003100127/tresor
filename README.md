<h1 align="center">
    <img src="https://github.com/2003100127/tresor/blob/main/img/Tresor-logo.png?raw=true" width="276" height="114">
    <br>
</h1>

[![Anaconda-Server Badge](https://anaconda.org/jianfeng_sun/tresor/badges/latest_release_date.svg)](https://anaconda.org/jianfeng_sun/tresor)
![PyPI](https://img.shields.io/pypi/v/tresor?logo=PyPI)
![Docker Image Version (latest)](https://img.shields.io/docker/v/2003100127/tresor)
![Docker Pulls](https://img.shields.io/docker/pulls/2003100127/tresor)
[![Anaconda-Server Badge](https://anaconda.org/jianfeng_sun/tresor/badges/version.svg)](https://anaconda.org/jianfeng_sun/tresor)
![](https://img.shields.io/docker/automated/2003100127/tresor.svg)
![](https://img.shields.io/github/stars/2003100127/tresor?logo=GitHub&color=blue)
[![Documentation Status](https://readthedocs.org/projects/tresor/badge/?version=latest)](https://tresor.readthedocs.io/en/latest/?badge=latest)
[![Downloads](https://pepy.tech/badge/tresor)](https://pepy.tech/project/tresor)

<hr>

#### Platform

![Python](https://img.shields.io/badge/-Python-000?&logo=Python)
![Docker](https://img.shields.io/badge/-Docker-000?&logo=Docker)
![Anaconda](https://img.shields.io/badge/-Anaconda-000?&logo=Anaconda)
![PyPI](https://img.shields.io/badge/-PyPI-000?&logo=PyPI)

###### tags: `computational biology`, `sequencing read simulation`

## Overview

```angular2html

```

Tresor is a Python toolkit for simulating sequencing reads at the single-locus, bulk RNA-seq, and single-cell levels. It is implemented based on phylogenetic tree-based methods, which allows for ultra-fast simulation read generation. Tresor allows both short-read and long-read sequencing read simulation, and substitution and indel (insertion and deletion) errors added to reads. Tresor implements a very flexible read generation framework, which allows users to design their simulated reads in any forms and structures. Tresor can vastly help both computational and experimental researchers to swiftly test their sequencing method ideas.

## Documentation

Please check how to use the full functionalities of Tresor in the documentation https://2003100127.github.io/tresor.

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
    url = {https://github.com/2003100127/tresor},
    year = {2023},
}
```

## Contact

Please report any questions on [issue](https://github.com/2003100127/tresor/issues) pages.
