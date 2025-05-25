## System Requirement

There is no requirement for Tresor. It is a cross-platform computational package.

## PyPI (**recommended**)

!!! info "Note"

    Please make sure to use the latest version of tresor, as earlier versions may contain bugs. If you do not include the `--upgrade` flag during installation, you might encounter issues that prevent tresor from functioning properly in your analysis.

[tresor homepage](https://pypi.org/project/tresor/)

```shell
# create a conda environment
conda create --name tresor python=3.11

# activate the conda environment
conda activate tresor

# the latest version
pip install tresor --upgrade
```

## Conda

[tresor homepage on Anaconda](https://anaconda.org/Jianfeng_Sun/tresor)

```shell
# create a conda environment
conda create --name tresor python=3.11

# activate the conda environment
conda activate tresor

# the latest version
conda install jianfeng_sun::tresor
```


## Docker

[tresor homepage on Docker](https://hub.docker.com/r/2003100127/tresor)

You can first choose which type of operating system (OS) you would like to install Docker software. Please refer to [https://docs.docker.com/engine/install](https://docs.docker.com/engine/install). For example, if your computational work is based on a Windows OS, you can choose to install a Desktop version of Docker. Please refer to [https://docs.docker.com/desktop/install/windows-install](https://docs.docker.com/desktop/install/windows-install).

```shell
docker pull 2003100127/tresor
```


## Github

[tresor homepage on Github](https://github.com/2003100127/tresor)

```shell
# create a conda environment
conda create --name tresor python=3.11

# activate the conda environment
conda activate tresor

# create a folder
mkdir project

# go to the folder
cd project

# fetch Tresor repository with the latest version
git clone https://github.com/2003100127/tresor.git

# enter this repository
cd tresor

# do the following command
pip install .
# or
python setup.py install
```
