__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


import os


def root_dict():
    """
    ..  @description:
        -------------
        abs file path.

    :return:
    """
    ROOT_DICT = os.path.dirname(os.path.abspath(__file__))
    return ROOT_DICT


def to(path):
    """

    Parameters
    ----------
    path

    Returns
    -------

    """
    return os.path.join(
        root_dict(),
        path
    )