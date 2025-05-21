__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"



class Hamming:

    def __init__(self, ):
        pass

    def general(self, s1, s2):
        return sum(i != j for i, j in zip(s1, s2))