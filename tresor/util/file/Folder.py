__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


import os


class Folder:

    def osmkdir(self, DIRECTORY):
        if not os.path.exists(DIRECTORY):
            os.makedirs(DIRECTORY)
        return 0