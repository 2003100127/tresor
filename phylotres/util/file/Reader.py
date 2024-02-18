__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import pandas as pd
from functools import wraps


class Reader(object):

    def __init__(self, ):
        pass

    def __call__(self, deal):
        generic = self.generic
        @wraps(deal)
        def read(ph, *args, **kwargs):
            deal(ph, **kwargs)
            keys = [*kwargs.keys()]
            if kwargs['type'] == 'generic':
                return generic(
                    df_fpn=kwargs['df_fpn'],
                    df_sep='\t' if 'df_sep' not in keys else kwargs['df_sep'],
                    skiprows=False if 'skiprows' not in keys else kwargs['skiprows'],
                    header=None if 'header' not in keys else kwargs['header'],
                    is_utf8=False if 'is_utf8' not in keys else kwargs['is_utf8'],
                )
        return read

    def generic(self, df_fpn, df_sep='\t', skiprows=None, header=None, is_utf8=False):
        if is_utf8:
            return pd.read_csv(
                df_fpn,
                sep=df_sep,
                header=header,
                encoding='utf-8',
                skiprows=skiprows,
            )
        else:
            return pd.read_csv(
                df_fpn,
                sep=df_sep,
                header=header,
                skiprows=skiprows
            )