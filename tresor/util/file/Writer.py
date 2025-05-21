__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


import pandas as pd
from functools import wraps


class Writer:

    def __init__(self, ):
        pass

    def __call__(self, deal):
        generic = self.generic
        @wraps(deal)
        def write(ph, *args, **kwargs):
            res = deal(ph, **kwargs)
            keys = [*kwargs.keys()]
            if kwargs['type'] == 'generic':
                generic(
                    df=kwargs['df'],
                    sv_fpn=kwargs['sv_fpn'],
                    df_sep='\t' if 'df_sep' not in keys else kwargs['df_sep'],
                    id_from=0 if 'id_from' not in keys else kwargs['id_from'],
                    header=None if 'header' not in keys else kwargs['header'],
                    index=False if 'index' not in keys else kwargs['index'],
                )
            return res
        return write

    def generic(self, df, sv_fpn, df_sep='\t', header=None, index=False, id_from=0):
        df_ = pd.DataFrame(df)
        # df_.index = df_.index + id_from
        return df_.to_csv(
            sv_fpn,
            sep=df_sep,
            header=header,
            index=index
        )