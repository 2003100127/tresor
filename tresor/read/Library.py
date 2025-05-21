__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


from functools import wraps


class Library:

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, deal):
        manage = self.kwargs['method']
        @wraps(deal)
        def build(ph, *args, **kwargs):
            res = deal(ph, **kwargs)
            if kwargs['is_sv'] is True:
                if manage == 'default':
                    with open(kwargs['lib_fpn'], 'a') as file:
                        file.write(res + "\n")
                elif manage == 'separate':
                    with open(kwargs['lib_fpn'], 'a') as file:
                        file.write(kwargs['res'] + "\n")
            return res
        return build