__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"



class Single:

    def __init__(self, ):
        pass

    def trim(action):
        def tube(deal):
            def wrapper(self, **kwargs):
                if action == 'leave_one_out':
                    res = deal(self, kwargs['ele_loo'])
                    # print(kwargs['ele_loo'])
                    res.remove(kwargs['ele_loo'])
                    return res
                elif action == 'leave_list_out':
                    import numpy as np
                    res = deal(self, kwargs['lis_loo'])
                    ref = np.arange(len(res))
                    c = list(set(ref).difference(set(kwargs['lis_loo'])))
                    # print(c)
                    res_ = np.array(res)
                    res_ = res_[c]
                    return res_
                else:
                    res = deal(self)
                    return res
            return wrapper
        return tube

    def _get(self, gap=False, universal=False):
        if universal:
            if gap:
                return ['A', 'C', 'G', 'T', 'B', '-']
            else:
                return ['A', 'C', 'G', 'T', 'B']
        else:
            if gap:
                return ['A', 'T', 'C', 'G', 'B', '-']
            else:
                return ['A', 'T', 'C', 'G', 'B']

    @trim(action='normal')
    def get(self, gap=False, universal=False):
        return self._get(gap=gap, universal=universal)

    @trim(action='leave_one_out')
    def getEleTrimmed(self, ele_loo, gap=False, universal=False):
        return self._get(gap=gap, universal=universal)

    @trim(action='leave_list_out')
    def getLisTrimmed(self, lis_loo=[], gap=False, universal=False):
        return self._get(gap=gap, universal=universal)

    def todict(self, nucleotides, reverse=False):
        aa_dict = {}
        for k, v in enumerate(nucleotides):
            aa_dict[v] = k
        if reverse:
            aa_dict = {v: k for k, v in aa_dict.items()}
        return aa_dict

    trim = staticmethod(trim)


if __name__ == "__main__":
    p = Single()
    # bs = p.get()
    bs = p.getEleTrimmed(ele_loo='A')
    # bs = p.getLisTrimmed(lis_loo=[0, 1])
    print(bs)
    print(p.todict(bs, reverse=True))