__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from phylotres.sequencing.Error import Error as seqerr
from phylotres.util.random.Sampling import Sampling as ranspl


class Calling:

    def __init__(self, seq_params):
        self.seq_params = seq_params

    def np(self, ):
        if self.seq_params['seq_sub_spl_number']:
            self.seq_params['spl_num'] = self.seq_params['seq_sub_spl_number']
        else:
            self.seq_params['spl_num'] = int(self.seq_params['data'].shape[0] * self.seq_params['seq_sub_spl_rate'])
        res_flow = self.flow(params=self.seq_params)
        # print(res_flow)
        return res_flow

    @seqerr(method='default')
    @ranspl(method='uniform')
    def flow(self, params):
        return params