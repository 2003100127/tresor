__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from phylotres.util.random.Ordering import Ordering as ranord
from phylotres.util.random.Sampling import Sampling as ranspl
from phylotres.util.random.Number import Number as rannum
from phylotres.pcr.Error import Error as pcrerr
from phylotres.util.Console import Console


class Amplify:

    def __init__(
            self,
            pcr_params,
    ):
        self.pcr_params = pcr_params
        self.console = Console()
        if 'verbose' in self.pcr_params.keys():
            self.console.verbose = self.pcr_params['verbose']
        else:
            self.console.verbose = True

    def np(self, ):
        for ipcr in range(self.pcr_params['pcr_num']):
            self.console.print('===>at PCR {}'.format(ipcr + 1))
            self.pcr_params['ipcr'] = ipcr
            self.console.print('===>Error assignment method: {}'.format(self.pcr_params['err_route']))
            if self.pcr_params['err_route'] == 'err2d':
                self.pcr_params = self.flow2D(params=self.pcr_params)
            elif self.pcr_params['err_route'] == 'tree':
                self.pcr_params = self.flowTree(params=self.pcr_params)
            elif self.pcr_params['err_route'] == 'minnow':
                self.pcr_params = self.flowMinnow(params=self.pcr_params)
            else:
                self.pcr_params = self.flow(params=self.pcr_params)
            # print(std_flow_params.keys())
        return self.pcr_params

    @pcrerr(method='default')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow(self, params):
        return params

    @pcrerr(method='err2d')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow2D(self, params):
        return params

    @pcrerr(method='minnow')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flowMinnow(self, params):
        return params

    @pcrerr(method='tree')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flowTree(self, params):
        return params