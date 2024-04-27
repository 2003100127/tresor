__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from tresor.pcr.Error import Error as pcrerr
from tresor.util.random.Ordering import Ordering as ranord
from tresor.util.random.Sampling import Sampling as ranspl
from tresor.util.random.Number import Number as rannum
from tresor.util.Console import Console


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
            if self.pcr_params['err_route'] == 'err1d':
                self.pcr_params = self.flow1D(params=self.pcr_params)
            elif self.pcr_params['err_route'] == 'err2d':
                self.pcr_params = self.flow2D(params=self.pcr_params)
            elif self.pcr_params['err_route'] == 'bftree':
                self.pcr_params = self.flow_bftree(params=self.pcr_params)
            elif self.pcr_params['err_route'] == 'sptree':
                self.pcr_params = self.flow_sptree(params=self.pcr_params)
            elif self.pcr_params['err_route'] == 'mutation_table_complete':
                self.pcr_params = self.flow_mutation_table_complete(params=self.pcr_params)
            elif self.pcr_params['err_route'] == 'mutation_table_minimum':
                self.pcr_params = self.flow_mutation_table_minimum(params=self.pcr_params)
            else:
                self.pcr_params = self.flow_bftree(params=self.pcr_params)
            # print(std_flow_params.keys())
        return self.pcr_params

    @pcrerr(method='err1d')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow1D(self, params):
        return params

    @pcrerr(method='err2d')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow2D(self, params):
        return params

    @pcrerr(method='mutation_table_complete')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow_mutation_table_complete(self, params):
        return params

    @pcrerr(method='mutation_table_minimum')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow_mutation_table_minimum(self, params):
        return params

    @pcrerr(method='sptree')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow_sptree(self, params):
        return params

    @pcrerr(method='bftree')
    @ranspl(method='uniform')
    @rannum(type='binomial')
    @ranord(method='uniform')
    def flow_bftree(self, params):
        return params