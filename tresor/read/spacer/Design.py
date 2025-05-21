__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


from tresor.read.inf.Pseudo import Pseudo as seqpseudo
from tresor.read.Library import Library as liblogginger


class Design(seqpseudo):

    def __init__(self, *args, **kwargs):
        super(Design, self).__init__(*args, **kwargs)
        self.args = args
        self.kwargs = kwargs

    @liblogginger(method='default')
    def general(self, lib_fpn='./spacer.txt', is_sv=True):
        return ''.join([
            self.kwargs['dna_map'][i] for i in
            self.kwargs['pseudorandom_num']
        ])

    @liblogginger(method='separate')
    def write(self, **kwargs):
        return 'written'