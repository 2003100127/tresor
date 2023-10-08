__version__ = "v1.0"
__copyright__ = "Copyright 2023"
__license__ = "MIT"
__lab__ = "cribbslab"

from phylotres.read.inf.Pseudo import Pseudo as seqpseudo
from phylotres.read.Library import Library as liblogginger


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