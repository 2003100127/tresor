__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import gzip
from Bio import SeqIO
import pyfastx


class read:

    def __init__(self):
        pass

    def fromgz(self, fastq_path, fastq_name, method='pyfastx'):
        names = []
        seqs = []
        placeholders = []
        qualities = []
        if method == 'biopython':
            with gzip.open(fastq_path + fastq_name + '.fastq.gz', 'rt') as handle:
                for record in SeqIO.parse(handle, 'fastq'):
                    # print()
                    seqs.append(''.join(record.seq))
                    names.append(''.join(record.name))
            return names, seqs, placeholders, qualities
        elif method == 'pyfastx':
            fq = pyfastx.Fastx(fastq_path + fastq_name + '.fastq.gz')
            for name, seq, qual, comment in fq:
                seqs.append(''.join(seq))
                names.append(''.join(name))
            return names, seqs, placeholders, qualities