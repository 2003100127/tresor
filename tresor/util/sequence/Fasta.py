__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import gzip
from Bio import SeqIO


class Fasta:

    def __init__(self):
        pass

    def get_from_gz(self, fasta_fpn):
        arr = []
        with gzip.open(fasta_fpn, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                arr.append([str(record.id), str(record.seq)])
            # print(arr)
        return arr

    def save_(self, list_2d, sv_fp):
        for i, e in enumerate(list_2d):
            prot_name = str(e[0])
            seq = str(e[1])
            print('No.{} saving {} in FASTA format.'.format(i+1, prot_name))
            f = open(sv_fp + prot_name + '.fasta', 'w')
            f.write('>' + prot_name + '\n')
            f.write(seq + '\n')
            f.close()
        return 0