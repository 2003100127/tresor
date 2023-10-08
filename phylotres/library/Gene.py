__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import time
import numpy as np
from scipy.sparse import coo_matrix
from phylotres.util.random.Sampling import Sampling as ranspl
from phylotres.util.random.Number import Number as rannum
from phylotres.util.file.read.Reader import Reader as pfreader
from phylotres.util.file.create.Folder import Folder as crtfolder
from phylotres.util.sequence.symbol.Single import Single as dnasgl
from phylotres.read.barcode.Design import Design as dbc
from phylotres.read.umi.Design import Design as dumi
from phylotres.read.seq.Design import Design as dseq
from phylotres.util.similarity.distance.Hamming import Hamming


class Bulk:

    def __init__(
            self,
            gspl,
            is_seed=False,
            umi_unit_pattern=3,
            umi_unit_len=12,
            seq_len=100,
            is_sv_umi_lib=True,
            is_sv_seq_lib=True,
            umi_lib_fpn='./umi.txt',
            seq_lib_fpn='./seq.txt',
            working_dir='./simu/',
            condis=['umi'],
            sim_thres=2,
            permutation=0,
    ):
        self.pfreader = pfreader()
        self.ranspl = ranspl()
        self.rannum = rannum()
        self.dnasgl = dnasgl()
        self.crtfolder = crtfolder()
        self.dumi = dumi
        self.dbc = dbc
        self.dseq = dseq
        self.is_seed = is_seed
        self.is_sv_umi_lib = is_sv_umi_lib
        self.umi_lib_fpn = umi_lib_fpn
        self.is_sv_seq_lib = is_sv_seq_lib
        self.seq_lib_fpn = seq_lib_fpn
        self.umi_unit_pattern = umi_unit_pattern
        self.umi_unit_len = umi_unit_len
        self.seq_len = seq_len
        self.condis = condis
        self.sim_thres = sim_thres
        self.permutation = permutation
        self.dna_map = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=True)
        self.crtfolder.osmkdir(working_dir)

        self.gspl = gspl
        self.gene_map = {k: v for k, v in enumerate(self.gspl.index)}
        # print(self.gene_map)
        csr_ = coo_matrix(self.gspl)
        self.gspl_arr = np.transpose([
            csr_.row.tolist(),
            csr_.col.tolist(),
            csr_.data.tolist(),
        ]).astype(int)
        self.gspl_arr = self.gspl_arr[self.gspl_arr[:, 0] == 0]
        print(self.gspl_arr)

    @property
    def umi_len(self, ):
        return self.umi_unit_pattern * self.umi_unit_len

    def pooling(self,):
        stime = time.time()
        seqs = []
        umi_pool = []
        umi_cnt = 0
        for x, gs in enumerate(self.gspl_arr):
            sample = gs[0]
            gene = gs[1]
            seq_num = gs[2]
            for id in np.arange(seq_num):
                read_struct_ref = {}
                if 'umi' in self.condis:
                    umi_flag = False
                    while not umi_flag:
                        umip = self.dumi(
                            dna_map=self.dna_map,
                            umi_unit_pattern=self.umi_unit_pattern,
                            pseudorandom_num=self.rannum.uniform(
                                low=0,
                                high=4,
                                num=self.umi_unit_len,
                                use_seed=self.is_seed,
                                seed=id + self.permutation * seq_num + umi_cnt,
                            ),
                        )
                        umi_i = umip.reoccur(is_sv=False)
                        edh = np.array([Hamming().general(umi_i, j) for j in umi_pool])
                        # for j in umi_pool:
                        #     if hamming().general(umi_i, j) < self.sim_thres:
                        #         print(umi_i, j)
                        if len(edh[edh < self.sim_thres]) == 0:
                            # print(len(edh[edh < self.sim_thres]))
                            umi_pool.append(umi_i)
                            read_struct_ref['umi'] = umi_i
                            umi_flag = True
                            umip.write(
                                res=umi_i,
                                lib_fpn=self.umi_lib_fpn + 'umi' + '_s_' + str(sample) + '_g_' + str(gene) + '.txt',
                                is_sv=self.is_sv_umi_lib)
                        else:
                            # print(id)
                            umi_cnt += 1
                # if 'seq' in self.condis:
                #     seq_i = self.dseq(
                #         dna_map=self.dna_map,
                #         pseudorandom_num=self.rannum.uniform(
                #             low=0,
                #             high=4,
                #             num=self.seq_len,
                #             use_seed=self.is_seed,
                #             seed=id + self.permutation * seq_num + 5000000,
                #         ),
                #     ).general(lib_fpn=self.seq_lib_fpn, is_sv=self.is_sv_seq_lib)
                #     read_struct_ref['seq'] = seq_i
                read_struct_pfd_order = {condi: read_struct_ref[condi] for condi in self.condis}
                seqs.append(
                    [self.paste(read_struct=[*read_struct_pfd_order.values()]),
                    str(id) + '*s*' + str(sample) + '*g*' + str(gene) + '*',
                    'init'
                ])
        # print(umi_cnt)
        # print(umi_pool)
        etime = time.time()
        print("===>time for generating initial pool of sequences: {:.3f}s".format(etime-stime))
        return seqs

    def paste(self, read_struct=[]):
        return ''.join(read_struct)


if __name__ == "__main__":
    from phylotres.path import to
    DEFINE = {
        '': '',
    }
    # print(DEFINE['cand_pool_fpn'])

    from phylotres.gspl.FromSimulator import fromSimulator
    gspl = fromSimulator(simulator='SPsimSeq').run()

    p = Bulk(
        gspl=gspl,
        umi_unit_pattern=1,
        umi_unit_len=10,
        seq_len=100 - 10,
        is_seed=True,

        is_sv_umi_lib=True,
        umi_lib_fpn=to('data/simu/umi_seq/permute_1/'),
        is_sv_seq_lib=True,
        seq_lib_fpn=to('data/simu/umi_seq/permute_1/seq.txt'),
        working_dir=to('data/simu/umi_seq/permute_1/'),

        condis=['umi'],
        sim_thres=3,
        permutation=0,
    )

    # print(p.umi_len)
    res = p.pooling()
    print(res)
    # print(p.asda())