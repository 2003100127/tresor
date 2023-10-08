__version__ = "v1.0"
__copyright__ = "Copyright 2023"
__license__ = "MIT"
__lab__ = "cribbslab"

import time
import numpy as np
from phylotres.util.random.Sampling import Sampling as ranspl
from phylotres.util.random.Number import Number as rannum
from phylotres.util.file.write.Writer import Writer as pfwriter
from phylotres.util.file.create.Folder import Folder as crtfolder
from phylotres.util.sequence.symbol.Single import Single as dnasgl
from phylotres.util.similarity.distance.Hamming import hamming
from phylotres.read.umi.Design import Design as dumi
from phylotres.read.seq.Design import Design as dseq
from phylotres.read.primer.Design import Design as dprimer
from phylotres.read.adapter.Design import Design as dadapter
from phylotres.read.spacer.Design import Design as dspacer
from phylotres.util.Console import Console
from phylotres.util.sequence.Fasta import Fasta as sfasta
from phylotres.util.Kit import tactic6


class general:

    def __init__(
            self,
            len_params,
            fasta_cdna_fpn,
            seq_num=50,
            is_seed=False,

            working_dir='./simu/',
            condis=['umi'],
            sim_thres=3,
            permutation=0,
            is_sv_umi_lib=True,
            is_sv_seq_lib=True,
            is_sv_primer_lib=True,
            is_sv_adapter_lib=True,
            is_sv_spacer_lib=True,
            verbose=True,
    ):
        self.pfwriter = pfwriter()
        self.ranspl = ranspl()
        self.rannum = rannum()
        self.dnasgl = dnasgl()
        self.crtfolder = crtfolder()
        self.dumi = dumi
        self.dseq = dseq
        self.dprimer = dprimer
        self.dadapter = dadapter
        self.dspacer = dspacer
        self.working_dir = working_dir
        self.len_params = len_params
        self.is_seed = is_seed
        self.fasta_cdna_fpn = fasta_cdna_fpn
        self.is_sv_umi_lib = is_sv_umi_lib
        self.is_sv_seq_lib = is_sv_seq_lib
        self.is_sv_primer_lib = is_sv_primer_lib
        self.is_sv_adapter_lib = is_sv_adapter_lib
        self.is_sv_spacer_lib = is_sv_spacer_lib
        self.seq_num = seq_num
        self.condis = condis
        self.sim_thres = sim_thres
        self.permutation = permutation
        self.dna_map = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=True)
        self.crtfolder.osmkdir(working_dir)

        self.console = Console()
        self.console.verbose = verbose

    def pooling(self,):
        stime = time.time()
        seqs = []
        umi_pool = []
        umi_cnt = 0

        condi_map = {}
        for condi in self.condis:
            condi_arr = condi.split("_")
            condi_map[condi_arr[0]] = []
        for condi in self.condis:
            condi_arr = condi.split("_")
            if len(condi_arr) == 1:
                condi_map[condi_arr[0]].append('alone')
            else:
                condi_map[condi_arr[0]].append(condi_arr[1])
        print(condi_map)
        condi_keys = condi_map.keys()

        cdna_seqs_sel_maps = {}
        seq_cdna_map = tactic6(
            arr_2d=sfasta().get_from_gz(
                self.fasta_cdna_fpn,
            )
        )
        cdna_ids = [*seq_cdna_map.keys()]
        for i, seq_i in enumerate(condi_map['seq']):
            cdna_ids_sel = self.ranspl.uniform(
                data=cdna_ids,
                num=self.seq_num,
                use_seed=self.is_seed,
                seed=i + 10000,
                replace=True,
            )
            cdna_seqs_sel_maps[seq_i] = [seq_cdna_map[i] for i in cdna_ids_sel]
            # print(cdna_seqs_sel)
            self.pfwriter.generic(df=cdna_ids_sel, sv_fpn=self.working_dir + 'cdna_ids_' + seq_i + '.txt')
        del seq_cdna_map

        for id in np.arange(self.seq_num):
            read_struct_ref = {}
            if 'umi' in condi_keys:
                for umi_mark_id, umi_mark_suffix in enumerate(condi_map['umi']):
                    umi_mark = '_' + umi_mark_suffix if umi_mark_suffix != 'alone' else ''
                    umi_flag = False
                    while not umi_flag:
                        umip = self.dumi(
                            dna_map=self.dna_map,
                            umi_unit_pattern=self.len_params['umi' + umi_mark]['umi_unit_pattern'],
                            pseudorandom_num=self.rannum.uniform(
                                low=0,
                                high=4,
                                num=self.len_params['umi' + umi_mark]['umi_unit_len'],
                                use_seed=self.is_seed,
                                seed=id + self.permutation * self.seq_num + umi_cnt + umi_mark_id*1000000,
                            ),
                        )
                        umi_i = umip.reoccur(is_sv=False)
                        edh = np.array([hamming().general(umi_i, j) for j in umi_pool])
                        # for j in umi_pool:
                        #     if hamming().general(umi_i, j) < self.sim_thres:
                        #         print(umi_i, j)
                        if len(edh[edh < self.sim_thres]) == 0:
                            # print(len(edh[edh < self.sim_thres]))
                            umi_pool.append(umi_i)
                            read_struct_ref['umi' + umi_mark] = umi_i
                            umi_flag = True
                            umip.write(res=umi_i, lib_fpn=self.working_dir + 'umi' + umi_mark + '.txt', is_sv=self.is_sv_umi_lib)
                        else:
                            # print(id)
                            umi_cnt += 1

            if 'seq' in condi_keys:
                for _, seq_mark_suffix in enumerate(condi_map['seq']):
                    seq_mark = '_' + seq_mark_suffix if seq_mark_suffix != 'alone' else ''
                    seq_i = self.dseq(
                        cdna_seq=cdna_seqs_sel_maps[seq_mark_suffix][id],
                    ).cdna(lib_fpn=self.working_dir + 'seq' + seq_mark + '.txt', is_sv=self.is_sv_seq_lib)
                    read_struct_ref['seq' + seq_mark] = seq_i

            if 'primer' in condi_keys:
                for primer_mark_id, primer_mark_suffix in enumerate(condi_map['primer']):
                    primer_mark = '_' + primer_mark_suffix if primer_mark_suffix != 'alone' else ''
                    primer_i = self.dprimer(
                        dna_map=self.dna_map,
                        pseudorandom_num=self.rannum.uniform(
                            low=0,
                            high=4,
                            num=self.len_params['primer' + primer_mark],
                            use_seed=self.is_seed,
                            seed=id + self.permutation * self.seq_num + 8000000 + primer_mark_id*1000000,
                        ),
                    ).general(lib_fpn=self.working_dir + 'primer' + primer_mark + '.txt', is_sv=self.is_sv_primer_lib)
                    read_struct_ref['primer' + primer_mark] = primer_i

            if 'adapter' in condi_keys:
                for adapter_mark_id, adapter_mark_suffix in enumerate(condi_map['adapter']):
                    adapter_mark = '_' + adapter_mark_suffix if adapter_mark_suffix != 'alone' else ''
                    adapter_i = self.dadapter(
                        dna_map=self.dna_map,
                        pseudorandom_num=self.rannum.uniform(
                            low=0,
                            high=4,
                            num=self.len_params['adapter' + adapter_mark],
                            use_seed=self.is_seed,
                            seed=id + self.permutation * self.seq_num + 8000000 + adapter_mark_id*1000000,
                        ),
                    ).general(lib_fpn=self.working_dir + 'adapter' + adapter_mark + '.txt', is_sv=self.is_sv_adapter_lib)
                    read_struct_ref['adapter' + adapter_mark] = adapter_i

            if 'spacer' in condi_keys:
                for spacer_mark_id, spacer_mark_suffix in enumerate(condi_map['spacer']):
                    spacer_mark = '_' + spacer_mark_suffix if spacer_mark_suffix != 'alone' else ''
                    spacer_i = self.dspacer(
                        dna_map=self.dna_map,
                        pseudorandom_num=self.rannum.uniform(
                            low=0,
                            high=4,
                            num=self.len_params['spacer' + spacer_mark],
                            use_seed=self.is_seed,
                            seed=id + self.permutation * self.seq_num + 8000000 + spacer_mark_id*1000000,
                        ),
                    ).general(lib_fpn=self.working_dir + 'spacer' + spacer_mark + '.txt', is_sv=self.is_sv_spacer_lib)
                    read_struct_ref['spacer' + spacer_mark] = spacer_i

            read_struct_pfd_order = {condi: read_struct_ref[condi] for condi in self.condis}
            seqs.append([self.paste([*read_struct_pfd_order.values()]), str(id), 'init'])
        # print(umi_cnt)
        # print(umi_pool)
        etime = time.time()
        self.console.print("===>time for generating initial pool of sequences: {:.3f}s".format(etime-stime))
        return seqs

    def paste(self, read_struct=[]):
        return ''.join(read_struct)


if __name__ == "__main__":
    from phylotres.path import to
    # print(DEFINE['cand_pool_fpn'])
    p = general(
        seq_num=50,

        len_params={
            'umi': {
                'umi_unit_pattern': 3,
                'umi_unit_len': 12,
            },
            'umi_1': {
                'umi_unit_pattern': 3,
                'umi_unit_len': 12,
            },
            'adapter': 10,
            'adapter_1': 10,
            'primer': 10,
            'primer_1': 10,
            'spacer': 10,
            'spacer_1': 10,
        },
        is_seed=True,
        is_sv_umi_lib=True,
        is_sv_seq_lib=True,

        working_dir=to('data/simu/'),
        fasta_cdna_fpn=to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),

        # condis=['umi'],
        # condis=['umi', 'seq'],
        condis=['umi', 'umi_1', 'primer', 'primer_1', 'spacer', 'spacer_1', 'adapter', 'adapter_1', 'seq', 'seq_2'],
        sim_thres=3,
        permutation=0,
    )

    # print(p.umi_len)
    res = p.pooling()
    print(res)