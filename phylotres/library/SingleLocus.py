__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import time
import numpy as np
from phylotres.util.random.Sampling import Sampling as ranspl
from phylotres.util.random.Number import Number as rannum
from phylotres.util.file.write.Writer import Writer as pfwriter
from phylotres.util.file.create.Folder import Folder as crtfolder
from phylotres.util.sequence.symbol.Single import Single as dnasgl
from phylotres.util.similarity.distance.Hamming import Hamming
from phylotres.read.umi.Design import Design as dumi
from phylotres.read.seq.Design import Design as dseq
from phylotres.read.primer.Design import Design as dprimer
from phylotres.read.adapter.Design import Design as dadapter
from phylotres.read.spacer.Design import Design as dspacer
from phylotres.util.Console import Console
from phylotres.util.sequence.Fasta import Fasta as sfasta
from phylotres.util.Kit import tactic6


class SingleLocus:

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

    def pooling(self, ):
        """

        Attributes
        ----------
        condi_map
            {'umi': ['alone', '1'], 'primer': ['alone', '1'], 'spacer': ['alone', '1'], 'adapter': ['alone', '1'], 'seq': ['alone', '2']}

        Returns
        -------

        """
        self.console.print("======>Sequencing library preparation starts")
        stime = time.time()
        sequencing_library = []
        umi_pool = []
        umi_cnt = 0

        ### +++++++++++++++ block: condition map +++++++++++++++
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
        self.console.print("======>Condition map: {}".format(condi_map))
        condi_keys = condi_map.keys()

        ### +++++++++++++++ block: select CDNA from a reference ome +++++++++++++++
        if self.fasta_cdna_fpn:
            self.console.print("======>Read CDNAs from a reference ome")
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
                    seed=i + 100000,
                    replace=True,
                )
                cdna_seqs_sel_maps[seq_i] = [seq_cdna_map[i] for i in cdna_ids_sel]
                # print(cdna_seqs_sel)
                self.pfwriter.generic(df=cdna_ids_sel, sv_fpn=self.working_dir + 'cdna_ids_' + seq_i + '.txt')
            del seq_cdna_map

        ### +++++++++++++++ block: generate each read +++++++++++++++
        for id in np.arange(self.seq_num):
            self.console.print("======>Read {} generation".format(id + 1))
            read_struct_ref = {}
            ### +++++++++++++++ block: generate umis +++++++++++++++
            if 'umi' in condi_keys:
                self.console.print("=========>UMI generation start")
                for umi_mark_id, umi_mark_suffix in enumerate(condi_map['umi']):
                    # print(umi_mark_id, id + self.permutation * self.seq_num + umi_cnt + umi_mark_id + 100000000)
                    umi_mark = '_' + umi_mark_suffix if umi_mark_suffix != 'alone' else ''
                    self.console.print("============>UMI condition {}: {}".format(umi_mark_id, 'umi' + umi_mark))
                    umi_flag = False
                    while not umi_flag:
                        umi_seed = id + self.permutation * self.seq_num + umi_cnt + (umi_mark_id + 1) * 100000000
                        umip = self.dumi(
                            dna_map=self.dna_map,
                            umi_unit_pattern=self.len_params['umi' + umi_mark]['umi_unit_pattern'],
                            pseudorandom_num=self.rannum.uniform(
                                low=0,
                                high=4,
                                num=self.len_params['umi' + umi_mark]['umi_unit_len'],
                                use_seed=self.is_seed,
                                seed=umi_seed,
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
                            read_struct_ref['umi' + umi_mark] = umi_i
                            umi_flag = True
                            umip.write(res=umi_i, lib_fpn=self.working_dir + 'umi' + umi_mark + '.txt', is_sv=self.is_sv_umi_lib)
                            umip.write(
                                res=str(umi_seed),
                                lib_fpn=self.working_dir + 'umi' + umi_mark + '_seeds.txt',
                                is_sv=self.is_sv_umi_lib,
                            )
                        else:
                            # print(id)
                            umi_cnt += 1

            ### +++++++++++++++ block: generate seqs +++++++++++++++
            if 'seq' in condi_keys:
                self.console.print("=========>Sequence generation start")
                for seq_mark_id, seq_mark_suffix in enumerate(condi_map['seq']):
                    seq_mark = '_' + seq_mark_suffix if seq_mark_suffix != 'alone' else ''
                    self.console.print("============>Sequence condition {}: {}".format(seq_mark_id, 'seq' + seq_mark))
                    if self.fasta_cdna_fpn:
                        seq_i = self.dseq(
                            cdna_seq=cdna_seqs_sel_maps[seq_mark_suffix][id],
                        ).cdna(lib_fpn=self.working_dir + 'seq' + seq_mark + '.txt', is_sv=self.is_sv_seq_lib)
                    else:
                        seq_seed = id + self.permutation * self.seq_num + 8000000 + (seq_mark_id+1) * 200000000
                        pseq = self.dseq(
                            dna_map=self.dna_map,
                            pseudorandom_num=self.rannum.uniform(
                                low=0,
                                high=4,
                                num=self.len_params['seq' + seq_mark],
                                use_seed=self.is_seed,
                                seed=seq_seed,
                            ),
                        )
                        seq_i = pseq.general(lib_fpn=self.working_dir + 'seq' + seq_mark + '.txt', is_sv=self.is_sv_seq_lib)
                        pseq.write(
                            res=str(seq_seed),
                            lib_fpn=self.working_dir + 'seq' + seq_mark + '_seeds.txt',
                            is_sv=self.is_sv_seq_lib,
                        )
                    read_struct_ref['seq' + seq_mark] = seq_i

            ### +++++++++++++++ block: generate primers +++++++++++++++
            if 'primer' in condi_keys:
                self.console.print("=========>Primer generation start")
                for primer_mark_id, primer_mark_suffix in enumerate(condi_map['primer']):
                    primer_mark = '_' + primer_mark_suffix if primer_mark_suffix != 'alone' else ''
                    self.console.print("============>Primer condition {}: {}".format(primer_mark_id, 'primer' + primer_mark))
                    primer_seed = id + self.permutation * self.seq_num + 8000000 + (primer_mark_id + 1) * 300000000
                    pprimer = self.dprimer(
                        dna_map=self.dna_map,
                        pseudorandom_num=self.rannum.uniform(
                            low=0,
                            high=4,
                            num=self.len_params['primer' + primer_mark],
                            use_seed=self.is_seed,
                            seed=primer_seed,
                        ),
                    )
                    primer_i = pprimer.general(lib_fpn=self.working_dir + 'primer' + primer_mark + '.txt', is_sv=self.is_sv_primer_lib)
                    pprimer.write(
                        res=str(primer_seed),
                        lib_fpn=self.working_dir + 'primer' + primer_mark + '_seeds.txt',
                        is_sv=self.is_sv_primer_lib,
                    )
                    read_struct_ref['primer' + primer_mark] = primer_i

            ### +++++++++++++++ block: generate adapters +++++++++++++++
            if 'adapter' in condi_keys:
                self.console.print("=========>Adapter generation start")
                for adapter_mark_id, adapter_mark_suffix in enumerate(condi_map['adapter']):
                    adapter_mark = '_' + adapter_mark_suffix if adapter_mark_suffix != 'alone' else ''
                    self.console.print("============>Adapter condition {}: {}".format(adapter_mark_id, 'adapter' + adapter_mark))
                    adapter_seed = id + self.permutation * self.seq_num + 8000000 + (adapter_mark_id+1) * 400000000
                    padapter = self.dadapter(
                        dna_map=self.dna_map,
                        pseudorandom_num=self.rannum.uniform(
                            low=0,
                            high=4,
                            num=self.len_params['adapter' + adapter_mark],
                            use_seed=self.is_seed,
                            seed=adapter_seed,
                        ),
                    )
                    adapter_i = padapter.general(lib_fpn=self.working_dir + 'adapter' + adapter_mark + '.txt', is_sv=self.is_sv_adapter_lib)
                    padapter.write(
                        res=str(adapter_seed),
                        lib_fpn=self.working_dir + 'adapter' + adapter_mark + '_seeds.txt',
                        is_sv=self.is_sv_adapter_lib,
                    )
                    read_struct_ref['adapter' + adapter_mark] = adapter_i

            ### +++++++++++++++ block: generate spacers +++++++++++++++
            if 'spacer' in condi_keys:
                self.console.print("=========>Spacer generation start")
                for spacer_mark_id, spacer_mark_suffix in enumerate(condi_map['spacer']):
                    spacer_mark = '_' + spacer_mark_suffix if spacer_mark_suffix != 'alone' else ''
                    self.console.print("============>Spacer condition {}: {}".format(spacer_mark_id, 'spacer' + spacer_mark))
                    spacer_seed = id + self.permutation * self.seq_num + 8000000 + (spacer_mark_id+1) * 500000000
                    pspacer = self.dspacer(
                        dna_map=self.dna_map,
                        pseudorandom_num=self.rannum.uniform(
                            low=0,
                            high=4,
                            num=self.len_params['spacer' + spacer_mark],
                            use_seed=self.is_seed,
                            seed=spacer_seed,
                        ),
                    )
                    spacer_i = pspacer.general(lib_fpn=self.working_dir + 'spacer' + spacer_mark + '.txt', is_sv=self.is_sv_spacer_lib)
                    pspacer.write(
                        res=str(spacer_seed),
                        lib_fpn=self.working_dir + 'spacer' + spacer_mark + '_seeds.txt',
                        is_sv=self.is_sv_spacer_lib,
                    )
                    read_struct_ref['spacer' + spacer_mark] = spacer_i

            read_struct_pfd_order = {condi: read_struct_ref[condi] for condi in self.condis}
            sequencing_library.append([self.paste([*read_struct_pfd_order.values()]), str(id), 'init'])
        # print(umi_cnt)
        # print(umi_pool)
        self.pfwriter.generic(df=sequencing_library, sv_fpn=self.working_dir + 'sequencing_library.txt')
        etime = time.time()
        self.console.print("===>Time for sequencing library preparation: {:.3f}s".format(etime-stime))
        return sequencing_library

    def paste(self, read_struct=[]):
        return ''.join(read_struct)


if __name__ == "__main__":
    from phylotres.path import to
    # print(DEFINE['cand_pool_fpn'])
    p = SingleLocus(
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
            'seq': 100,
            'seq_2': 100,
            'adapter': 10,
            'adapter_1': 10,
            'primer': 10,
            'primer_1': 10,
            'spacer': 10,
            'spacer_1': 10,
        },
        is_seed=True,

        working_dir=to('data/simu/'),
        fasta_cdna_fpn=False,
        # fasta_cdna_fpn=to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),

        # condis=['umi'],
        # condis=['umi', 'seq'],
        condis=['umi', 'primer', 'primer_1', 'spacer', 'spacer_1', 'adapter', 'adapter_1', 'seq', 'seq_2', 'umi_1'],
        sim_thres=3,
        permutation=0,
        verbose=False,
    )

    # print(p.umi_len)
    p.pooling()