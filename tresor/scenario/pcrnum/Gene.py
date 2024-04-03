__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import time
import numpy as np
from tresor.library.Gene import Gene as bulksimulib
from tresor.pcr.Amplify import Amplify as pcr
from tresor.pcr.Subsampling import Subsampling
from tresor.sequencing.Calling import Calling as seq
from tresor.util.sequence.fastq.Write import write as wfastq
from tresor.util.Console import Console


class Gene:

    def __init__(
            self,
            gspl,

            len_params,
            seq_num,
            is_sv_umi_lib,
            is_sv_seq_lib,
            is_sv_primer_lib,
            is_sv_adapter_lib,
            is_sv_spacer_lib,
            working_dir,
            condis,
            sim_thres,
            permutation,

            err_route,
            pcr_error,
            seq_error,
            pcr_nums,
            err_num_met,
            use_seed,
            seed,

            ampl_rate,
            sv_fastq_fp,

            seq_sub_spl_number=None,
            seq_sub_spl_rate=1/3,

            verbose=True,
            **kwargs,
    ):
        self.len_params = len_params
        self.seq_num = seq_num
        self.is_sv_umi_lib = is_sv_umi_lib
        self.is_sv_seq_lib = is_sv_seq_lib
        self.is_sv_primer_lib = is_sv_primer_lib
        self.is_sv_adapter_lib = is_sv_adapter_lib
        self.is_sv_spacer_lib = is_sv_spacer_lib

        self.working_dir = working_dir
        self.condis = condis
        self.sim_thres = sim_thres
        self.permutation = permutation

        self.gspl = gspl

        self.err_route = err_route
        self.ampl_rate = ampl_rate
        self.pcr_error = pcr_error
        self.seq_error = seq_error
        self.err_num_met = err_num_met
        self.use_seed = use_seed
        self.seed = seed

        self.pcr_nums = pcr_nums
        self.seq_sub_spl_number = seq_sub_spl_number
        self.seq_sub_spl_rate = seq_sub_spl_rate
        self.sv_fastq_fp = sv_fastq_fp

        self.pcr = pcr
        self.seq = seq
        self.wfastq = wfastq
        self.subsampling = Subsampling()

        self.kwargs = kwargs

        self.console = Console()
        self.console.verbose = verbose
        self.verbose = verbose

        ### +++++++++++++++ block: generate sequencing library +++++++++++++++
        self.console.print('===>Sequencing library generation starts')
        self.sequencing_library = bulksimulib(
            gspl=self.gspl,
            len_params=self.len_params,
            seq_num=self.seq_num,
            is_seed=self.use_seed,
            working_dir=self.working_dir,
            condis=self.condis,
            sim_thres=self.sim_thres,
            permutation=self.permutation,
            is_sv_umi_lib=self.is_sv_umi_lib,
            is_sv_seq_lib=self.is_sv_seq_lib,
            is_sv_primer_lib=self.is_sv_primer_lib,
            is_sv_adapter_lib=self.is_sv_adapter_lib,
            is_sv_spacer_lib=self.is_sv_spacer_lib,
            verbose=self.verbose,
            mode=self.kwargs['mode'],
            material_params=self.kwargs['material_params'] if 'material_params' in self.kwargs.keys() else None,
            seq_params=self.kwargs['seq_params'] if 'seq_params' in self.kwargs.keys() else None,
        ).pooling()
        ### self.sequencing_library
        # [['GGGGGGCCCTTTTTTTTTTTTAAAAAACCCAAAGGG', '0*s*0*g*3*', 'init'],
        # ['GGGTTTGGGAAAAAAAAAAAATTTTTTGGGCCCGGG', '1*s*0*g*3*', 'init'],
        # ['CCCTTTCCCGGGGGGTTTGGGCCCCCCTTTAAAAAA', '2*s*0*g*3*', 'init'],
        # ...
        # ['TTTGGGGGGGGGGGGGGGGGGAAAAAAAAAAAATTT', '421*s*0*g*5*', 'init'],
        # ['AAATTTTTTCCCCCCCCCGGGAAAAAAGGGAAATTT', '422*s*0*g*5*', 'init']]
        self.console.print('===>Sequencing library has been generated')

    def generate(self, ):
        """

        [['1' '1' '0']
         ['1' '1' '1']
         ['1' '1' '2']
         ...
         ['1' '1' '48']
         ['1' '1' '49']]
        Returns
        -------

        """
        ### +++++++++++++++ block: PCR amplification: Preparation +++++++++++++++
        self.console.print('===>PCR amplification starts...')
        self.console.print('======>Assign parameters...')
        # print(np.array(self.sequencing_library))
        for i, pcr_num_i in enumerate(self.pcr_nums):
            self.console.print('======>{}. PCR cycle number: {}'.format(i, pcr_num_i))
            pcr_ampl_params = {
                'read_lib_fpn': self.working_dir + 'sequencing_library.txt',

                'data': np.array(self.sequencing_library),
                'ampl_rate': self.ampl_rate,
                'pcr_error': self.pcr_error,
                'pcr_num': pcr_num_i,

                'err_route': self.err_route,
                'err_num_met': self.err_num_met,
                'use_seed': self.use_seed,
                'seed': self.seed,
                'recorder_nucleotide_num': [],
                'recorder_pcr_err_num': [],
                'recorder_pcr_read_num': [],

                'seq_sub_spl_number': self.seq_sub_spl_number,
                'seq_sub_spl_rate': self.seq_sub_spl_rate,

                'pcr_deletion': True,
                'pcr_insertion': True, # False True
                'pcr_del_rate': 2.4*10e-6,
                'pcr_ins_rate': 7.1*10e-7,

                'verbose': self.verbose,
            }
            # pcr_ampl_params['data']
            # [['GGGGGGCCCTTTTTTTTTTTTAAAAAACCCAAAGGG' '0*s*0*g*3*' 'init']
            #  ['GGGTTTGGGAAAAAAAAAAAATTTTTTGGGCCCGGG' '1*s*0*g*3*' 'init']
            #  ['CCCTTTCCCGGGGGGTTTGGGCCCCCCTTTAAAAAA' '2*s*0*g*3*' 'init']
            #  ...
            #  ['GGGAAAAAAAAAAAACCCAAACCCTTTAAATTTGGG' '420*s*0*g*5*' 'init']
            #  ['TTTGGGGGGGGGGGGGGGGGGAAAAAAAAAAAATTT' '421*s*0*g*5*' 'init']
            #  ['AAATTTTTTCCCCCCCCCGGGAAAAAAGGGAAATTT' '422*s*0*g*5*' 'init']]
            # print(pcr_ampl_params['data'][:, 1:3])
            # print(pcr_ampl_params)
            if pcr_ampl_params['err_route'] == 'tree':
                pcr_ampl_params['data'] = pcr_ampl_params['data'][:, 1:3]
            if pcr_ampl_params['err_route'] == 'minnow':
                pcr_ampl_params['data'] = pcr_ampl_params['data'][:, 1:3]
            if pcr_ampl_params['err_route'] == 'mutation_table_minimum' or pcr_ampl_params['err_route'] == 'mutation_table_complete':
                # print(pcr_ampl_params['data'][:, 0])
                def calc_len(a):
                    return len(a)
                vfunc = np.vectorize(calc_len)
                # [[36] vfunc(pcr_ampl_params['data'][:, 0])[:, np.newaxis]
                #  [36]
                #  [36]
                #  ...
                #  [36]
                #  [36]]
                pcr_ampl_params['data'] = np.hstack((
                    vfunc(pcr_ampl_params['data'][:, 0])[:, np.newaxis],
                    pcr_ampl_params['data'][:, 1:3],
                ))
                # pcr_ampl_params['data']
                # [['36' '0*s*0*g*3*' 'init']
                #  ['36' '1*s*0*g*3*' 'init']
                #  ['36' '2*s*0*g*3*' 'init']
                #  ...
                #  ['36' '421*s*0*g*5*' 'init']
                #  ['36' '422*s*0*g*5*' 'init']]
                col_0 = np.array([[1] for _ in range(pcr_ampl_params['data'].shape[0])])
                mut_info_table = np.hstack((col_0, col_0))
                col_2 = pcr_ampl_params['data'][:, 1].astype(str)[:, np.newaxis]
                mut_info_table = np.hstack((mut_info_table, col_2))
                pcr_ampl_params['mut_info'] = mut_info_table
                # print(mut_info_table)
                # pcr_ampl_params['mut_info'] = np.empty(shape=[0, 3])
                # print(pcr_ampl_params['mut_info'])

            ### +++++++++++++++ block: PCR amplification: simulation +++++++++++++++
            pcr_stime = time.time()
            pcr = self.pcr(pcr_params=pcr_ampl_params).np()
            # print(pcr.keys())
            self.console.print('======>PCR amplification completes in {}s'.format(time.time() - pcr_stime))

            ### +++++++++++++++ block: Subsampling: sequencing depth or rate +++++++++++++++
            if pcr_ampl_params['err_route'] == 'tree':
                pcr['data'] = self.subsampling.pcrtree(pcr_dict=pcr)

            # pcr['data']
            # [['68*s*0*g*3*_1_3_9' 'pcr-9']
            #  ['334*s*0*g*5*_1_2_4_5_6_7_8' 'pcr-8']
            #  ['54*s*0*g*5*_2_3_4_5_6' 'pcr-6']
            #  ...
            #  ['55*s*0*g*3*_2_3_4_5_6_7_9_10' 'pcr-10']
            #  ['175*s*0*g*5*_1_3_8_9_10' 'pcr-10']
            #  ['8*s*0*g*5*_2_6_7_10' 'pcr-10']]
            # pcr['data'].shape
            # 396201, 2
            if pcr_ampl_params['err_route'] == 'minnow':
                pcr['data'] = self.subsampling.minnow(pcr_dict=pcr)
            # pcr['data']
            # [['TTTTTTTTTCCCGGGGGGCCCGGGAAAGGGAAAGGG' '261*s*0*g*4*_1_3_5_7_8_9' 'pcr-9']
            #  ['AAAGGGTTTGGGCCCTTTAAAGGGGGGGGGAAAAAA' '175*s*0*g*4*_1_4_7' 'pcr-7']
            #  ['AAAGGGTTTGGGCCCTTTAAAGGGGGGGGGAAAAAA' '175*s*0*g*4*_1_4_5_8_9' 'pcr-9']
            # ...
            # ['GGGCCCAAACCCTTTGGGAAACCCGGGAAACCCGGG' '171*s*0*g*4*_2_3_4_6_8_9' 'pcr-9']
            #  ['GGGCCCCCCTTTCCCCCCTTTCCCTTTGGGAAAAAA' '226*s*0*g*4*_3_5_6_9' 'pcr-9']]
            # pcr['data'].shape
            # (200, 3)

            if pcr_ampl_params['err_route'] == 'mutation_table_minimum':
                pcr['data'] = self.subsampling.mutation_table_minimum(pcr_dict=pcr)

            if pcr_ampl_params['err_route'] == 'mutation_table_complete':
                pcr['data'] = self.subsampling.mutation_table_complete(pcr_dict=pcr)

            ### +++++++++++++++ block: Sequencing: parameters +++++++++++++++
            self.console.print('======>Sequencing starts')
            seq_params = {
                'data': pcr['data'],

                'seq_error': self.seq_error,
                'err_num_met': self.err_num_met,
                'use_seed': self.use_seed,
                'seed': self.seed,
                'verbose': self.verbose,

                'seq_sub_spl_number': self.seq_sub_spl_number,
                'seq_sub_spl_rate': self.seq_sub_spl_rate,

                'seq_deletion': True,
                'seq_insertion': True,  # False True
                'seq_del_rate': 2.4 * 10e-6,
                'seq_ins_rate': 7.1 * 10e-7,
            }
            seq = self.seq(seq_params=seq_params).np()
            self.console.print('=========>Sequencing has completed')
            self.console.print('=========>Reads write to files in FastQ format')
            self.wfastq().togz(
                list_2d=seq['data'],
                sv_fp=self.sv_fastq_fp,
                fn='pcr_num_' + str(i),
                # fn='pcr_num_' + "{:.2f}".format(pcr_num_i),
                symbol='-',
            )
            del seq
            self.console.print('=========>FastQ file is saved')
        self.console.print('======>Simulation completes')
        return


if __name__ == "__main__":
    from tresor.path import to

    from tresor.gsample.FromSimulator import fromSimulator

    gspl = fromSimulator(
        R_root='D:/Programming/R/R-4.3.2/',
        num_samples=2,
        num_genes=6,
        simulator='spsimseq',
    ).run()

    p = Gene(
        # initial sequence generation
        gspl=gspl,

        len_params={
            'umi': {
                'umi_unit_pattern': 3,
                'umi_unit_len': 12,
            },
            'umi_1': {
                'umi_unit_pattern': 3,
                'umi_unit_len': 12,
            },
            'barcode': 16,
            'seq': 100,
            'seq_2': 100,
            'adapter': 10,
            'adapter_1': 10,
            'primer': 10,
            'primer_1': 10,
            'spacer': 10,
            'spacer_1': 10,
        },
        seq_params={
            'custom': 'AAGC',
            'custom_1': 'A',
        },
        material_params={
            'fasta_cdna_fpn': to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),  # None False
        },
        seq_num=50,
        working_dir=to('data/simu/'),

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        condis=['umi'],
        # condis=['umi', 'seq'],
        # condis=['umi', 'custom', 'seq', 'custom_1'],
        sim_thres=3,
        permutation=0,

        # PCR amplification
        ampl_rate=0.9,
        err_route='minnow', # tree minnow err1d err2d mutation_table_minimum mutation_table_complete
        pcr_error=1e-04,
        pcr_nums=np.arange(1, 18 + 1, 1),
        err_num_met='nbinomial',
        seq_error=0.01,
        seq_sub_spl_number=200, # None
        # seq_sub_spl_rate=0.333,
        use_seed=True,
        seed=1,

        verbose=False, # True False

        mode='short_read',  # long_read short_read

        sv_fastq_fp=to('data/simu/'),
    )
    print(p.generate())