__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import time
import numpy as np
from tresor.library.SingleCell import SingleCell as scsimulib
from tresor.pcr.Amplify import Amplify as pcr
from tresor.pcr.Subsampling import Subsampling
from tresor.sequencing.Calling import Calling as seq
from tresor.util.sequence.fastq.Write import write as wfastq
from tresor.util.Console import Console


class SingleCell:

    def __init__(
            self,
            gmat,
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

            ampl_rate,
            err_route,
            pcr_error,
            pcr_num,
            err_num_met,
            use_seed,
            seed,

            sv_fastq_fp,

            seq_errors,
            seq_sub_spl_number=None,
            seq_sub_spl_rate=1/3,

            verbose=True,
            **kwargs,
    ):
        self.len_params = len_params
        self.seq_num = seq_num
        self.is_sv_umi_lib = is_sv_umi_lib
        self.is_sv_seq_lib = is_sv_seq_lib
        self.working_dir = working_dir
        self.condis = condis
        self.sim_thres = sim_thres
        self.permutation = permutation

        self.gmat = gmat

        self.err_route = err_route
        self.ampl_rate = ampl_rate
        self.pcr_error = pcr_error
        self.pcr_num = pcr_num
        self.err_num_met = err_num_met
        self.use_seed = use_seed
        self.seed = seed

        self.seq_errors = seq_errors
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
        self.sequencing_library = scsimulib(
            len_params=len_params,
            gmat=self.gmat,
            seq_num=seq_num,
            is_seed=use_seed,
            working_dir=working_dir,
            condis=condis,
            sim_thres=sim_thres,
            permutation=permutation,
            is_sv_umi_lib=is_sv_umi_lib,
            is_sv_seq_lib=is_sv_seq_lib,
            is_sv_primer_lib=is_sv_primer_lib,
            is_sv_adapter_lib=is_sv_adapter_lib,
            is_sv_spacer_lib=is_sv_spacer_lib,
            verbose=self.verbose,
            mode=self.kwargs['mode'],
            material_params=self.kwargs['material_params'] if 'material_params' in self.kwargs.keys() else None,
            seq_params=self.kwargs['seq_params'] if 'seq_params' in self.kwargs.keys() else None,
        ).pooling()
        # print(self.sequencing_library)
        # print(len(self.sequencing_library))
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
        pcr_ampl_params = {
            'read_lib_fpn': self.working_dir + 'sequencing_library.txt',

            'data': np.array(self.sequencing_library),
            'ampl_rate': self.ampl_rate,
            'pcr_error': self.pcr_error,
            'pcr_num': self.pcr_num,

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
        # single-cell
        # [['GAAATCATGTAGTTCGGGGGGGCCCTTTTTTTTTTTTAAAAAACCCAAAGGG' '0*c*0*g*1*' 'init']
        #  ['GAAATCATGTAGTTCGCCCCCCAAATTTAAACCCAAATTTCCCAAAAAAAAA' '0*c*0*g*4*' 'init']
        #  ['GAAATCATGTAGTTCGAAACCCGGGGGGCCCTTTCCCTTTGGGAAATTTCCC' '1*c*0*g*4*' 'init']
        # ...
        # ['CGCGTTAGTAATTCATAAAGGGGGGCCCAAACCCGGGGGGGGGGGGTTTCCC' '221*c*1*g*4*' 'init']
        #  ['CGCGTTAGTAATTCATAAATTTCCCGGGCCCCCCGGGAAAGGGCCCCCCTTT' '222*c*1*g*4*' 'init']]
        if pcr_ampl_params['err_route'] == 'bftree':
            pcr_ampl_params['data'] = pcr_ampl_params['data'][:, 1:3]
        if pcr_ampl_params['err_route'] == 'sptree':
            pcr_ampl_params['data'] = pcr_ampl_params['data'][:, 1:3]
        if pcr_ampl_params['err_route'] == 'mutation_table_minimum' or pcr_ampl_params['err_route'] == 'mutation_table_complete':
            # print(pcr_ampl_params['data'][:, 0])
            def calc_len(a):
                return len(a)
            vfunc = np.vectorize(calc_len)
            print(vfunc(pcr_ampl_params['data'][:, 0])[:, np.newaxis])
            # [[52] vfunc(pcr_ampl_params['data'][:, 0])[:, np.newaxis]
            #  [52]
            #  [52]
            #  ...
            #  [52]
            #  [52]]
            pcr_ampl_params['data'] = np.hstack((
                vfunc(pcr_ampl_params['data'][:, 0])[:, np.newaxis],
                pcr_ampl_params['data'][:, 1:3],
            ))
            # pcr_ampl_params['data']
            # [['52' '0*c*0*g*1*' 'init']
            #  ['52' '0*c*0*g*4*' 'init']
            #  ['52' '1*c*0*g*4*' 'init']
            #  ...
            #   ['52' '221*c*1*g*4*' 'init']
            #  ['52' '222*c*1*g*4*' 'init']]
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
        # print(pcr['data'])
        # print(pcr['data'].shape)
        if pcr_ampl_params['err_route'] == 'bftree':
            pcr['data'] = self.subsampling.bftree(pcr_dict=pcr)
        # print(pcr['data'])
        # print(pcr['data'].shape)

        if pcr_ampl_params['err_route'] == 'sptree':
            pcr['data'] = self.subsampling.sptree(pcr_dict=pcr)

        if pcr_ampl_params['err_route'] == 'mutation_table_minimum':
            pcr['data'] = self.subsampling.mutation_table_minimum(pcr_dict=pcr)

        if pcr_ampl_params['err_route'] == 'mutation_table_complete':
            pcr['data'] = self.subsampling.mutation_table_complete(pcr_dict=pcr)

        ### +++++++++++++++ block: Sequencing: parameters +++++++++++++++
        self.console.print('======>Sequencing starts')
        for i, seq_err_i in enumerate(self.seq_errors):
            self.console.print('======>{}. Sequencing error rate: {}'.format(i, seq_err_i))
            seq_params = {
                'data': pcr['data'],

                'seq_error': seq_err_i,
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
                fn='seq_err_' + str(i),
                symbol='-',
            )
            del seq
            self.console.print('=========>FastQ file is saved')
        self.console.print('======>Simulation completes')
        return


if __name__ == "__main__":
    from tresor.path import to

    from tresor.gcell.FromSimulator import fromSimulator

    gmat, _, _ = fromSimulator(
        simulator='spsimseq',
        R_root='D:/Programming/R/R-4.3.2/',
        num_genes=6,
        num_cells=2,
    ).run()

    p = SingleCell(
        # initial sequence generation
        gmat=gmat,

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
        # condis=['umi'],
        # condis=['umi', 'seq'],
        condis=['barcode', 'umi'],
        # condis=['barcode', 'custom', 'umi', 'custom_1'],
        sim_thres=3,
        permutation=0,

        # PCR amplification
        ampl_rate=0.85,
        err_route='mutation_table_complete', # bftree sptree err1d err2d mutation_table_minimum mutation_table_complete
        pcr_error=1e-4,
        pcr_num=10,
        err_num_met='nbinomial',
        seq_errors=[1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3],
        seq_sub_spl_number=200, # None
        seq_sub_spl_rate=0.333,
        use_seed=True,
        seed=1,

        verbose=False, # True False

        mode='short_read',  # long_read short_read

        sv_fastq_fp=to('data/simu/'),
    )
    print(p.generate())