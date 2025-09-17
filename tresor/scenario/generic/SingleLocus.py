__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


import time
import numpy as np
import pandas as pd
from tresor.library.SingleLocus import SingleLocus as simuip
from tresor.pcr.Amplify import Amplify as pcr
from tresor.sequencing.Calling import Calling as seq
from tresor.pcr.Subsampling import Subsampling
from tresor.util.sequence.fastq.Write import write as wfastq
from tresor.util.file.Writer import Writer as fwriter
from tresor.util.Console import Console


class SingleLocus:

    def __init__(
            self,
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

            seq_error,
            err_route,
            pcr_error,
            pcr_num,
            err_num_met,
            use_seed,
            seed,

            ampl_rate,
            sv_fastq_fp,
            sv_fastq_fn,

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

        # self.bead_mutation = kwargs['bead_mutation']
        # self.bead_mut_rate = kwargs['bead_mut_rate']
        # self.bead_deletion = kwargs['bead_deletion']
        # self.bead_del_rate = kwargs['bead_del_rate']
        # self.bead_insertion = kwargs['bead_insertion']
        # self.bead_ins_rate = kwargs['bead_ins_rate']

        self.err_route = err_route
        self.ampl_rate = ampl_rate
        self.pcr_num = pcr_num
        self.pcr_error = pcr_error
        self.seq_error = seq_error
        self.err_num_met = err_num_met
        self.use_seed = use_seed
        self.seed = seed

        self.seq_sub_spl_number = seq_sub_spl_number
        self.seq_sub_spl_rate = seq_sub_spl_rate
        self.sv_fastq_fp = sv_fastq_fp
        self.sv_fastq_fn = sv_fastq_fn

        self.pcr = pcr
        self.seq = seq
        self.wfastq = wfastq
        self.subsampling = Subsampling(verbose=verbose)
        self.fwriter = fwriter()

        self.kwargs = kwargs
        print(self.kwargs)

        self.console = Console()
        self.console.verbose = verbose
        self.verbose = verbose

        ### +++++++++++++++ block: generate sequencing library +++++++++++++++
        self.console.print('===>Sequencing library generation starts')
        self.sequencing_library, self.lib_err_mark = simuip(
            len_params=self.len_params,
            seq_num=self.seq_num,
            is_seed=self.use_seed,
            working_dir=self.working_dir,
            condis=self.condis,
            sim_thres=self.sim_thres,
            # bead_mutation=self.bead_mutation,
            # bead_mut_rate=self.bead_mut_rate,
            # bead_deletion=self.bead_deletion,
            # bead_del_rate=self.bead_del_rate,
            # bead_insertion=self.bead_insertion,
            bead_mutation=self.kwargs['bead_mutation'] if 'bead_mutation' in self.kwargs.keys() else False,
            bead_mut_rate=self.kwargs['bead_mut_rate'] if 'bead_mut_rate' in self.kwargs.keys() else False,
            bead_deletion=self.kwargs['bead_deletion'] if 'bead_deletion' in self.kwargs.keys() else False,
            bead_del_rate=self.kwargs['bead_del_rate'] if 'bead_del_rate' in self.kwargs.keys() else False,
            bead_insertion=self.kwargs['bead_insertion'] if 'bead_insertion' in self.kwargs.keys() else False,
            bead_ins_rate=self.kwargs['bead_del_rate'] if 'bead_del_rate' in self.kwargs.keys() else False,
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
        self.lib_err_mark['pcr_err_mark'] = False
        # print(self.lib_err_mark)
        self.sequencing_library = pd.concat([
            pd.DataFrame(self.sequencing_library),
            self.lib_err_mark,
        ], axis=1).values
        # print(self.sequencing_library)
        self.console.print('===>Sequencing library has been generated')
        # print(self.lib_err_mark)

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
        time_arr = []
        satime = time.time()

        self.console.print('===>PCR amplification starts...')
        self.console.print('======>Assign parameters...')
        # print(np.array(self.sequencing_library))

        self.console.print('======>The criterion: {}'.format('12'))
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

            # 'pcr_deletion': False,
            # 'pcr_insertion': False, # False True
            # 'pcr_del_rate': 2.4*10e-6,
            # 'pcr_ins_rate': 7.1*10e-7,

            'pcr_deletion': self.kwargs['pcr_deletion'] if 'pcr_deletion' in self.kwargs.keys() else False,
            'pcr_insertion': self.kwargs['pcr_insertion'] if 'pcr_insertion' in self.kwargs.keys() else False,
            'pcr_del_rate': self.kwargs['pcr_del_rate'] if 'pcr_del_rate' in self.kwargs.keys() else 0,
            'pcr_ins_rate': self.kwargs['pcr_ins_rate'] if 'pcr_ins_rate' in self.kwargs.keys() else 0,

            'verbose': self.verbose,
        }
        # print(pcr_ampl_params['data'])
        # if checking bead synthesis errors, you will get the pcr_ampl_params['data'] like this
        # [['GAT...A' '0' 'init' False False False False]
        #  ['GTA...G' '1' 'init' False False False False]
        #  ['ATG...T' '2' 'init' False False True False]
        #  ['CAA...T' '3' 'init' False False False False]
        # ...
        #  ['TTA...C' '49' 'init' False False False False]]
        # in the normal state, pcr_ampl_params['data'] looks like:
        # pcr_ampl_params['data']
        # [['GGGAAATTTAAACCCTTTAAAGGGAAAAAAGGGCCC' '0' 'init']
        #  ['GGGTTTAAACCCCCCCCCGGGAAATTTTTTGGGTTT' '1' 'init']
        #  ['AAATTTGGGCCCGGGAAAGGGCCCAAAAAAGGGAAA' '2' 'init']
        #  ...
        #  ['CCCAAAAAAGGGCCCAAAGGGCCCTTTGGGGGGCCC' '48' 'init']
        #  ['TTTTTTAAATTTAAAAAAGGGAAAGGGGGGGGGCCC' '49' 'init']]
        # print(pcr_ampl_params['data'][:, 1:3])
        # print(pcr_ampl_params)
        if pcr_ampl_params['err_route'] == 'err1d':
            pcr_ampl_params['data'] = pcr_ampl_params['data'][:, 1:3]
        if pcr_ampl_params['err_route'] == 'err2d':
            pcr_ampl_params['data'] = pcr_ampl_params['data'][:, 1:3]
        if pcr_ampl_params['err_route'] == 'bftree':
            pcr_ampl_params['data'] = pcr_ampl_params['data'][:, 1:3]
        if pcr_ampl_params['err_route'] == 'sptree':
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
            # print(pcr_ampl_params)
            # print(pcr_ampl_params['data'])
            # pcr_ampl_params['data']
            # [['36' '0' 'init']
            #  ['36' '1' 'init']
            #  ['36' '2' 'init']
            #  ...
            #  ['36' '48' 'init']
            #  ['36' '49' 'init']]
            col_0 = np.array([[1] for _ in range(pcr_ampl_params['data'].shape[0])])
            mut_info_table = np.hstack((col_0, col_0))
            col_2 = pcr_ampl_params['data'][:, 1].astype(str)[:, np.newaxis]
            mut_info_table = np.hstack((mut_info_table, col_2))
            pcr_ampl_params['mut_info'] = mut_info_table
            # print(mut_info_table)
            # print(self.lib_err_mark.values)
            # pcr_ampl_params['mut_info'] = np.empty(shape=[0, 3])
            # print(pcr_ampl_params['mut_info'])
        if pcr_ampl_params['err_route'] == 'mutation_table_complete':
            def calc_len(a):
                return len(a)
            vfunc = np.vectorize(calc_len)
            pcr_ampl_params['data'] = np.hstack((
                vfunc(pcr_ampl_params['data'][:, 0])[:, np.newaxis],
                pcr_ampl_params['data'][:, 1:7],
            ))
            col_0 = np.array([[1] for _ in range(pcr_ampl_params['data'].shape[0])])
            mut_info_table = np.hstack((col_0, col_0))
            col_2 = pcr_ampl_params['data'][:, 1].astype(str)[:, np.newaxis]
            mut_info_table = np.hstack((mut_info_table, col_2))
            mut_info_table = pd.concat([pd.DataFrame(mut_info_table), self.lib_err_mark], axis=1).values
            pcr_ampl_params['mut_info'] = mut_info_table

        ### +++++++++++++++ block: PCR amplification: simulation +++++++++++++++
        pcr_stime = time.time()
        pcr = self.pcr(pcr_params=pcr_ampl_params).np()
        # print(pcr.keys())
        self.console.print('======>PCR amplification completes in {}s'.format(time.time() - pcr_stime))
        # print(pcr['data'])
        # print(pcr['data'].shape)

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
            # print(pcr['data'])
            # print(pcr['data'].shape)

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

            # 'seq_deletion': False,
            # 'seq_insertion': False, # False True
            # 'seq_del_rate': 2.4 * 10e-6,
            # 'seq_ins_rate': 7.1 * 10e-7,

            'seq_deletion': self.kwargs['seq_deletion'] if 'seq_deletion' in self.kwargs.keys() else False,
            'seq_insertion': self.kwargs['seq_insertion'] if 'seq_insertion' in self.kwargs.keys() else False,
            'seq_del_rate': self.kwargs['seq_del_rate'] if 'seq_del_rate' in self.kwargs.keys() else 0,
            'seq_ins_rate': self.kwargs['seq_ins_rate'] if 'seq_ins_rate' in self.kwargs.keys() else 0,
        }
        seq = self.seq(seq_params=seq_params).np()
        self.console.print('=========>Sequencing has completed')
        print('======>simulation completes in {}s'.format(time.time() - satime))
        time_arr.append(time.time() - satime)

        self.console.print('=========>Reads write to files in FastQ format')
        self.wfastq().togz(
            list_2d=seq['data'],
            sv_fp=self.sv_fastq_fp,
            fn=self.sv_fastq_fn,
            symbol='-',
        )
        del seq

        self.console.print('=========>FastQ file is saved')
        self.console.print('======>Simulation completes')
        return {
            'time_arr': time_arr,
        }


if __name__ == "__main__":
    from tresor.path import to

    p = SingleLocus(
        # initial sequence generation
        len_params={
            'umi': {
                'umi_unit_pattern': 1,
                'umi_unit_len': 12,
            },
            'seq': 100,
        },
        seq_params={
            'custom': 'AAGC',
            'custom_1': 'A',
        },
        material_params={
            # 'fasta_cdna_fpn': to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),  # None False
        },
        seq_num=50,
        working_dir=to('data/simu/'),

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        # condis=['umi'],
        condis=['umi', 'seq'],
        # condis=['umi', 'custom', 'seq', 'custom_1'],
        sim_thres=3,
        permutation=0,

        # PCR amplification
        ampl_rate=0.85,
        err_route='mutation_table_complete', # bftree sptree err1d err2d mutation_table_minimum mutation_table_complete
        pcr_error=1e-4,
        pcr_num=10,
        err_num_met='nbinodmial',
        seq_error=0.01,
        seq_sub_spl_number=500, # None 200
        seq_sub_spl_rate=False, # 0.333
        use_seed=True,
        seed=1,

        bead_mutation=True,  # True False
        bead_mut_rate=1e-4,  # 0.016 0.00004
        bead_deletion=True,  # True False
        bead_insertion=True,
        bead_del_rate=0.1/112,  # 0.016 0.00004, 2.4e-7
        bead_ins_rate=7.1e-7,  # 0.011 0.00001, 7.1e-7

        pcr_deletion=True,  # True False
        pcr_insertion=True,
        pcr_del_rate=2.4e-6,  # 0.016 0.00004
        pcr_ins_rate=7.1e-7,  # 0.011 0.00001
        seq_deletion=False,
        seq_insertion=False,
        seq_del_rate=2.4e-6,
        seq_ins_rate=7.1e-7,

        verbose=True, # True False

        mode='short_read',  # long_read short_read

        sv_fastq_fp=to('data/simu/'),
        sv_fastq_fn='example',
    )
    print(p.generate())