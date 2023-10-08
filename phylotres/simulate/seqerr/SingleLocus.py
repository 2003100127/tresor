__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import time
import numpy as np
from phylotres.library.SingleLocus import SingleLocus as simuip
from phylotres.pcr.Amplify import Amplify as pcr
from phylotres.seq.Calling import Calling as seq
from phylotres.pcr.Subsampling import Subsampling
from phylotres.util.sequence.fastq.Write import write as wfastq
from phylotres.util.Console import Console


class SingleLocus:

    def __init__(
            self,
            len_params,
            seq_num,
            seq_len,
            is_sv_umi_lib,
            is_sv_seq_lib,
            is_sv_primer_lib,
            is_sv_adapter_lib,
            is_sv_spacer_lib,
            fasta_cdna_fpn,
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

            seq_errors,
            seq_sub_spl_rate,

            sv_fastq_fp,

            verbose=True,
    ):
        self.seq_num = seq_num
        self.seq_len = seq_len
        self.is_sv_umi_lib = is_sv_umi_lib
        self.is_sv_seq_lib = is_sv_seq_lib
        self.working_dir = working_dir
        self.fasta_cdna_fpn = fasta_cdna_fpn
        self.condis = condis
        self.sim_thres = sim_thres
        self.permutation = permutation

        self.err_route = err_route
        self.ampl_rate = ampl_rate
        self.pcr_error = pcr_error
        self.pcr_num = pcr_num
        self.err_num_met = err_num_met
        self.use_seed = use_seed
        self.seed = seed

        self.seq_errors = seq_errors
        self.seq_sub_spl_rate = seq_sub_spl_rate
        self.sv_fastq_fp = sv_fastq_fp

        self.pcr = pcr
        self.seq = seq
        self.wfastq = wfastq
        self.subsampling = Subsampling()
        self.console = Console()
        self.console.verbose = verbose
        self.verbose = verbose

        # ### /*** block. Init a pool of sequences ***/
        self.console.print('===>Generating an init pool of sequences starts...')
        self.sequencing_library = simuip(
            len_params=len_params,
            fasta_cdna_fpn=fasta_cdna_fpn,
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
        ).pooling()
        self.console.print('===>Init pool of sequences has completed')

    def generate(self, ):
        # ### /*** block 1. PCR amplification: Preparation ***/
        self.console.print('===>PCR amplification starts...')
        self.console.print('======>Assign parameters...')
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
            'verbose': self.verbose,

        }
        # print(pcr_ampl_params['data'])
        if pcr_ampl_params['err_route'] == 'tree':
            pcr_ampl_params['data'] = pcr_ampl_params['data'][:, 1:3]
        if pcr_ampl_params['err_route'] == 'minnow':
            def calclen(a):
                return len(a)
            vfunc = np.vectorize(calclen)
            pcr_ampl_params['data'] = np.hstack((vfunc(pcr_ampl_params['data'][:, 0])[:, np.newaxis], pcr_ampl_params['data'][:, 1:3]))
            # print(pcr_ampl_params['data'])
            col_0 = np.array([[1] for _ in range(pcr_ampl_params['data'].shape[0])])
            cc = np.hstack((col_0, col_0))
            col_2 = pcr_ampl_params['data'][:, 1].astype(str)[:, np.newaxis]
            # print(col_2)
            cc = np.hstack((cc, col_2))
            # print(cc)
            pcr_ampl_params['mut_info'] = cc
            # pcr_ampl_params['mut_info'] = np.empty(shape=[0, 3])
            # print(pcr_ampl_params['mut_info'])

        # ### /*** block 2. PCR amplification: simulation ***/
        pcr_stime = time.time()
        pcr = self.pcr(pcr_params=pcr_ampl_params).np()
        # print(pcr.keys())
        self.console.print('======>PCR amplification completes in {}s'.format(time.time() - pcr_stime))

        # ### /*** block 3. Subsampling: sequencing depth or rate ***/
        # print(pcr['data'])
        # print(pcr['data'].shape)
        if pcr_ampl_params['err_route'] == 'tree':
            pcr['data'] = self.subsampling.pcrtree(pcr_dict=pcr)

        # print(pcr['data'])
        # print(pcr['data'].shape)

        if pcr_ampl_params['err_route'] == 'minnow':
            pcr['data'] = self.subsampling.minnow(pcr_dict=pcr)

        # ### /*** block 4. Sequencing: parameters ***/
        self.console.print('======>Sequencing has started...')
        for id, iseq_err in enumerate(self.seq_errors):
            seq_params = {
                'data': pcr['data'],
                'seq_sub_spl_rate': self.seq_sub_spl_rate,
                'seq_error': iseq_err,
                'err_num_met': self.err_num_met,
                'use_seed': self.use_seed,
                'seed': self.seed,
                'verbose': False,
            }
            seq = self.seq(seq_params=seq_params).np()
            self.console.print('=========>Sequencing has completed')
            self.console.print('=========>Write seqs in fastq format')
            self.wfastq().togz(
                list_2d=seq['data'],
                sv_fp=self.sv_fastq_fp,
                fn='seq_err_' + str(id),
                symbol='-',
            )
            del seq
            self.console.print('=========>Write seqs finished')
        self.console.print('======>Simulation ends...')
        return


if __name__ == "__main__":
    from phylotres.path import to

    p = SingleLocus(
        # initial sequence generation
        len_params={
            'umi': {
                'umi_unit_pattern': 3,
                'umi_unit_len': 12,
            },
            'seq': 100,
        },
        seq_num=50,
        seq_len=100,
        working_dir=to('data/simu/'),
        fasta_cdna_fpn=False,
        # fasta_cdna_fpn=to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        condis=['umi'],
        # condis=['umi', 'seq'],
        sim_thres=3,
        permutation=0,

        # PCR amplification
        ampl_rate=0.85,
        err_route='minnow', # tree minnow default err2d
        pcr_error=1e-4,
        pcr_num=14,
        err_num_met='nbinomial',
        seq_errors=[1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3],
        seq_sub_spl_rate=0.333,
        use_seed=True,
        seed=None,

        verbose=True,

        sv_fastq_fp=to('data/simu/'),
    )
    print(p.generate())