__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


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
            pcr_errors,
            seq_error,
            pcr_num,
            err_num_met,
            use_seed,
            seed,

            ampl_rate,
            sv_fastq_fp,

            vis_num_err_per_read=False,

            seq_sub_spl_number=None,
            seq_sub_spl_rate=1/3,

            verbose=True,
            **kwargs,
    ):
        print(kwargs)
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
        self.vis_num_err_per_read = vis_num_err_per_read

        self.gspl = gspl

        self.err_route = err_route
        self.ampl_rate = ampl_rate
        self.seq_error = seq_error
        self.pcr_num = pcr_num
        self.err_num_met = err_num_met
        self.use_seed = use_seed
        self.seed = seed

        self.pcr_errors = pcr_errors
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
        satime = time.time()
        for i, pcr_error_i in enumerate(self.pcr_errors):
            self.console.print('======>{}. PCR error rate: {}'.format(i, pcr_error_i))
            pcr_ampl_params = {
                'read_lib_fpn': self.working_dir + 'sequencing_library.txt',

                'data': np.array(self.sequencing_library),
                'ampl_rate': self.ampl_rate,
                'pcr_error': pcr_error_i,
                'pcr_num': self.pcr_num,

                'err_route': self.err_route,
                'err_num_met': self.err_num_met,
                'use_seed': self.use_seed,
                'seed': self.seed,
                'recorder_nucleotide_num': [],
                'recorder_pcr_err_num': [],
                'recorder_pcr_read_num': [],

                'num_err_per_read_dict': {},

                'seq_sub_spl_number': self.seq_sub_spl_number,
                'seq_sub_spl_rate': self.seq_sub_spl_rate,

                'vis_num_err_per_read': self.vis_num_err_per_read,

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
            print(pcr)

            from tresor.vis.Binomial import plot_binomial_panels
            from tresor.vis.NBinomial import plot_nbinom_panels
            from tresor.vis.BetaBinomial import plot_betabinom_panels
            from tresor.vis.Poisson import plot_poisson_panels
            from tresor.vis.PlotError import analyze_tresor_pcr_counts

            if pcr['num_err_per_read_dict']:
                if pcr['err_num_met'] == 'binomial':
                    plot_binomial_panels(
                        pcr['num_err_per_read_dict'],
                        n_trials=None,
                        cols=4,
                        read_len=80, # 10x r1
                        suptitle="Binomial",
                    )
                    # plot_from_error_series(
                    #     counts_series=pcr['num_err_per_read_dict'],
                    #     readlen=100,  # 若要与 polyester(100bp) 对比；不清楚就留 None
                    #     p_theory=None,  # 有理论 p 就填；否则自动用 mean/n 估
                    #     title="PCR errors"
                    # )
                    analyze_tresor_pcr_counts(
                        num_err_per_read_dict=pcr['num_err_per_read_dict'],
                        readlen=100,
                        p_theory=None,
                        title="PCR"
                    )
                elif pcr['err_num_met'] == 'nbinomial':
                    plot_nbinom_panels(
                        pcr['num_err_per_read_dict'],
                        cols=4,
                        read_len=80,
                        show_poisson=True, # 10x r1
                        suptitle="Negative binomial",
                    )
                    # plot_from_error_series(
                    #     counts_series=pcr['num_err_per_read_dict'],
                    #     readlen=100,  # 若要与 polyester(100bp) 对比；不清楚就留 None
                    #     p_theory=None,  # 有理论 p 就填；否则自动用 mean/n 估
                    #     title="PCR errors"
                    # )
                    analyze_tresor_pcr_counts(
                        num_err_per_read_dict=pcr['num_err_per_read_dict'],
                        readlen=100,
                        p_theory=None,
                        title="PCR"
                    )
                elif pcr['err_num_met'] == 'betabinomial':
                    plot_betabinom_panels(
                        pcr['num_err_per_read_dict'],
                        # n_trials=8,
                        cols=4,
                        read_len=80,
                        show_binomial=False,
                        suptitle="Beta-binomial",
                    )
                    analyze_tresor_pcr_counts(
                        num_err_per_read_dict=pcr['num_err_per_read_dict'],
                        readlen=100,
                        p_theory=None,
                        title="PCR"
                    )
                elif pcr['err_num_met'] == 'poisson':
                    plot_poisson_panels(
                        pcr['num_err_per_read_dict'],
                        cols=4,
                        read_len=80,
                        suptitle="poisson",
                    )
                    analyze_tresor_pcr_counts(
                        num_err_per_read_dict=pcr['num_err_per_read_dict'],
                        readlen=100,
                        p_theory=None,
                        title="PCR"
                    )
                import matplotlib.pyplot as plt
                plt.show()

            self.console.print('======>PCR amplification completes in {}s'.format(time.time() - pcr_stime))

            ### +++++++++++++++ block: Subsampling: sequencing depth or rate +++++++++++++++
            if pcr_ampl_params['err_route'] == 'bftree':
                pcr['data'] = self.subsampling.bftree(pcr_dict=pcr)

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
            if pcr_ampl_params['err_route'] == 'sptree':
                pcr['data'] = self.subsampling.sptree(pcr_dict=pcr)
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

                'num_err_per_read_dict': {},

                'seq_sub_spl_number': self.seq_sub_spl_number,
                'seq_sub_spl_rate': self.seq_sub_spl_rate,

                'vis_num_err_per_read': self.vis_num_err_per_read,

                'seq_deletion': True,
                'seq_insertion': True,  # False True
                'seq_del_rate': 2.4 * 10e-6,
                'seq_ins_rate': 7.1 * 10e-7,
            }
            seq = self.seq(seq_params=seq_params).np()
            self.console.print('=========>Sequencing has completed')
            self.console.print('=========>Reads write to files in FastQ format')
            print('======>simulation completes in {}s'.format(time.time() - satime))

            if seq['num_err_per_read_dict']:
                if seq['err_num_met'] == 'binomial':
                    plot_binomial_panels(
                        seq['num_err_per_read_dict'],
                        n_trials=None,
                        cols=1,
                        read_len=80, # 10x r1
                        fs_sub_row=5.2,
                        fs_sub_col=4.2,
                        suptitle="Binomial",
                    )
                elif seq['err_num_met'] == 'nbinomial':
                    plot_nbinom_panels(
                        seq['num_err_per_read_dict'],
                        cols=1,
                        read_len=80, # 10x r1
                        fs_sub_row=5.2,
                        fs_sub_col=4.2,
                        show_poisson=False,
                        suptitle="Negative binomial",
                    )
                import matplotlib.pyplot as plt
                plt.show()

            self.wfastq().togz(
                list_2d=seq['data'],
                sv_fp=self.sv_fastq_fp,
                fn='pcr_err_' + str(i),
                # fn='pcr_err_' + "{:.2f}".format(pcr_err_i),
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
        err_route='sptree', # bftree sptree err1d err2d mutation_table_minimum mutation_table_complete
        pcr_errors = [1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01],
        pcr_num=10,
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