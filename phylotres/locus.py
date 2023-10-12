__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from phylotres.scenario.seqerr.SingleLocus import SingleLocus as seqerr
from phylotres.scenario.pcrerr.SingleLocus import SingleLocus as pcrerr
from phylotres.scenario.pcrnum.SingleLocus import SingleLocus as pcrnum
from phylotres.scenario.amplrate.SingleLocus import SingleLocus as amplrate
from phylotres.scenario.umilen.SingleLocus import SingleLocus as umilen
from phylotres.scenario.seqdep.SingleLocus import SingleLocus as seqdep


def simu_seq_err(
        len_params,
        seq_num,
        seq_len,
        working_dir,

        condis,
        sim_thres,
        permutation,

        ampl_rate,
        err_route,
        pcr_error,
        pcr_num,
        err_num_met,
        seq_errors,
        seq_sub_spl_number,
        seq_sub_spl_rate,

        use_seed=True,
        seed=1,

        fasta_cdna_fpn=None,

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        sv_fastq_fp=True,

        verbose=True,

        **kwargs,
):
    seqerr(
        # initial sequence generation
        len_params=len_params,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        seq_len=seq_len,
        working_dir=working_dir,
        fasta_cdna_fpn=fasta_cdna_fpn,

        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,

        # PCR amplification
        ampl_rate=ampl_rate,
        err_route=err_route,
        pcr_error=pcr_error,
        pcr_num=pcr_num,
        err_num_met=err_num_met,
        seq_errors=seq_errors,
        seq_sub_spl_number=seq_sub_spl_number,
        seq_sub_spl_rate=seq_sub_spl_rate,
        use_seed=use_seed,
        seed=seed,

        verbose=verbose,

        sv_fastq_fp=sv_fastq_fp,
    ).generate()
    return


def simu_pcr_err(
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
        err_route,
        pcr_errors,
        seq_error,
        pcr_num,
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
    pcrerr(
        len_params=len_params,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        seq_len=seq_len,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        fasta_cdna_fpn=fasta_cdna_fpn,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        pcr_errors=pcr_errors,
        seq_error=seq_error,
        pcr_num=pcr_num,
        err_num_met=err_num_met,
        use_seed=use_seed,
        seed=seed,
        ampl_rate=ampl_rate,
        sv_fastq_fp=sv_fastq_fp,
        seq_sub_spl_number=seq_sub_spl_number,
        seq_sub_spl_rate=seq_sub_spl_rate,
        verbose=verbose,
    ).generate()
    return 'Finished!'


def simu_pcr_num(
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
    pcrnum(
        len_params=len_params,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        seq_len=seq_len,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        fasta_cdna_fpn=fasta_cdna_fpn,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        pcr_error=pcr_error,
        seq_error=seq_error,
        pcr_nums=pcr_nums,
        err_num_met=err_num_met,
        use_seed=use_seed,
        seed=seed,
        ampl_rate=ampl_rate,
        sv_fastq_fp=sv_fastq_fp,
        seq_sub_spl_number=seq_sub_spl_number,
        seq_sub_spl_rate=seq_sub_spl_rate,
        verbose=verbose,
    ).generate()
    return 'Finished!'


def simu_ampl_rate(
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
        err_route,
        pcr_error,
        seq_error,
        pcr_num,
        err_num_met,
        use_seed,
        seed,
        ampl_rates,
        sv_fastq_fp,
        seq_sub_spl_number=None,
        seq_sub_spl_rate=1/3,
        verbose=True,
        **kwargs,
):
    amplrate(
        len_params=len_params,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        seq_len=seq_len,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        fasta_cdna_fpn=fasta_cdna_fpn,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        pcr_error=pcr_error,
        seq_error=seq_error,
        pcr_num=pcr_num,
        err_num_met=err_num_met,
        use_seed=use_seed,
        seed=seed,
        ampl_rates=ampl_rates,
        sv_fastq_fp=sv_fastq_fp,
        seq_sub_spl_number=seq_sub_spl_number,
        seq_sub_spl_rate=seq_sub_spl_rate,
        verbose=verbose,
    ).generate()
    return 'Finished!'


def simu_umi_len(
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
        err_route,
        ampl_rate,
        pcr_error,
        seq_error,
        pcr_num,
        err_num_met,
        use_seed,
        seed,
        sv_fastq_fp,
        seq_sub_spl_number=None,
        seq_sub_spl_rate=1/3,
        verbose=True,
        **kwargs,
):
    umilen(
        len_params=len_params,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        seq_len=seq_len,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        fasta_cdna_fpn=fasta_cdna_fpn,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        ampl_rate=ampl_rate,
        pcr_error=pcr_error,
        seq_error=seq_error,
        pcr_num=pcr_num,
        err_num_met=err_num_met,
        use_seed=use_seed,
        seed=seed,
        sv_fastq_fp=sv_fastq_fp,
        seq_sub_spl_number=seq_sub_spl_number,
        seq_sub_spl_rate=seq_sub_spl_rate,
        verbose=verbose,
    ).generate()
    return 'Finished!'


def simu_seq_dep(
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
        err_route,
        pcr_error,
        seq_error,
        pcr_num,
        err_num_met,
        use_seed,
        seed,
        ampl_rate,
        sv_fastq_fp,
        seq_sub_spl_numbers=None,
        seq_sub_spl_rate=1/3,
        verbose=True,
        **kwargs,
):
    seqdep(
        len_params=len_params,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        seq_len=seq_len,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        fasta_cdna_fpn=fasta_cdna_fpn,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        pcr_error=pcr_error,
        seq_error=seq_error,
        pcr_num=pcr_num,
        err_num_met=err_num_met,
        use_seed=use_seed,
        seed=seed,
        ampl_rate=ampl_rate,
        sv_fastq_fp=sv_fastq_fp,
        seq_sub_spl_numbers=seq_sub_spl_numbers,
        seq_sub_spl_rate=seq_sub_spl_rate,
        verbose=verbose,
    ).generate()
    return 'Finished!'


if __name__ == "__main__":
    from phylotres.path import to
    import numpy as np

    # simu_seq_err(
    #     len_params={
    #         'umi': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
    #         'seq': 100,
    #     },
    #     seq_params={
    #         'custom': 'AAGC',
    #         'custom_1': 'A',
    #     },
    #     seq_num=50,
    #     seq_len=100,
    #     working_dir=to('data/simu/'),
    #     fasta_cdna_fpn=False,
    #     # fasta_cdna_fpn=to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),
    #
    #     is_sv_umi_lib=True,
    #     is_sv_seq_lib=True,
    #     is_sv_primer_lib=True,
    #     is_sv_adapter_lib=True,
    #     is_sv_spacer_lib=True,
    #     # condis=['umi'],
    #     # condis=['umi', 'seq'],
    #     condis=['umi', 'custom', 'seq', 'custom_1'],
    #     sim_thres=3,
    #     permutation=0,
    #
    #     # PCR amplification
    #     ampl_rate=0.85,
    #     err_route='minnow',  # tree minnow err1d err2d mutation_table_minimum mutation_table_complete
    #     pcr_error=1e-4,
    #     pcr_num=10,
    #     err_num_met='nbinomial',
    #     seq_errors=[
    #         1e-05,
    #         2.5e-05,
    #         5e-05,
    #         7.5e-05,
    #         0.0001,
    #         0.00025,
    #         0.0005,
    #         0.00075,
    #         0.001,
    #         0.0025,
    #         0.005,
    #         0.0075,
    #         0.01,
    #         0.025,
    #         0.05,
    #         0.075,
    #         0.1,
    #         0.2,
    #         0.3,
    #     ],
    #     seq_sub_spl_number=200,  # None
    #     seq_sub_spl_rate=0.333,
    #     use_seed=True,
    #     seed=1,
    #
    #     verbose=False,  # True
    #
    #     sv_fastq_fp=to('data/simu/'),
    # )

    # simu_pcr_err(
    #     # initial sequence generation
    #     len_params={
    #         'umi': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
    #         'umi_1': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
    #         'barcode': 16,
    #         'seq': 100,
    #         'seq_2': 100,
    #         'adapter': 10,
    #         'adapter_1': 10,
    #         'primer': 10,
    #         'primer_1': 10,
    #         'spacer': 10,
    #         'spacer_1': 10,
    #     },
    #     seq_params={
    #         'custom': 'AAGC',
    #         'custom_1': 'A',
    #     },
    #     seq_num=50,
    #     seq_len=100,
    #     working_dir=to('data/simu/'),
    #     fasta_cdna_fpn=False,
    #     # fasta_cdna_fpn=to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),
    #
    #     is_sv_umi_lib=True,
    #     is_sv_seq_lib=True,
    #     is_sv_primer_lib=True,
    #     is_sv_adapter_lib=True,
    #     is_sv_spacer_lib=True,
    #     condis=['umi'],
    #     # condis=['umi', 'seq'],
    #     # condis=['umi', 'custom', 'seq', 'custom_1'],
    #     sim_thres=3,
    #     permutation=0,
    #
    #     # PCR amplification
    #     ampl_rate=0.9,
    #     err_route='minnow',  # tree minnow err1d err2d mutation_table_minimum mutation_table_complete
    #     pcr_errors=[1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075,
    #                 0.01],
    #     pcr_num=10,
    #     err_num_met='nbinomial',
    #     seq_error=0.01,
    #     seq_sub_spl_number=200,  # None
    #     # seq_sub_spl_rate=0.333,
    #     use_seed=True,
    #     seed=1,
    #
    #     verbose=False,  # True False
    #
    #     sv_fastq_fp=to('data/simu/'),
    # )

    # simu_pcr_num(
    #     len_params={
    #         'umi': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
    #         'umi_1': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
    #         'barcode': 16,
    #         'seq': 100,
    #         'seq_2': 100,
    #         'adapter': 10,
    #         'adapter_1': 10,
    #         'primer': 10,
    #         'primer_1': 10,
    #         'spacer': 10,
    #         'spacer_1': 10,
    #     },
    #     seq_params={
    #         'custom': 'AAGC',
    #         'custom_1': 'A',
    #     },
    #     seq_num=50,
    #     seq_len=100,
    #     working_dir=to('data/simu/'),
    #     fasta_cdna_fpn=False,
    #     # fasta_cdna_fpn=to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),
    #     is_sv_umi_lib=True,
    #     is_sv_seq_lib=True,
    #     is_sv_primer_lib=True,
    #     is_sv_adapter_lib=True,
    #     is_sv_spacer_lib=True,
    #     condis=['umi'],
    #     # condis=['umi', 'seq'],
    #     # condis=['umi', 'custom', 'seq', 'custom_1'],
    #     sim_thres=3,
    #     permutation=0,
    #
    #     # PCR amplification
    #     ampl_rate=0.9,
    #     err_route='minnow',  # tree minnow err1d err2d mutation_table_minimum mutation_table_complete
    #     pcr_error=1e-04,
    #     pcr_nums=np.arange(1, 18 + 1, 1),
    #     err_num_met='nbinomial',
    #     seq_error=0.01,
    #     seq_sub_spl_number=90,  # None 200
    #     # seq_sub_spl_rate=0.333,
    #     use_seed=True,
    #     seed=1,
    #
    #     verbose=False,  # True False
    #
    #     sv_fastq_fp=to('data/simu/'),
    # )

    # simu_ampl_rate(
    #     # initial sequence generation
    #
    #     len_params={
    #         'umi': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
    #         'umi_1': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
    #         'barcode': 16,
    #         'seq': 100,
    #         'seq_2': 100,
    #         'adapter': 10,
    #         'adapter_1': 10,
    #         'primer': 10,
    #         'primer_1': 10,
    #         'spacer': 10,
    #         'spacer_1': 10,
    #     },
    #     seq_params={
    #         'custom': 'AAGC',
    #         'custom_1': 'A',
    #     },
    #     seq_num=50,
    #     seq_len=100,
    #     working_dir=to('data/simu/'),
    #     fasta_cdna_fpn=False,
    #     # fasta_cdna_fpn=to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),
    #
    #     is_sv_umi_lib=True,
    #     is_sv_seq_lib=True,
    #     is_sv_primer_lib=True,
    #     is_sv_adapter_lib=True,
    #     is_sv_spacer_lib=True,
    #     condis=['umi'],
    #     # condis=['umi', 'seq'],
    #     # condis=['umi', 'custom', 'seq', 'custom_1'],
    #     sim_thres=3,
    #     permutation=0,
    #
    #     # PCR amplification
    #     ampl_rates=np.linspace(0.1, 1, 10),
    #     err_route='minnow',  # tree minnow err1d err2d mutation_table_minimum mutation_table_complete
    #     pcr_error=1e-4,
    #     pcr_num=10,
    #     err_num_met='nbinomial',
    #     seq_error=0.01,
    #     # seq_sub_spl_number=200, # None
    #     seq_sub_spl_rate=0.333,
    #     use_seed=True,
    #     seed=1,
    #
    #     verbose=False,  # True False
    #
    #     sv_fastq_fp=to('data/simu/'),
    # )

    # simu_umi_len(
    #     # initial sequence generation
    #
    #     len_params={
    #         'umi': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_lens': np.arange(7, 36 + 1, 1),
    #         },
    #         'umi_1': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
    #         'barcode': 16,
    #         'seq': 100,
    #         'seq_2': 100,
    #         'adapter': 10,
    #         'adapter_1': 10,
    #         'primer': 10,
    #         'primer_1': 10,
    #         'spacer': 10,
    #         'spacer_1': 10,
    #     },
    #     seq_params={
    #         'custom': 'AAGC',
    #         'custom_1': 'A',
    #     },
    #     seq_num=50,
    #     seq_len=100,
    #     working_dir=to('data/simu/'),
    #     fasta_cdna_fpn=False,
    #     # fasta_cdna_fpn=to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),
    #
    #     is_sv_umi_lib=True,
    #     is_sv_seq_lib=True,
    #     is_sv_primer_lib=True,
    #     is_sv_adapter_lib=True,
    #     is_sv_spacer_lib=True,
    #     condis=['umi'],
    #     # condis=['umi', 'seq'],
    #     # condis=['umi', 'custom', 'seq', 'custom_1'],
    #     sim_thres=3,
    #     permutation=0,
    #
    #     # PCR amplification
    #     ampl_rate=0.85,
    #     err_route='minnow',  # tree minnow err1d err2d mutation_table_minimum mutation_table_complete
    #     pcr_error=1e-4,
    #     pcr_num=10,
    #     err_num_met='nbinomial',
    #     seq_error=0.01,
    #     seq_sub_spl_number=200,  # None
    #     seq_sub_spl_rate=0.333,
    #     use_seed=True,
    #     seed=1,
    #
    #     verbose=False,  # True False
    #
    #     sv_fastq_fp=to('data/simu/'),
    # )

    simu_seq_dep(
        # initial sequence generation

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
        # condis=['umi', 'custom', 'seq', 'custom_1'],
        sim_thres=3,
        permutation=0,

        # PCR amplification
        ampl_rate=0.9,
        err_route='minnow',  # tree minnow err1d err2d mutation_table_minimum mutation_table_complete
        pcr_error=1e-04,
        pcr_num=10,
        err_num_met='nbinomial',
        seq_error=0.01,
        seq_sub_spl_numbers=[100, 500, 1000, 10000],  # None 200
        # seq_sub_spl_rate=0.333,
        use_seed=True,
        seed=1,

        verbose=False,  # True False

        sv_fastq_fp=to('data/simu/'),
    )