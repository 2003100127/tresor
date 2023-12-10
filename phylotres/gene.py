__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from phylotres.library.Gene import Gene as libgene
from phylotres.scenario.seqerr.Gene import Gene as seqerr
from phylotres.scenario.pcrerr.Gene import Gene as pcrerr
from phylotres.scenario.pcrnum.Gene import Gene as pcrnum
from phylotres.scenario.amplrate.Gene import Gene as amplrate
from phylotres.scenario.umilen.Gene import Gene as umilen
from phylotres.scenario.seqdep.Gene import Gene as seqdep


def library(
        working_dir,
        seq_num,

        # gspl
        r_root,
        num_samples,
        num_genes,
        simulator,

        sim_thres,
        permutation,
        is_seed,
        is_sv_umi_lib,
        is_sv_seq_lib,
        is_sv_primer_lib,
        is_sv_adapter_lib,
        is_sv_spacer_lib,
        mode,

        len_params=None,
        seq_params=None,
        material_params=None,
        condis=None,

        config_fpn=None,
        verbose=True,
):
    if config_fpn:
        import yaml
        with open(config_fpn, "r") as f:
            configs = yaml.safe_load(f)
            # for k, item in configs.items():
            #     print(k, item)
        len_params = configs['len_params']
        seq_params = configs['seq_params']
        material_params = configs['material_params']
        condis = configs['condis']

    from phylotres.gsample.FromSimulator import fromSimulator

    gspl = fromSimulator(
        R_root=r_root,
        num_samples=num_samples,
        num_genes=num_genes,
        simulator=simulator,
    ).run()
    print(gspl)

    libgene(
        gspl=gspl,
        seq_num=seq_num,
        len_params=len_params,
        seq_params=seq_params,
        material_params=material_params,
        condis=condis,

        working_dir=working_dir,

        sim_thres=sim_thres,
        permutation=permutation,

        mode=mode,  # long_read short_read

        is_seed=is_seed,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,

        verbose=verbose,  # False True
    ).pooling()
    return 'Finished'


def simu_seq_err(
        gspl,
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
        gspl=gspl,
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
        **kwargs
    ).generate()
    return 'Finished!'


def simu_pcr_err(
        gspl,
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
        gspl=gspl,
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
        gspl,
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
        gspl=gspl,
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
        gspl,
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
        gspl=gspl,
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
        gspl,
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
        gspl=gspl,
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
        gspl,
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
        gspl=gspl,
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
    from phylotres.gmat import spsimseq_bulk
    import numpy as np

    gspl = spsimseq_bulk(
        R_root='D:/Programming/R/R-4.3.2/',
        num_genes=6,
        num_samples=2,
    )

    print(gspl)


    # print(simu_seq_err(
    #     gspl=gspl,
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
    #     ampl_rate=0.85,
    #     err_route='minnow',  # tree minnow err1d err2d mutation_table_minimum mutation_table_complete
    #     pcr_error=1e-4,
    #     pcr_num=10,
    #     err_num_met='nbinomial',
    #     seq_errors=[1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075,
    #                 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3],
    #     seq_sub_spl_number=200,  # None
    #     seq_sub_spl_rate=0.333,
    #     use_seed=True,
    #     seed=1,
    #
    #     verbose=True,  # True False
    #
    #     sv_fastq_fp=to('data/simu/'),
    # ))

    # simu_pcr_err(
    #     # initial sequence generation
    #     gspl=gspl,
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
    #     gspl=gspl,
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
    #     seq_sub_spl_number=200,  # None
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
    #     gspl=gspl,
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
    #     gspl=gspl,
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