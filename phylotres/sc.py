__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from phylotres.simulate.seqerr.SingleCell import SingleCell as seqerrsc


def simu_seq_err(
        gmat,
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
    seqerrsc(
        # initial sequence generation
        gmat=gmat,
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


if __name__ == "__main__":
    from phylotres.path import to
    from phylotres.gmat import spsimseq_sc

    gmat, _, _ = spsimseq_sc(
        R_root='D:/Programming/R/R-4.3.1/',
        num_genes=6,
        num_cells=2,
    )

    print(gmat)

    print(simu_seq_err(
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
        # condis=['umi'],
        # condis=['umi', 'seq'],
        condis=['barcode', 'umi'],
        # condis=['barcode', 'custom', 'umi', 'custom_1'],
        sim_thres=3,
        permutation=0,

        # PCR amplification
        ampl_rate=0.85,
        err_route='minnow',  # tree minnow err1d err2d mutation_table_minimum mutation_table_complete
        pcr_error=1e-4,
        pcr_num=10,
        err_num_met='nbinomial',
        seq_errors=[1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075,
                    0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3],
        seq_sub_spl_number=200,  # None
        seq_sub_spl_rate=0.333,
        use_seed=True,
        seed=1,

        verbose=True,  # True False

        sv_fastq_fp=to('data/simu/'),
    ))