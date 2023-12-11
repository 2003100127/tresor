__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import click
from pyfiglet import Figlet
from phylotres.locus import library as lib_sl
from phylotres.gene import library as lib_bulk
from phylotres.sc import library as lib_sc

from phylotres.locus import simu_seq_err as seqerr_sl
from phylotres.locus import simu_pcr_err as pcrerr_sl

from phylotres.util.Console import Console


vignette1 = Figlet(font='slant')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

console = Console()
console.verbose = True


class HelpfulCmd(click.Command):
    def format_help(self, ctx, formatter):
        click.echo(vignette1.renderText('PhyloTres'))
        click.echo(
            """
            tool 
                library_sl | library_bulk | library_sc
                
                @@@ library_sl
                phylotres library_sl -cfpn ./phylotres/data/libslocus.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
                
                @@@ library_bulk
                phylotres library_bulk -cfpn ./phylotres/data/libgene.yml -snum 50 -rfpn D:/Programming/R/R-4.3.2/ -nspl 2 -ngene 20 -bsimulator spsimseq -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True 

                @@@ library_sc
                phylotres library_sc -cfpn ./phylotres/data/libsc.yml -snum 50 -rfpn D:/Programming/R/R-4.3.2/ -ncell 10 -ngene 10 -bsimulator spsimseq -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True 

                @@@ seqerr_sl
                phylotres seqerr_sl -cfpn ./phylotres/data/seqerr_sl.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ pcrerr_sl
                phylotres pcrerr_sl -cfpn ./phylotres/data/pcrerr_sl.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
            
            """
        )


@click.command(cls=HelpfulCmd, context_settings=CONTEXT_SETTINGS)
@click.argument('tool', type=str)
@click.option(
    '-cfpn', '--config_fpn', type=str, required=True,
    help="""
        Path to a YMAL file
    """
)
@click.option(
    '-wd', '--working_dir', type=str, required=True,
    help="""
        Path to store results in the working directory
    """
)
@click.option(
    '-snum', '--seq_num', type=int,
    help="""
        read/UMI number
    """
)
@click.option(
    '-sthres', '--sim_thres', type=int,
    help="""
        Similarity threshold between UMIs
    """
)
@click.option(
    '-permut', '--permutation', type=int,
    help="""
        permutation 
    """
)
@click.option(
    '-md', '--mode', type=str,
    help="""
        short_read or long_read
    """
)
@click.option(
    '-is', '--is_seed', type=bool, default=True,
    help="""
        permutation 
    """
)
@click.option(
    '-isv_umi', '--is_sv_umi_lib', type=bool, default=True,
    help="""
        if is it save UMI library 
    """
)
@click.option(
    '-isv_seq', '--is_sv_seq_lib', type=bool, default=True,
    help="""
        if is it save sequence library 
    """
)
@click.option(
    '-isv_primer', '--is_sv_primer_lib', type=bool, default=True,
    help="""
        if is it save primer library 
    """
)
@click.option(
    '-isv_adapter', '--is_sv_adapter_lib', type=bool, default=True,
    help="""
        if is it save adapter library 
    """
)
@click.option(
    '-isv_spacer', '--is_sv_spacer_lib', type=bool, default=True,
    help="""
        if is it save spacer library 
    """
)
@click.option(### @@@ library bulk-RNA-seq simulation params
    '-rfpn', '--r_root', type=str,
    help="""
        R root directory path
    """
)
@click.option(
    '-nspl', '--num_samples', type=int,
    help="""
        number of samples
    """
)
@click.option(
    '-ngene', '--num_genes', type=int,
    help="""
        number of genes
    """
)
@click.option(
    '-bsimulator', '--bulk_simulator', type=str, default="spsimseq",
    help="""
        bulk simulator
    """
)

@click.option(
    '-ncell', '--num_cells', type=int,
    help="""
        number of cells
    """
)

@click.option(
    '-vb', '--verbose', type=bool, default=True,
    help="""
        Print verbose output
    """
)
def main(
        tool,
        config_fpn,
        working_dir,
        seq_num,

        # gspl
        r_root,
        num_samples,
        num_genes,
        bulk_simulator,

        # gmat
        num_cells,

        sim_thres,
        permutation,
        is_seed,
        is_sv_umi_lib,
        is_sv_seq_lib,
        is_sv_primer_lib,
        is_sv_adapter_lib,
        is_sv_spacer_lib,
        mode,
        verbose,
):
    print(vignette1.renderText('PhyloTres'))
    if tool == "library_sl":
        console.print("=============>Tool {} is being used...".format(tool))
        lib_sl(
            config_fpn=config_fpn,
            working_dir=working_dir,
            seq_num=seq_num,
            sim_thres=sim_thres,
            permutation=permutation,
            is_seed=is_seed,
            is_sv_umi_lib=is_sv_umi_lib,
            is_sv_seq_lib=is_sv_seq_lib,
            is_sv_primer_lib=is_sv_primer_lib,
            is_sv_adapter_lib=is_sv_adapter_lib,
            is_sv_spacer_lib=is_sv_spacer_lib,
            mode=mode,
            verbose=verbose,
        )
    elif tool == "library_bulk":
        console.print("=============>Tool {} is being used...".format(tool))
        lib_bulk(
            config_fpn=config_fpn,
            working_dir=working_dir,
            seq_num=seq_num,
            sim_thres=sim_thres,
            r_root=r_root,
            num_samples=num_samples,
            num_genes=num_genes,
            simulator=bulk_simulator,
            permutation=permutation,
            is_seed=is_seed,
            is_sv_umi_lib=is_sv_umi_lib,
            is_sv_seq_lib=is_sv_seq_lib,
            is_sv_primer_lib=is_sv_primer_lib,
            is_sv_adapter_lib=is_sv_adapter_lib,
            is_sv_spacer_lib=is_sv_spacer_lib,
            mode=mode,
            verbose=verbose,
        )
    elif tool == "library_sc":
        console.print("=============>Tool {} is being used...".format(tool))
        lib_sc(
            config_fpn=config_fpn,
            working_dir=working_dir,
            seq_num=seq_num,
            sim_thres=sim_thres,
            r_root=r_root,
            num_cells=num_cells,
            num_genes=num_genes,
            simulator=bulk_simulator,
            permutation=permutation,
            is_seed=is_seed,
            is_sv_umi_lib=is_sv_umi_lib,
            is_sv_seq_lib=is_sv_seq_lib,
            is_sv_primer_lib=is_sv_primer_lib,
            is_sv_adapter_lib=is_sv_adapter_lib,
            is_sv_spacer_lib=is_sv_spacer_lib,
            mode=mode,
            verbose=verbose,
        )
    elif tool == "seqerr_sl":
        console.print("=============>Tool {} is being used...".format(tool))
        seqerr_sl(
            config_fpn=config_fpn,
            working_dir=working_dir,
            seq_num=seq_num,
            sim_thres=sim_thres,
            permutation=permutation,
            is_seed=is_seed,
            is_sv_umi_lib=is_sv_umi_lib,
            is_sv_seq_lib=is_sv_seq_lib,
            is_sv_primer_lib=is_sv_primer_lib,
            is_sv_adapter_lib=is_sv_adapter_lib,
            is_sv_spacer_lib=is_sv_spacer_lib,
            mode=mode,
            verbose=verbose,
            sv_fastq_fp=working_dir,
        )
    elif tool == "pcrerr_sl":
        console.print("=============>Tool {} is being used...".format(tool))
        pcrerr_sl(
            config_fpn=config_fpn,
            working_dir=working_dir,
            seq_num=seq_num,
            sim_thres=sim_thres,
            permutation=permutation,
            is_seed=is_seed,
            is_sv_umi_lib=is_sv_umi_lib,
            is_sv_seq_lib=is_sv_seq_lib,
            is_sv_primer_lib=is_sv_primer_lib,
            is_sv_adapter_lib=is_sv_adapter_lib,
            is_sv_spacer_lib=is_sv_spacer_lib,
            mode=mode,
            verbose=verbose,
            sv_fastq_fp=working_dir,
        )