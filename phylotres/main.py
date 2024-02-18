__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import click
from pyfiglet import Figlet

### @@@ gmat_bulk gmat_sc
from phylotres.gmat import spsimseq_bulk as gmat_spsimseq_bulk
from phylotres.gmat import spsimseq_sc as gmat_spsimseq_sc

### library_sl | library_bulk | library_sc
from phylotres.locus import library as lib_sl
from phylotres.gene import library as lib_bulk
from phylotres.sc import library as lib_sc

### seqerr_sl | pcrerr_sl | pcrnum_sl | amplrate_sl | umilen_sl | seqdep_sl | generic_sl
from phylotres.locus import simu_seq_err as seqerr_sl
from phylotres.locus import simu_pcr_err as pcrerr_sl
from phylotres.locus import simu_pcr_num as pcrnum_sl
from phylotres.locus import simu_ampl_rate as amplrate_sl
from phylotres.locus import simu_umi_len as umilen_sl
from phylotres.locus import simu_seq_dep as seqdep_sl
from phylotres.locus import simu_generic as generic_sl

### seqerr_gene | pcrerr_gene | pcrnum_gene | amplrate_gene | umilen_gene | seqdep_gene
from phylotres.gene import simu_seq_err as seqerr_gene
from phylotres.gene import simu_pcr_err as pcrerr_gene
from phylotres.gene import simu_pcr_num as pcrnum_gene
from phylotres.gene import simu_ampl_rate as amplrate_gene
from phylotres.gene import simu_umi_len as umilen_gene
from phylotres.gene import simu_seq_dep as seqdep_gene

### seqerr_sc | pcrerr_sc | pcrnum_sc | amplrate_sc | umilen_sc | seqdep_sc
from phylotres.sc import simu_seq_err as seqerr_sc
from phylotres.sc import simu_pcr_err as pcrerr_sc
from phylotres.sc import simu_pcr_num as pcrnum_sc
from phylotres.sc import simu_ampl_rate as amplrate_sc
from phylotres.sc import simu_umi_len as umilen_sc
from phylotres.sc import simu_seq_dep as seqdep_sc

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
                gmat_bulk | gmat_sc
                
                @@@ gmat_bulk
                phylotres gmat_bulk -rfpn D:/Programming/R/R-4.3.2/ -nspl 2 -ngene 10 -gsimulator spsimseq -wd ./phylotres/data/spsimseq_bulk.h5 -is True -vb True
            
                @@@ gmat_sc
                phylotres gmat_sc -rfpn D:/Programming/R/R-4.3.2/ -ncell 10 -ngene 10 -gsimulator spsimseq -wd ./phylotres/data/spsimseq_sc.h5 -is True -vb True
            
            
                library_sl | library_bulk | library_sc
                
                @@@ library_sl
                phylotres library_sl -cfpn ./phylotres/data/libslocus.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
                
                @@@ library_bulk
                phylotres library_bulk -cfpn ./phylotres/data/libgene.yml -snum 50 -rfpn D:/Programming/R/R-4.3.2/ -nspl 2 -ngene 20 -gsimulator spsimseq -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True 

                @@@ library_sc
                phylotres library_sc -cfpn ./phylotres/data/libsc.yml -snum 50 -rfpn D:/Programming/R/R-4.3.2/ -ncell 10 -ngene 10 -gsimulator spsimseq -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True 

                
                seqerr_sl | pcrerr_sl | pcrnum_sl | amplrate_sl | umilen_sl | seqdep_sl | generic_sl
                
                @@@ seqerr_sl
                phylotres seqerr_sl -cfpn ./phylotres/data/seqerr_sl.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ pcrerr_sl
                phylotres pcrerr_sl -cfpn ./phylotres/data/pcrerr_sl.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ pcrnum_sl
                phylotres pcrnum_sl -cfpn ./phylotres/data/pcrnum_sl.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ amplrate_sl
                phylotres amplrate_sl -cfpn ./phylotres/data/amplrate_sl.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ umilen_sl
                phylotres umilen_sl -cfpn ./phylotres/data/umilen_sl.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ seqdep_sl
                phylotres seqdep_sl -cfpn ./phylotres/data/seqdep_sl.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ generic_sl
                phylotres generic_sl -cfpn ./phylotres/data/generic_sl.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
            
                seqerr_gene | pcrerr_gene | pcrnum_gene | amplrate_gene | umilen_gene | seqdep_gene
                @@@ seqerr_gene
                phylotres seqerr_gene -cfpn ./phylotres/data/seqerr_gene.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ pcrerr_gene
                phylotres pcrerr_gene -cfpn ./phylotres/data/pcrerr_gene.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ pcrnum_gene
                phylotres pcrnum_gene -cfpn ./phylotres/data/pcrnum_gene.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ amplrate_gene
                phylotres amplrate_gene -cfpn ./phylotres/data/amplrate_gene.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ umilen_gene
                phylotres umilen_gene -cfpn ./phylotres/data/umilen_gene.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ seqdep_gene
                phylotres seqdep_gene -cfpn ./phylotres/data/seqdep_gene.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
            
                seqerr_sc | pcrerr_sc | pcrnum_sc | amplrate_sc | umilen_sc | seqdep_sc
                @@@ seqerr_sc
                phylotres seqerr_sc -cfpn ./phylotres/data/seqerr_sc.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ pcrerr_sc
                phylotres pcrerr_sc -cfpn ./phylotres/data/pcrerr_sc.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ pcrnum_sc
                phylotres pcrnum_sc -cfpn ./phylotres/data/pcrnum_sc.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ amplrate_sc
                phylotres amplrate_sc -cfpn ./phylotres/data/amplrate_sc.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ umilen_sc
                phylotres umilen_sc -cfpn ./phylotres/data/umilen_sc.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                @@@ seqdep_sc
                phylotres seqdep_sc -cfpn ./phylotres/data/seqdep_sc.yml -snum 50 -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
            
                
            """
        )


@click.command(cls=HelpfulCmd, context_settings=CONTEXT_SETTINGS)
@click.argument('tool', type=str)
@click.option(
    '-cfpn', '--config_fpn', type=str,
    # required=True,
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
    '-gsimulator', '--gmat_simulator', type=str, default="spsimseq",
    help="""
        gmat simulator
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
        gmat_simulator,

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
    ### @@@ simu library
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
            simulator=gmat_simulator,
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
            simulator=gmat_simulator,
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

    ### @@@ simu gmat
    elif tool == "gmat_bulk":
        console.print("=============>Tool {} is being used...".format(tool))
        if gmat_simulator == 'spsimseq':
            gmat_spsimseq_bulk(
                R_root=r_root,
                num_samples=num_samples,
                num_genes=num_genes,
                simulator=gmat_simulator,
                sv_fpn=working_dir,
            )
    elif tool == "gmat_sc":
        console.print("=============>Tool {} is being used...".format(tool))
        if gmat_simulator == 'spsimseq':
            gmat_spsimseq_sc(
                R_root=r_root,
                num_cells=num_cells,
                num_genes=num_genes,
                simulator=gmat_simulator,
                sv_fpn=working_dir,
            )
    ### @@@ simu single locus
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
    elif tool == "pcrnum_sl":
        console.print("=============>Tool {} is being used...".format(tool))
        pcrnum_sl(
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
    elif tool == "amplrate_sl":
        console.print("=============>Tool {} is being used...".format(tool))
        amplrate_sl(
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
    elif tool == "umilen_sl":
        console.print("=============>Tool {} is being used...".format(tool))
        umilen_sl(
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
    elif tool == "seqdep_sl":
        console.print("=============>Tool {} is being used...".format(tool))
        seqdep_sl(
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
    elif tool == "generic_sl":
        console.print("=============>Tool {} is being used...".format(tool))
        generic_sl(
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

    ### @@@ simu bulk
    elif tool == "seqerr_gene":
        console.print("=============>Tool {} is being used...".format(tool))
        seqerr_gene(
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
    elif tool == "pcrerr_gene":
        console.print("=============>Tool {} is being used...".format(tool))
        pcrerr_gene(
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
    elif tool == "pcrnum_gene":
        console.print("=============>Tool {} is being used...".format(tool))
        pcrnum_gene(
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
    elif tool == "amplrate_gene":
        console.print("=============>Tool {} is being used...".format(tool))
        amplrate_gene(
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
    elif tool == "umilen_gene":
        console.print("=============>Tool {} is being used...".format(tool))
        umilen_gene(
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
    elif tool == "seqdep_gene":
        console.print("=============>Tool {} is being used...".format(tool))
        seqdep_gene(
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

    ### @@@ simu sc
    elif tool == "seqerr_sc":
        console.print("=============>Tool {} is being used...".format(tool))
        seqerr_sc(
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
    elif tool == "pcrerr_sc":
        console.print("=============>Tool {} is being used...".format(tool))
        pcrerr_sc(
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
    elif tool == "pcrnum_sc":
        console.print("=============>Tool {} is being used...".format(tool))
        pcrnum_sc(
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
    elif tool == "amplrate_sc":
        console.print("=============>Tool {} is being used...".format(tool))
        amplrate_sc(
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
    elif tool == "umilen_sc":
        console.print("=============>Tool {} is being used...".format(tool))
        umilen_sc(
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
    elif tool == "seqdep_sc":
        console.print("=============>Tool {} is being used...".format(tool))
        seqdep_sc(
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