__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import click
from pyfiglet import Figlet
from phylotres.library.SingleLocus import SingleLocus as libslocus
from phylotres.library.Gene import Gene as libgene
from phylotres.library.SingleCell import SingleCell as libsc
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
                
                library_sl
                phylotres library_sl -cfpn ./phylotres/data/libslocus.yml -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True
                
                @@@ library_bulk
                phylotres library_bulk -cfpn ./phylotres/data/libslocus.yml -rfpn D:/Programming/R/R-4.3.2/ -nspl 2 -ngene 20 -bsimulator spsimseq -permut 0 -sthres 3 -wd ./phylotres/data/simu/ -md short_read -is True -vb True 


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
    '-vb', '--verbose', type=bool, default=True,
    help="""
        Print verbose output
    """
)
def main(
        tool,
        config_fpn,
        working_dir,

        # gspl
        r_root,
        num_samples,
        num_genes,
        bulk_simulator,

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
        library_sl(
            config_fpn=config_fpn,
            working_dir=working_dir,
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
        library_bulk(
            config_fpn=config_fpn,
            working_dir=working_dir,
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


def library_sl(
        config_fpn,
        working_dir,
        sim_thres,
        permutation,
        is_seed,
        mode,
        is_sv_umi_lib,
        is_sv_seq_lib,
        is_sv_primer_lib,
        is_sv_adapter_lib,
        is_sv_spacer_lib,
        verbose,
):
    import yaml

    with open(config_fpn, "r") as f:
        configs = yaml.safe_load(f)
        # for k, item in configs.items():
        #     print(k, item)

    libslocus(
        seq_num=50,
        len_params=configs['len_params'],
        seq_params=configs['seq_params'],
        material_params=configs['material_params'],
        condis=configs['condis'],

        working_dir=working_dir,

        sim_thres=sim_thres,
        permutation=permutation,

        mode=mode,  # long_read short_read

        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,

        is_seed=is_seed,
        verbose=verbose,  # False True
    ).pooling()
    return 'Finished'


def library_bulk(
        config_fpn,
        working_dir,

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
        verbose,
):
    import yaml

    with open(config_fpn, "r") as f:
        configs = yaml.safe_load(f)
        # for k, item in configs.items():
        #     print(k,item)

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
        seq_num=50,
        len_params=configs['len_params'],
        seq_params=configs['seq_params'],
        material_params=configs['material_params'],
        condis=configs['condis'],

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