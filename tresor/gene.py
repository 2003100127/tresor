__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


import time
from tresor.library.Gene import Gene as libgene
from tresor.scenario.seqerr.Gene import Gene as seqerr
from tresor.scenario.pcrerr.Gene import Gene as pcrerr
from tresor.scenario.pcrnum.Gene import Gene as pcrnum
from tresor.scenario.amplrate.Gene import Gene as amplrate
from tresor.scenario.umilen.Gene import Gene as umilen
from tresor.scenario.seqdep.Gene import Gene as seqdep


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
        is_seed=True,
        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,

        len_params=None,
        seq_params=None,
        material_params=None,
        condis=None,
        mode=None,

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

    from tresor.gsample.FromSimulator import fromSimulator

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
        seq_num,
        working_dir,
        sim_thres,
        permutation,

        gspl=None,
        len_params=None,
        condis=None,

        ampl_rate=None,
        err_route=None,
        pcr_error=None,
        pcr_num=None,
        err_num_met=None,
        seq_errors=None,
        seq_sub_spl_number=None,
        seq_sub_spl_rate=None,

        use_seed=True,
        seed=1,

        pcr_deletion=False,
        pcr_insertion=False,
        pcr_del_rate=0,
        pcr_ins_rate=0,
        seq_deletion=False,
        seq_insertion=False,
        seq_del_rate=0,
        seq_ins_rate=0,

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        sv_fastq_fp=True,

        config_fpn=None,
        verbose=True,

        **kwargs,
):
    if config_fpn:
        import yaml
        with open(config_fpn, "r") as f:
            configs = yaml.safe_load(f)
            # for k, item in configs.items():
            #     print(k, item)
        len_params = configs['len_params']
        kwargs['seq_params'] = configs['seq_params']
        kwargs['material_params'] = configs['material_params']
        condis = configs['condis']

        if configs['gspl_fpn']:
            from tresor.gmat import read as gmatread
            gspl = gmatread(gmat_fpn=configs['gspl_fpn'])
            print(gspl)

        ampl_rate = configs['ampl_rate']
        err_route = configs['err_route']
        pcr_error = configs['pcr_error']
        pcr_num = configs['pcr_num']
        err_num_met = configs['err_num_met']
        seq_errors = configs['seq_errors']
        seq_sub_spl_number = configs['seq_sub_spl_number']
        seq_sub_spl_rate = configs['seq_sub_spl_rate']
        if "pcr_deletion" in configs.keys():
            pcr_deletion = configs['pcr_deletion']
        if "pcr_insertion" in configs.keys():
            pcr_insertion = configs['pcr_insertion']
        if "pcr_del_rate" in configs.keys():
            pcr_del_rate = configs['pcr_del_rate']
        if "pcr_ins_rate" in configs.keys():
            pcr_ins_rate = configs['pcr_ins_rate']
        if "seq_deletion" in configs.keys():
            seq_deletion = configs['seq_deletion']
        if "seq_insertion" in configs.keys():
            seq_insertion = configs['seq_insertion']
        if "seq_del_rate" in configs.keys():
            seq_del_rate = configs['seq_del_rate']
        if "seq_ins_rate" in configs.keys():
            seq_ins_rate = configs['seq_ins_rate']

    seqerr(
        # initial sequence generation
        gspl=gspl,
        len_params=len_params,
        mode=kwargs['mode'] if 'mode' in kwargs.keys() else None,
        material_params=kwargs['material_params'] if 'material_params' in kwargs.keys() else None,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        working_dir=working_dir,
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

        # indels
        pcr_deletion=pcr_deletion,
        pcr_insertion=pcr_insertion,
        pcr_del_rate=pcr_del_rate,
        pcr_ins_rate=pcr_ins_rate,
        seq_deletion=seq_deletion,
        seq_insertion=seq_insertion,
        seq_del_rate=seq_del_rate,
        seq_ins_rate=seq_ins_rate,

        use_seed=use_seed,
        seed=seed,
        verbose=verbose,
        sv_fastq_fp=sv_fastq_fp,
        # **kwargs
    ).generate()
    return 'Finished!'


def simu_pcr_err(
        seq_num,
        working_dir,
        sim_thres,
        permutation,

        gspl=None,
        len_params=None,
        condis=None,
        ampl_rate=None,
        err_route=None,
        pcr_errors=None,
        seq_error=None,
        pcr_num=None,
        err_num_met=None,

        use_seed=True,
        seed=1,

        pcr_deletion=False,
        pcr_insertion=False,
        pcr_del_rate=0,
        pcr_ins_rate=0,
        seq_deletion=False,
        seq_insertion=False,
        seq_del_rate=0,
        seq_ins_rate=0,

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        sv_fastq_fp=True,

        vis_num_err_per_read=False,

        seq_sub_spl_number=None,
        seq_sub_spl_rate=1/3,

        config_fpn=None,

        verbose=True,
        **kwargs,
):
    if config_fpn:
        import yaml
        with open(config_fpn, "r") as f:
            configs = yaml.safe_load(f)
            # for k, item in configs.items():
            #     print(k, item)
        len_params=configs['len_params']
        kwargs['seq_params']=configs['seq_params']
        kwargs['material_params']=configs['material_params']
        condis=configs['condis']

        if configs['gspl_fpn']:
            from tresor.gmat import read as gmatread
            gspl = gmatread(gmat_fpn=configs['gspl_fpn'])
            print(gspl)

        ampl_rate = configs['ampl_rate']
        err_route = configs['err_route']
        seq_error = configs['seq_error']
        pcr_num = configs['pcr_num']
        err_num_met = configs['err_num_met']
        pcr_errors = configs['pcr_errors']
        seq_sub_spl_number = configs['seq_sub_spl_number']
        seq_sub_spl_rate = configs['seq_sub_spl_rate']
        if "pcr_deletion" in configs.keys():
            pcr_deletion = configs['pcr_deletion']
        if "pcr_insertion" in configs.keys():
            pcr_insertion = configs['pcr_insertion']
        if "pcr_del_rate" in configs.keys():
            pcr_del_rate = configs['pcr_del_rate']
        if "pcr_ins_rate" in configs.keys():
            pcr_ins_rate = configs['pcr_ins_rate']
        if "seq_deletion" in configs.keys():
            seq_deletion = configs['seq_deletion']
        if "seq_insertion" in configs.keys():
            seq_insertion = configs['seq_insertion']
        if "seq_del_rate" in configs.keys():
            seq_del_rate = configs['seq_del_rate']
        if "seq_ins_rate" in configs.keys():
            seq_ins_rate = configs['seq_ins_rate']

    pcrerr(
        gspl=gspl,
        len_params=len_params,
        mode=kwargs['mode'] if 'mode' in kwargs.keys() else None,
        material_params=kwargs['material_params'] if 'material_params' in kwargs.keys() else None,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        pcr_errors=pcr_errors,
        seq_error=seq_error,
        pcr_num=pcr_num,
        err_num_met=err_num_met,

        # indels
        pcr_deletion=pcr_deletion,
        pcr_insertion=pcr_insertion,
        pcr_del_rate=pcr_del_rate,
        pcr_ins_rate=pcr_ins_rate,
        seq_deletion=seq_deletion,
        seq_insertion=seq_insertion,
        seq_del_rate=seq_del_rate,
        seq_ins_rate=seq_ins_rate,

        vis_num_err_per_read=vis_num_err_per_read,

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
        seq_num,
        working_dir,
        sim_thres,
        permutation,

        gspl=None,
        len_params=None,
        condis=None,
        ampl_rate=None,
        err_route=None,
        pcr_error=None,
        seq_error=None,
        pcr_nums=None,
        err_num_met=None,

        use_seed=True,
        seed=1,

        pcr_deletion=False,
        pcr_insertion=False,
        pcr_del_rate=0,
        pcr_ins_rate=0,
        seq_deletion=False,
        seq_insertion=False,
        seq_del_rate=0,
        seq_ins_rate=0,

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        sv_fastq_fp=True,

        seq_sub_spl_number=None,
        seq_sub_spl_rate=1 / 3,

        config_fpn=None,

        verbose=True,
        **kwargs,
):
    if config_fpn:
        import yaml
        with open(config_fpn, "r") as f:
            configs = yaml.safe_load(f)
            # for k, item in configs.items():
            #     print(k, item)
        len_params = configs['len_params']
        kwargs['seq_params'] = configs['seq_params']
        kwargs['material_params'] = configs['material_params']
        condis = configs['condis']

        if configs['gspl_fpn']:
            from tresor.gmat import read as gmatread
            gspl = gmatread(gmat_fpn=configs['gspl_fpn'])
            print(gspl)

        ampl_rate = configs['ampl_rate']
        err_route = configs['err_route']
        seq_error = configs['seq_error']
        pcr_nums = configs['pcr_nums']
        err_num_met = configs['err_num_met']
        pcr_error = configs['pcr_error']
        seq_sub_spl_number = configs['seq_sub_spl_number']
        seq_sub_spl_rate = configs['seq_sub_spl_rate']
        if "pcr_deletion" in configs.keys():
            pcr_deletion = configs['pcr_deletion']
        if "pcr_insertion" in configs.keys():
            pcr_insertion = configs['pcr_insertion']
        if "pcr_del_rate" in configs.keys():
            pcr_del_rate = configs['pcr_del_rate']
        if "pcr_ins_rate" in configs.keys():
            pcr_ins_rate = configs['pcr_ins_rate']
        if "seq_deletion" in configs.keys():
            seq_deletion = configs['seq_deletion']
        if "seq_insertion" in configs.keys():
            seq_insertion = configs['seq_insertion']
        if "seq_del_rate" in configs.keys():
            seq_del_rate = configs['seq_del_rate']
        if "seq_ins_rate" in configs.keys():
            seq_ins_rate = configs['seq_ins_rate']

    stime = time.time()
    res_dict = pcrnum(
        gspl=gspl,
        len_params=len_params,
        mode=kwargs['mode'] if 'mode' in kwargs.keys() else None,
        material_params=kwargs['material_params'] if 'material_params' in kwargs.keys() else None,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        pcr_error=pcr_error,
        seq_error=seq_error,
        pcr_nums=pcr_nums,
        err_num_met=err_num_met,

        # indels
        pcr_deletion=pcr_deletion,
        pcr_insertion=pcr_insertion,
        pcr_del_rate=pcr_del_rate,
        pcr_ins_rate=pcr_ins_rate,
        seq_deletion=seq_deletion,
        seq_insertion=seq_insertion,
        seq_del_rate=seq_del_rate,
        seq_ins_rate=seq_ins_rate,

        use_seed=use_seed,
        seed=seed,
        ampl_rate=ampl_rate,
        sv_fastq_fp=sv_fastq_fp,
        seq_sub_spl_number=seq_sub_spl_number,
        seq_sub_spl_rate=seq_sub_spl_rate,
        verbose=verbose,
    ).generate()
    etime = time.time()
    print("===>Time: {:.3f}s".format(etime - stime))
    print('Finished!')
    return res_dict


def simu_ampl_rate(
        seq_num,
        working_dir,
        sim_thres,
        permutation,

        gspl=None,
        len_params=None,
        condis=None,
        ampl_rates=None,
        err_route=None,
        pcr_error=None,
        seq_error=None,
        pcr_num=None,
        err_num_met=None,

        use_seed=True,
        seed=1,

        pcr_deletion=False,
        pcr_insertion=False,
        pcr_del_rate=0,
        pcr_ins_rate=0,
        seq_deletion=False,
        seq_insertion=False,
        seq_del_rate=0,
        seq_ins_rate=0,

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        sv_fastq_fp=True,

        seq_sub_spl_number=None,
        seq_sub_spl_rate=1 / 3,

        config_fpn=None,

        verbose=True,
        **kwargs,
):
    if config_fpn:
        import yaml
        with open(config_fpn, "r") as f:
            configs = yaml.safe_load(f)
            # for k, item in configs.items():
            #     print(k, item)
        len_params = configs['len_params']
        kwargs['seq_params'] = configs['seq_params']
        kwargs['material_params'] = configs['material_params']
        condis = configs['condis']

        if configs['gspl_fpn']:
            from tresor.gmat import read as gmatread
            gspl = gmatread(gmat_fpn=configs['gspl_fpn'])
            print(gspl)

        ampl_rates = configs['ampl_rates']
        err_route = configs['err_route']
        seq_error = configs['seq_error']
        pcr_num = configs['pcr_num']
        err_num_met = configs['err_num_met']
        pcr_error = configs['pcr_error']
        seq_sub_spl_number = configs['seq_sub_spl_number']
        seq_sub_spl_rate = configs['seq_sub_spl_rate']
        if "pcr_deletion" in configs.keys():
            pcr_deletion = configs['pcr_deletion']
        if "pcr_insertion" in configs.keys():
            pcr_insertion = configs['pcr_insertion']
        if "pcr_del_rate" in configs.keys():
            pcr_del_rate = configs['pcr_del_rate']
        if "pcr_ins_rate" in configs.keys():
            pcr_ins_rate = configs['pcr_ins_rate']
        if "seq_deletion" in configs.keys():
            seq_deletion = configs['seq_deletion']
        if "seq_insertion" in configs.keys():
            seq_insertion = configs['seq_insertion']
        if "seq_del_rate" in configs.keys():
            seq_del_rate = configs['seq_del_rate']
        if "seq_ins_rate" in configs.keys():
            seq_ins_rate = configs['seq_ins_rate']

    amplrate(
        gspl=gspl,
        len_params=len_params,
        mode=kwargs['mode'] if 'mode' in kwargs.keys() else None,
        material_params=kwargs['material_params'] if 'material_params' in kwargs.keys() else None,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        pcr_error=pcr_error,
        seq_error=seq_error,
        pcr_num=pcr_num,
        err_num_met=err_num_met,

        # indels
        pcr_deletion=pcr_deletion,
        pcr_insertion=pcr_insertion,
        pcr_del_rate=pcr_del_rate,
        pcr_ins_rate=pcr_ins_rate,
        seq_deletion=seq_deletion,
        seq_insertion=seq_insertion,
        seq_del_rate=seq_del_rate,
        seq_ins_rate=seq_ins_rate,

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
        seq_num,
        working_dir,
        sim_thres,
        permutation,

        gspl=None,
        len_params=None,
        condis=None,
        ampl_rate=None,
        err_route=None,
        pcr_error=None,
        seq_error=None,
        pcr_num=None,
        err_num_met=None,

        use_seed=True,
        seed=1,

        pcr_deletion=False,
        pcr_insertion=False,
        pcr_del_rate=0,
        pcr_ins_rate=0,
        seq_deletion=False,
        seq_insertion=False,
        seq_del_rate=0,
        seq_ins_rate=0,

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        sv_fastq_fp=True,

        seq_sub_spl_number=None,
        seq_sub_spl_rate=1 / 3,

        config_fpn=None,

        verbose=True,
        **kwargs,
):
    if config_fpn:
        import yaml
        with open(config_fpn, "r") as f:
            configs = yaml.safe_load(f)
            # for k, item in configs.items():
            #     print(k, item)
        len_params = configs['len_params']
        kwargs['seq_params'] = configs['seq_params']
        kwargs['material_params'] = configs['material_params']
        condis = configs['condis']

        if configs['gspl_fpn']:
            from tresor.gmat import read as gmatread
            gspl = gmatread(gmat_fpn=configs['gspl_fpn'])
            print(gspl)

        ampl_rate = configs['ampl_rate']
        err_route = configs['err_route']
        seq_error = configs['seq_error']
        pcr_num = configs['pcr_num']
        err_num_met = configs['err_num_met']
        pcr_error = configs['pcr_error']
        seq_sub_spl_number = configs['seq_sub_spl_number']
        seq_sub_spl_rate = configs['seq_sub_spl_rate']
        if "pcr_deletion" in configs.keys():
            pcr_deletion = configs['pcr_deletion']
        if "pcr_insertion" in configs.keys():
            pcr_insertion = configs['pcr_insertion']
        if "pcr_del_rate" in configs.keys():
            pcr_del_rate = configs['pcr_del_rate']
        if "pcr_ins_rate" in configs.keys():
            pcr_ins_rate = configs['pcr_ins_rate']
        if "seq_deletion" in configs.keys():
            seq_deletion = configs['seq_deletion']
        if "seq_insertion" in configs.keys():
            seq_insertion = configs['seq_insertion']
        if "seq_del_rate" in configs.keys():
            seq_del_rate = configs['seq_del_rate']
        if "seq_ins_rate" in configs.keys():
            seq_ins_rate = configs['seq_ins_rate']

    umilen(
        gspl=gspl,
        len_params=len_params,
        mode=kwargs['mode'] if 'mode' in kwargs.keys() else None,
        material_params=kwargs['material_params'] if 'material_params' in kwargs.keys() else None,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        ampl_rate=ampl_rate,
        pcr_error=pcr_error,
        seq_error=seq_error,
        pcr_num=pcr_num,

        # indels
        pcr_deletion=pcr_deletion,
        pcr_insertion=pcr_insertion,
        pcr_del_rate=pcr_del_rate,
        pcr_ins_rate=pcr_ins_rate,
        seq_deletion=seq_deletion,
        seq_insertion=seq_insertion,
        seq_del_rate=seq_del_rate,
        seq_ins_rate=seq_ins_rate,

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
        seq_num,
        working_dir,
        sim_thres,
        permutation,

        gspl=None,
        len_params=None,
        condis=None,
        ampl_rate=None,
        err_route=None,
        pcr_error=None,
        seq_error=None,
        pcr_num=None,
        err_num_met=None,

        use_seed=True,
        seed=1,

        pcr_deletion=False,
        pcr_insertion=False,
        pcr_del_rate=0,
        pcr_ins_rate=0,
        seq_deletion=False,
        seq_insertion=False,
        seq_del_rate=0,
        seq_ins_rate=0,

        is_sv_umi_lib=True,
        is_sv_seq_lib=True,
        is_sv_primer_lib=True,
        is_sv_adapter_lib=True,
        is_sv_spacer_lib=True,
        sv_fastq_fp=True,

        seq_sub_spl_numbers=None,
        seq_sub_spl_rate=1 / 3,

        config_fpn=None,

        verbose=True,
        **kwargs,
):
    if config_fpn:
        import yaml
        with open(config_fpn, "r") as f:
            configs = yaml.safe_load(f)
            # for k, item in configs.items():
            #     print(k, item)
        len_params = configs['len_params']
        kwargs['seq_params'] = configs['seq_params']
        kwargs['material_params'] = configs['material_params']
        condis = configs['condis']

        if configs['gspl_fpn']:
            from tresor.gmat import read as gmatread
            gspl = gmatread(gmat_fpn=configs['gspl_fpn'])
            print(gspl)

        ampl_rate = configs['ampl_rate']
        err_route = configs['err_route']
        seq_error = configs['seq_error']
        pcr_num = configs['pcr_num']
        err_num_met = configs['err_num_met']
        pcr_error = configs['pcr_error']
        seq_sub_spl_numbers = configs['seq_sub_spl_numbers']
        seq_sub_spl_rate = configs['seq_sub_spl_rate']
        if "pcr_deletion" in configs.keys():
            pcr_deletion = configs['pcr_deletion']
        if "pcr_insertion" in configs.keys():
            pcr_insertion = configs['pcr_insertion']
        if "pcr_del_rate" in configs.keys():
            pcr_del_rate = configs['pcr_del_rate']
        if "pcr_ins_rate" in configs.keys():
            pcr_ins_rate = configs['pcr_ins_rate']
        if "seq_deletion" in configs.keys():
            seq_deletion = configs['seq_deletion']
        if "seq_insertion" in configs.keys():
            seq_insertion = configs['seq_insertion']
        if "seq_del_rate" in configs.keys():
            seq_del_rate = configs['seq_del_rate']
        if "seq_ins_rate" in configs.keys():
            seq_ins_rate = configs['seq_ins_rate']

    seqdep(
        gspl=gspl,
        len_params=len_params,
        mode=kwargs['mode'] if 'mode' in kwargs.keys() else None,
        material_params=kwargs['material_params'] if 'material_params' in kwargs.keys() else None,
        seq_params=kwargs['seq_params'] if 'seq_params' in kwargs.keys() else None,
        seq_num=seq_num,
        is_sv_umi_lib=is_sv_umi_lib,
        is_sv_seq_lib=is_sv_seq_lib,
        is_sv_primer_lib=is_sv_primer_lib,
        is_sv_adapter_lib=is_sv_adapter_lib,
        is_sv_spacer_lib=is_sv_spacer_lib,
        working_dir=working_dir,
        condis=condis,
        sim_thres=sim_thres,
        permutation=permutation,
        err_route=err_route,
        pcr_error=pcr_error,
        seq_error=seq_error,
        pcr_num=pcr_num,
        err_num_met=err_num_met,

        # indels
        pcr_deletion=pcr_deletion,
        pcr_insertion=pcr_insertion,
        pcr_del_rate=pcr_del_rate,
        pcr_ins_rate=pcr_ins_rate,
        seq_deletion=seq_deletion,
        seq_insertion=seq_insertion,
        seq_del_rate=seq_del_rate,
        seq_ins_rate=seq_ins_rate,

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
    from tresor.path import to
    from tresor.gmat import spsimseq_bulk
    import numpy as np

    gspl = spsimseq_bulk(
        R_root='D:/Programming/R/R-4.3.2/',
        num_genes=6,
        num_samples=2,
    )
    print(gspl)

    # print(library(
    #     r_root='D:/Programming/R/R-4.3.2/',
    #     num_genes=6,
    #     num_samples=2,
    #     simulator='spsimseq',
    #
    #     seq_num=50,
    #     len_params={
    #         'umi': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
    #         'umi_1': {
    #             'umi_unit_pattern': 3,
    #             'umi_unit_len': 12,
    #         },
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
    #         'custom': 'BAGC',
    #         'custom_1': 'V',
    #     },
    #     material_params={
    #         'fasta_cdna_fpn': to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),  # None False
    #     },
    #     is_seed=True,
    #
    #     working_dir=to('data/simu/docs/'),
    #
    #     # condis=['umi'],
    #     # condis=['umi', 'seq'],
    #     condis=['umi', 'custom', 'seq', 'custom_1'],
    #
    #     # condis=['umi', 'primer', 'primer_1', 'spacer', 'spacer_1', 'adapter', 'adapter_1', 'seq', 'seq_2', 'umi_1'],
    #     sim_thres=3,
    #     permutation=0,
    #
    #     mode='short_read',  # long_read short_read
    #
    #     verbose=False,
    # ))

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
    #     material_params={
    #         'fasta_cdna_fpn': to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),  # None False
    #     },
    #     seq_num=50,
    #     working_dir=to('data/simu/docs/'),
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
    #     err_route='sptree',  # bftree sptree err1d err2d mutation_table_minimum mutation_table_complete
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
    #     verbose=False,  # True False
    #     mode='short_read',  # long_read short_read
    #
    #     sv_fastq_fp=to('data/simu/docs/'),
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
    #     material_params={
    #         'fasta_cdna_fpn': to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),  # None False
    #     },
    #     seq_num=50,
    #     working_dir=to('data/simu/docs/'),
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
    #     err_route='sptree',  # bftree sptree err1d err2d mutation_table_minimum mutation_table_complete
    #     pcr_errors=[1e-05, 2.5e-05, 5e-05, 7.5e-05, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.0025, 0.005, 0.0075, 0.01],
    #     pcr_num=10,
    #     err_num_met='nbinomial',
    #     seq_error=0.01,
    #     seq_sub_spl_number=200,  # None
    #     # seq_sub_spl_rate=0.333,
    #     use_seed=True,
    #     seed=1,
    #
    #     verbose=False,  # True False
    #     mode='short_read',  # long_read short_read
    #
    #     sv_fastq_fp=to('data/simu/docs/'),
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
    #     material_params={
    #         'fasta_cdna_fpn': to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),  # None False
    #     },
    #     seq_num=50,
    #     working_dir=to('data/simu/docs/'),
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
    #     err_route='sptree',  # bftree sptree err1d err2d mutation_table_minimum mutation_table_complete
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
    #     mode='short_read',  # long_read short_read
    #
    #     sv_fastq_fp=to('data/simu/docs/'),
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
    #     material_params={
    #         'fasta_cdna_fpn': to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),  # None False
    #     },
    #     seq_num=50,
    #     working_dir=to('data/simu/docs/'),
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
    #     err_route='sptree',  # bftree sptree err1d err2d mutation_table_minimum mutation_table_complete
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
    #     mode='short_read',  # long_read short_read
    #
    #     sv_fastq_fp=to('data/simu/docs/'),
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
    #     material_params={
    #         'fasta_cdna_fpn': to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),  # None False
    #     },
    #     seq_num=50,
    #     working_dir=to('data/simu/docs/'),
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
    #     err_route='sptree',  # bftree sptree err1d err2d mutation_table_minimum mutation_table_complete
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
    #     mode='short_read',  # long_read short_read
    #
    #     sv_fastq_fp=to('data/simu/docs/'),
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
        material_params={
            'fasta_cdna_fpn': to('data/Homo_sapiens.GRCh38.cdna.all.fa.gz'),  # None False
        },
        seq_num=50,
        working_dir=to('data/simu/docs/'),

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
        err_route='sptree',  # bftree sptree err1d err2d mutation_table_minimum mutation_table_complete
        pcr_error=1e-04,
        pcr_num=10,
        err_num_met='nbinomial',
        seq_error=0.01,
        seq_sub_spl_numbers=[100, 500, 1000, 10000],  # None 200
        # seq_sub_spl_rate=0.333,
        use_seed=True,
        seed=1,

        verbose=False,  # True False

        mode='short_read',  # long_read short_read

        sv_fastq_fp=to('data/simu/docs/'),
    )