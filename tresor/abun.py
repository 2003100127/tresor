from tresor.tool import Abundance
from tresor.tool.TPM import TPM


def est_from_gtf(
        gtf_fpn,
        num_selected_mols,
        select_mode="random",
        num_total_mols=10,
        strip_version=True,
        dist="dirichlet",
        alpha=1.0,
        min_count=0,
        sort_output_by="count",
        include_noncoding=False,
        output_tsv="./abundance.pure.tsv",
        versioned_out="./abundance.versioned.tsv",
        subset_ids_txt="./refs.sub.ids.gtf.txt",
        seed=1,
):
    param_gtf = Abundance.AbundanceParams(
        gtf_fpn=gtf_fpn,
        num_selected_mols=num_selected_mols,
        select_mode=select_mode,
        num_total_mols=num_total_mols,
        strip_version=strip_version,
        dist=dist,
        alpha=alpha,
        min_count=min_count,
        sort_output_by=sort_output_by,
        include_noncoding=include_noncoding,
        output_tsv=output_tsv,
        versioned_out=versioned_out,
        subset_ids_txt=subset_ids_txt,
        seed=seed,
    )
    gen = Abundance.FromGTF(param_gtf)
    ids, counts = gen.build()
    gen.write(ids, counts)
    return {
        "ids": ids,
        "counts": counts,
    }


def est_from_fasta(
        fasta_fpn,
        num_selected_mols=50,
        select_mode="random",
        num_total_mols=10000,
        use_pipe_if_present=True,
        pipe_index=0,
        strip_version=True,
        dist="dirichlet",
        alpha=1.0,
        min_count=0,
        sort_output_by="count",
        output_tsv="./abundance.tsv",
        subset_ids_txt="./refs.sub.ids.fa.txt",
        seed=1,
):
    param_gtf = Abundance.AbundanceParams(
        fasta_fpn=fasta_fpn,
        output_tsv=output_tsv,
        num_selected_mols=num_selected_mols,
        select_mode=select_mode,
        num_total_mols=num_total_mols,
        use_pipe_if_present=use_pipe_if_present,
        pipe_index=pipe_index,
        strip_version=strip_version,
        dist=dist,
        alpha=alpha,
        min_count=min_count,
        sort_output_by=sort_output_by,
        subset_ids_txt=subset_ids_txt,
        seed=seed,
    )
    gen = Abundance.FromFasta(param_gtf)
    ids, counts = gen.build()
    gen.write(ids, counts)
    return {
        "ids": ids,
        "counts": counts,
    }


def to_tpm(
        counts_tsv,
        out_tsv,
        cell_barcode,
        include_header=False,
        float_format="{:.6f}",
):
    return TPM(
        counts_tsv=counts_tsv,
        out_tsv=out_tsv,
        cell_barcode=cell_barcode,
        include_header=include_header,
        float_format=float_format,
    ).convert_from_cnt()




if __name__ == "__main__":
    # res = est_from_fasta(
    #     fasta_fpn="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/gencode.v48.transcripts.fa.gz",
    #     num_selected_mols=50,
    #     select_mode="random",
    #     num_total_mols=10000,
    #     use_pipe_if_present=True,
    #     pipe_index=0,
    #     strip_version=True,
    #     dist="dirichlet",
    #     alpha=1.0,
    #     min_count=0,
    #     sort_output_by="count",
    #     output_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.tsv",
    #     subset_ids_txt="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/refs.sub.ids.fa.txt",
    #     seed=1,
    # )
    # print("Top5:")
    # for i in range(min(5, len(res['ids']))):
    #     print(res['ids'][i], res['counts'][i])

    res = est_from_gtf(
        gtf_fpn="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/gencode.v48.annotation.gtf",
        num_selected_mols=5,
        select_mode="random",
        num_total_mols=10,
        strip_version=True,
        dist="dirichlet",
        alpha=1.0,
        min_count=0,
        sort_output_by="count",
        include_noncoding=False,
        output_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.pure.tsv",
        versioned_out="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.versioned.tsv",
        subset_ids_txt="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/refs.sub.ids.gtf.txt",
        seed=1,
    )
    print("Top5:")
    for i in range(min(5, len(res['ids']))):
        print(res['ids'][i], res['counts'][i])

    to_tpm(
        counts_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.versioned.tsv",
        out_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.tpm.tsv",
        cell_barcode="ATGCATGCATGCATGC",
        include_header=False,
        float_format="{:.6f}"
    )