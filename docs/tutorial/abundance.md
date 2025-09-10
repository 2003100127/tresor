
The `ts.abun.est_from_fasta` and `ts.abun.est_from_gtf` modules estimate the relative abundance of each transcript in order to proportionally sample molecules from a Fasta/GTF file. The total number of molecules is controlled by `--molecule-count`. It makes “one cell” abundance (the true molecule count setting for this sample).

- Sampling methods
    * [x] dirichlet
    * [x] lognormal
    * [x] zipf
    * [x] uniform
- Selecting methods
    * [x] random (uniform)
    * [x] first (first occurrence)


# Estimation from Fasta

Abundance can be estimated using a transcript fasta file through the following function.

Download the fasta file.

``` shell
wget -O gencode.v48.transcripts.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz
```

If you want to generate **50** transcripts from the file and the total expression value is **10000**, do

``` py
import tresor as ts

res = ts.abun.est_from_fasta(
    fasta_fpn="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/gencode.v48.transcripts.fa.gz",
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
    output_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.tsv",
    subset_ids_txt="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/refs.sub.ids.fa.txt",
    seed=1,
)
print("Top5:")
for i in range(min(5, len(res['ids']))):
    print(res['ids'][i], res['counts'][i])
```

The format of `abundance.tsv` file: tab-delimited, must contain at least Transcript ID and Abundance as two columns.

`abundance.tsv` looks like

``` shell
transcript	count
ENST00000769524	935
ENST00000816532	780
ENST00000820931	767
...
ENST00000754605	6
ENST00000672107	4
```

# Estimation from GTF

Download a GTF file from Genecode.

``` shell
wget -O gencode.v48.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
  
# wget -O gencode.v48.basic.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.basic.annotation.gtf.gz
```

``` py
import tresor as ts

res = ts.abun.est_from_gtf(
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
```

`abundance.pure.tsv` looks like

``` shell
ENST00000477535	4
ENST00000304081	3
ENST00000509522	2
ENST00000593404	1
ENST00000610462	0

```

`abundance.versioned.tsv` looks like

``` shell
ENST00000593404.5	878
ENST00000620405.1	732
ENST00000530097.1	720
...
ENST00000542002.5	317
ENST00000519650.5	286
```

# TPM
Simple TPM normalization

$$
TPM_i = \frac{count_i}{\sum count} \times 10^6
$$

!!! info "Explanation"

    * $TPM_i$: TPM value of the i-th transcript
    * $count_i$: Read count of the i-th transcript
    * $Σ(count)$: Sum of counts across all transcripts
    * $1e6$: Scaling factor (per million)

You can add a cell barcode (e.g., _ATGCATGCATGCATGC_) to the sampled transcripts in `cell_barcode` to denote from which cells they are derived.

``` py
import tresor as ts

res = ts.abun.to_tpm(
    counts_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.versioned.tsv",
    out_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.versioned.tpm.tsv",
    cell_barcode="ATGCATGCATGCATGC",
    include_header=False,
    float_format="{:.6f}"
)
```

`abundance.versioned.tpm.tsv` looks like

``` shell
ENST00000477535.5	400000.000000	ATGCATGCATGCATGC
ENST00000304081.9	300000.000000	ATGCATGCATGCATGC
ENST00000509522.6	200000.000000	ATGCATGCATGCATGC
ENST00000593404.5	100000.000000	ATGCATGCATGCATGC
ENST00000610462.1	0.000000	ATGCATGCATGCATGC
```