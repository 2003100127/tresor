# 统一放到一个工作目录
WORK=/path/to/work
REF_FA=$WORK/ref/transcriptome.fa              # 转录本FASTA（如 GENCODE transcripts）
TX2GENE=$WORK/ref/tx2gene.tsv                  # transcript_id \t gene_id
WHITELIST=$WORK/ref/737K-august-2016.txt       # 10x v2/v3 白名单
INDEX_DIR=$WORK/mn_index
EST_DIR=$WORK/mn_estimate
SIM_DIR=$WORK/mn_sim
THREADS=16
READLEN=101   # 与 index -k 一致