is a module that can simulate reads consisting of only
UMIs per each, or UMI+Genomic sequence per each. The general-purpose
design gives the module this name. To achieve this purpose, a case-study
CLI should look like below:

## Gene-by-cell matrix
=== "Python"

    ``` py
    import phylotres as pts
    
    gbycell, _, _ = pts.spsimseq_sc(
        R_root='D:/Programming/R/R-4.3.2/',
        num_genes=10,
        num_cells=10,
        simulator='spsimseq',
        sv_fpn=to('data/spsimseq_sc.h5'),
    )
    print(gbycell)
    
    ```

=== "Shell"

    ```shell
    phylotres gmat_sc \
    -rfpn D:/Programming/R/R-4.3.2/ \ 
    -ncell 10 \ 
    -ngene 10 \ 
    -gsimulator spsimseq \ 
    -wd ./phylotres/data/spsimseq_sc.h5 \ 
    -is True \ 
    -vb True
    ```




## Output


``` py
12/12/2023 00:37:15 logger: =========>spsimseq is being used
SPsimSeq package version 1.12.0 
R[write to console]: Estimating featurewise correlations ...

R[write to console]: Selecting candidate DE genes ...

R[write to console]: Fitting zero probability model ...

R[write to console]: Estimating densities ...

R[write to console]: Configuring design ...

R[write to console]: Simulating data ...

R[write to console]:  ...1 of 1

12/12/2023 00:37:23 logger: =========>spsimseq completes simulation
12/12/2023 00:37:23 logger: =========>spsimseq simu result:
 $a
        Sample_1 Sample_2 Sample_3 Sample_4 Sample_5 Sample_6 Sample_7 Sample_8
Gene_1         6       66       15       16       28       12      253       57
Gene_2         0       96        0        0        0        0      122       38
Gene_3         0        0        0       40        0        0        0        0
Gene_4        52       34       20       15       18       14       14       43
Gene_5         0       98       37        0        0        0       17        0
Gene_6         0        0        0        0        0        0        0        0
Gene_7         0        0        0        0        0        0        0        0
Gene_8         0        0        0        0        0        0        0        0
Gene_9         0      188        7      153       23        0      201       20
Gene_10        1        0        0        0        0        0       31        0
        Sample_9 Sample_10
Gene_1         5         6
Gene_2         4         6
Gene_3         0         7
Gene_4         1         1
Gene_5         0         0
Gene_6         0         0
Gene_7         0         0
Gene_8         0         0
Gene_9         0        60
Gene_10        0         3

$b
          Batch Group sim.Lib.Size
Sample_1      1     1          543
Sample_2      1     1         1980
Sample_3      1     1          950
Sample_4      1     1         1098
Sample_5      1     1          711
Sample_6      1     2          147
Sample_7      1     2         1568
Sample_8      1     2         1425
Sample_9      1     2          933
Sample_10     1     2          475

$c
        DE.ind       source.ID
Gene_1    TRUE ENSG00000138160
Gene_2    TRUE ENSG00000138160
Gene_3   FALSE ENSG00000066322
Gene_4   FALSE ENSG00000101608
Gene_5   FALSE ENSG00000130066
Gene_6   FALSE ENSG00000196358
Gene_7   FALSE ENSG00000273674
Gene_8   FALSE ENSG00000269694
Gene_9   FALSE ENSG00000162992
Gene_10  FALSE ENSG00000142634


        Gene_1  Gene_2  Gene_3  Gene_4  ...  Gene_7  Gene_8  Gene_9  Gene_10
Cell_0     6.0     0.0     0.0    52.0  ...     0.0     0.0     0.0      1.0
Cell_1    66.0    96.0     0.0    34.0  ...     0.0     0.0   188.0      0.0
Cell_2    15.0     0.0     0.0    20.0  ...     0.0     0.0     7.0      0.0
Cell_3    16.0     0.0    40.0    15.0  ...     0.0     0.0   153.0      0.0
Cell_4    28.0     0.0     0.0    18.0  ...     0.0     0.0    23.0      0.0
Cell_5    12.0     0.0     0.0    14.0  ...     0.0     0.0     0.0      0.0
Cell_6   253.0   122.0     0.0    14.0  ...     0.0     0.0   201.0     31.0
Cell_7    57.0    38.0     0.0    43.0  ...     0.0     0.0    20.0      0.0
Cell_8     5.0     4.0     0.0     1.0  ...     0.0     0.0     0.0      0.0
Cell_9     6.0     6.0     7.0     1.0  ...     0.0     0.0    60.0      3.0

[10 rows x 10 columns]
```