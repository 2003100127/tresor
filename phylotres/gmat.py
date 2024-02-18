__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from phylotres.gcell.FromSimulator import fromSimulator as scgmatsimu
from phylotres.gsample.FromSimulator import fromSimulator as bulkgmatsimu
from phylotres.util.file.Writer import Writer as pfwriter



def deepcvaesc_sc(

):
    from phylotres.gcell.deepcvaesc.Predict import sample
    sample(model_fpn, num_labels=num_labels, batch_size=32)
    return


def spsimseq_sc(
        R_root=None,
        num_genes=10,
        num_cells=10,
        simulator='spsimseq',
        sv_fpn=False,
):
    gbycell, df_cells, df_genes = scgmatsimu(
        simulator=simulator,
        R_root=R_root,
        num_cells=num_genes,
        num_genes=num_cells,
    ).run()
    if sv_fpn:
        # pfwriter().generic(df=gbycell, sv_fpn=sv_fpn)
        gbycell.to_hdf(sv_fpn, key='df', mode='w')
    return gbycell, df_cells, df_genes


def spsimseq_bulk(
        R_root=None,
        num_samples=2,
        num_genes=10,
        simulator='spsimseq',
        sv_fpn=False,
):
    gspl = bulkgmatsimu(
        simulator=simulator,
        R_root=R_root,
        num_samples=num_samples,
        num_genes=num_genes,
    ).run()
    if sv_fpn:
        # pfwriter().generic(df=gspl, sv_fpn=sv_fpn)
        gspl.to_hdf(sv_fpn, key='df', mode='w')
    return gspl


def read(
        gmat_fpn,
):
    import pandas as pd
    df_gmat = pd.read_hdf(gmat_fpn, 'df')
    if 'Y' in df_gmat.columns:
        df_gmat = df_gmat.drop(columns=['Y'])
    return df_gmat


if __name__ == "__main__":
    from phylotres.path import to

    # gbycell, _, _ = spsimseq_sc(
    #     R_root='D:/Programming/R/R-4.3.2/',
    #     num_genes=10,
    #     num_cells=10,
    #     simulator='spsimseq',
    #     sv_fpn=to('data/spsimseq_sc.h5'),
    # )
    # print(gbycell)

    # gspl = spsimseq_bulk(
    #     R_root='D:/Programming/R/R-4.3.2/',
    #     num_samples=2,
    #     num_genes=10,
    #     simulator='spsimseq',
    #     sv_fpn=to('data/spsimseq_bulk.h5'),
    # )
    # print(gspl)

    print(read(
        # gmat_fpn=to('data/gmat_customized.h5'),
        gmat_fpn=to('data/spsimseq_sc.h5'),
        # gmat_fpn=to('data/spsimseq_bulk.h5'),
    ))