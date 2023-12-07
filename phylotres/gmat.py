__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from phylotres.gcell.FromSimulator import fromSimulator as scgmatsimu
from phylotres.gcell.deepcvaesc.Predict import sample
from phylotres.gsample.FromSimulator import fromSimulator as bulkgmatsimu


#################################


def deepcvaesc_sc(

):
    sample(model_fpn, num_labels=num_labels, batch_size=32)
    return


def spsimseq_sc(
        R_root=None,
        num_genes=10,
        num_cells=10,
):
    gbycell, df_cells, df_genes = scgmatsimu(
        simulator='spsimseq',
        R_root=R_root,
        num_cells=num_genes,
        num_genes=num_cells,
    ).run()
    return gbycell, df_cells, df_genes


def spsimseq_bulk(
        R_root=None,
        num_samples=2,
        num_genes=10,
):
    gspl = bulkgmatsimu(
        simulator='spsimseq',
        R_root=R_root,
        num_samples=num_samples,
        num_genes=num_genes,
    ).run()
    return gspl


if __name__ == "__main__":
    # gbycell, _, _ = spsimseq_sc(
    #     R_root='D:/Programming/R/R-4.3.1/',
    #     num_genes=10,
    #     num_cells=10,
    # )

    # print(gbycell)

    gspl = spsimseq_bulk(
        R_root='D:/Programming/R/R-4.3.1/',
        num_samples=2,
        num_genes=10,
    )

    print(gspl)