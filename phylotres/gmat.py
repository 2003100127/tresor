__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from phylotres.gmat.FromSimulator import fromSimulator


def spsimseq_sc(
        R_root=None,
        num_genes=10,
        num_cells=10,
):
    gbycell, _, _ = fromSimulator(
        simulator='spsimseq',
        R_root=R_root,
        num_cells=num_genes,
        num_genes=num_cells,
    ).run()
    return gbycell, _, _



if __name__ == "__main__":
    gbycell = spsimseq_sc(
        R_root='D:/Programming/R/R-4.3.1/',
        num_genes=10,
        num_cells=10,
    )