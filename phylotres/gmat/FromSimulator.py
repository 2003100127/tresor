__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from phylotres.util.Console import Console


class fromSimulator:

    def __init__(
            self,
            simulator,
            R_root,
            num_cells=10,
            num_genes=10,
            verbose=True,
    ):
        self.simulator = simulator
        self.num_cells = num_cells
        self.num_genes = num_genes
        self.R_root = R_root

        self.console = Console()
        self.console.verbose = verbose

    def SPsimSeqFixSM(self, ):
        import os
        os.environ['R_HOME'] = self.R_root
        import rpy2.robjects as rob
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
        self.console.print('=========>SPsimSeqFixSM is being used')
        res = rob.r(
            """
            suppressPackageStartupMessages(library(SPsimSeq))
            cat("SPsimSeq package version", 
                as.character(packageVersion("SPsimSeq")), "\n")

            suppressPackageStartupMessages(library(SingleCellExperiment))
            # load the NGP nutlin data (availabl with the package, processed with 
            # SMARTer/C1 protocol, and contains read-counts)

            data("scNGP.data")
            # filter genes with sufficient expression level (important step) 
            scNGP.data2 <- scNGP.data[rowSums(counts(scNGP.data) > 0)>=5, ]  
            treatment <- ifelse(scNGP.data2$characteristics..treatment=="nutlin",2,1)
            set.seed(6543)
            scNGP.data2 <- scNGP.data2[sample(nrow(scNGP.data2), 20), ]
            # simulate data (we simulate here only a single data, n.sim = 1)
            sim.data.sc <- SPsimSeq(n.sim = 1, s.data = scNGP.data2,
                                    group = treatment, n.genes = {}, 
                                    """.format(self.num_genes)+"""
                                    batch.config = 1,
                                    group.config = c(0.5, 0.5), tot.samples = {},
                                    """.format(self.num_cells)+"""
                                    pDE = 0.2, lfc.thrld = 0.5, model.zero.prob = TRUE,
                                    result.format = "SCE")

            sim.data.sc1 <- sim.data.sc[[1]]
            class(sim.data.sc1)
            sum(counts(sim.data.sc1))
            dd <- list(
                "a"=data.frame(counts(sim.data.sc1)),
                "b"=data.frame(colData(sim.data.sc1)),
                "c"=data.frame(rowData(sim.data.sc1))
            )
            return (dd)
            """
        )
        self.console.print('=========>SPsimSeqFixSM completes simulation')
        a, b, c = res
        print(res)
        with localconverter(rob.default_converter + pandas2ri.converter):
            df = rob.conversion.rpy2py(a)
            df_cells = rob.conversion.rpy2py(b)
            df_genes = rob.conversion.rpy2py(c)
        df.columns = ['Cell_' + str(i) for i in range(self.num_cells)]
        df = df.T
        return df, df_cells, df_genes

    def tool(self, ):
        return {
            'SPsimSeqFixSM': self.SPsimSeqFixSM()
        }

    def run(self, ):
        return self.tool()[self.simulator]


if __name__ == "__main__":

    p = gmat, _, _ = fromSimulator(simulator='SPsimSeqFixSM').run()
    from scipy.sparse import coo_matrix
    csr_ = coo_matrix(gmat)
    print(csr_)