__version__ = "v1.0"
__copyright__ = "Copyright 2023"
__license__ = "MIT"
__lab__ = "cribbslab"

print(1)
import os
# os.environ['R_HOME'] = 'D:/Programming/anaconda3/envs/umi/Lib/R'
os.environ['R_HOME'] = 'D:/Programming/R/R-4.3.1/'
import rpy2.robjects as rob
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


class fromSimulator:

    def __init__(self, simulator):
        self.simulator = simulator

    def ifinstalled(self, ):
        from rpy2.robjects.packages import importr
        # try:
        #     utils = importr("SPsimSeq")
        #     print(utils)
        # except Exception:
        #     base = importr('base')
        #     print(base.source("http://www.bioconductor.org/biocLite.R"))
        #     # biocinstaller = importr("BiocInstaller")
        #     # biocinstaller.biocLite("SPsimSeq")
        #     print(1)
        return 1

    def SPsimSeqFixSM(self, ):
        self.ifinstalled()
        res = rob.r(
            '''
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
                                    group = treatment, n.genes = 10, batch.config = 1,
                                    group.config = c(0.5, 0.5), tot.samples = 10, 
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
            '''
        )
        a, b, c = res
        with localconverter(rob.default_converter + pandas2ri.converter):
            df = rob.conversion.rpy2py(a)
            df_cells = rob.conversion.rpy2py(b)
            df_genes = rob.conversion.rpy2py(c)
        df.columns = ['Cell_' + str(i) for i in range(10)]
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