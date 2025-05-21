__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


from tresor.util.Console import Console


class fromSimulator:

    def __init__(
            self,
            simulator,
            R_root,
            num_samples=2,
            num_genes=20,
            verbose=True,
    ):
        self.simulator = simulator
        self.num_samples = num_samples
        self.num_genes = num_genes
        self.R_root = R_root

        self.console = Console()
        self.console.verbose = verbose

    def spsimseq(self, ):
        import os
        os.environ['R_HOME'] = self.R_root
        import rpy2.robjects as rob
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
        self.console.print('=========>spsimseq is being used')
        res = rob.r(
            """
            suppressPackageStartupMessages(library(SPsimSeq))
            cat("SPsimSeq package version", as.character(packageVersion("SPsimSeq")), "\n")
            
            # load the Zhang bulk RNA-seq data
            data("zhang.data.sub")
            # filter genes with sufficient expression (important step) 
            zhang.counts <- zhang.data.sub$counts[rowSums(zhang.data.sub$counts > 0)>=5, ]
            
            set.seed(6452)
            zhang.counts2 <- zhang.counts[sample(nrow(zhang.counts), 20), ]
            sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts2,
                                      group = zhang.data.sub$MYCN.status, 
                                      n.genes = {},
            """.format(self.num_genes)+""" 
                                      batch.config = 1,
                                      group.config = c(0.5, 0.5), tot.samples = {},
                                      """.format(self.num_samples)+"""
                                      pDE = 0.5, lfc.thrld = 0.5, 
                                      result.format = "list")
            sim.data.bulk1 <- sim.data.bulk[[1]]                              
            return (data.frame(sim.data.bulk1$counts))
            """
        )
        self.console.print('=========>spsimseq completes simulation')
        with localconverter(rob.default_converter + pandas2ri.converter):
            df = rob.conversion.rpy2py(res)
        return df.T

    def tool(self, ):
        return {
            'spsimseq': self.spsimseq()
        }

    def run(self, ):
        return self.tool()[self.simulator]