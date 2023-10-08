__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import os
# os.environ['R_HOME'] = 'D:/Programming/anaconda3/envs/umi/Lib/R'
os.environ['R_HOME'] = 'D:/Programming/R/R-4.3.1/'
import rpy2.robjects as rob
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


class fromSimulator:

    def __init__(
            self,
            simulator,
            num_samples=2,
            num_genes=20,
    ):
        self.simulator = simulator
        self.num_samples = num_samples
        self.num_genes = num_genes

    def SPsimSeqFixSM(self, ):
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
        with localconverter(rob.default_converter + pandas2ri.converter):
            df = rob.conversion.rpy2py(res)
        return df.T

    def tool(self, ):
        return {
            'SPsimSeqFixSM': self.SPsimSeqFixSM()
        }

    def run(self, ):
        return self.tool()[self.simulator]