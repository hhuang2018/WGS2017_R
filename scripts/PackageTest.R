#source("https://bioconductor.org/biocLite.R")
#biocLite("GeneGeneInteR")
library(GeneGeneInteR)

data("gene.pair")
head(gene.pair$Y)
gene.pair$G1
gene.pair$G2

# Multidimensional methods at the gene level
# 1. PCA
PCA.test(Y=gene.pair$Y, G1=gene.pair$G1, G2=gene.pair$G2,threshold=0.7,
         method="GenFreq")
PCA.test(Y=gene.pair$Y, G1=gene.pair$G1, G2=gene.pair$G2,threshold=0.7,
          method="Std")

# 2.  Canonical Correlation Analysis (CCA)
set.seed(1234)
CCA.test(Y=gene.pair$Y, G1=gene.pair$G1, G2=gene.pair$G2,n.boot=500)

# 3.  Kernel Canonical Correlation Analysis (KCCA)
set.seed(1234)
KCCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,
          kernel="rbfdot",sigma = 0.05,n.boot=500)

KCCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,
          kernel="polydot",degree = 1, scale = 1, offset = 1)

# 4. Partial Least Square Path Modeling (PLSPM)
set.seed(1234)
PLSPM.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,n.perm=1000)
          
# 5. Composite Linkage Disequilibrium (CLD)
CLD.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,n.perm=2000)

# 6. Gene-Based Information Gain Method (GBIGM)
GBIGM.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,n.perm=2000)


##  From SNP-SNP interaction to Gene-Gene interaction testing
# 1. minP
set.seed(1234)
minP.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2)

# 2. GATES
set.seed(1234)
gates.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,me.est="ChevNy")

set.seed(1234)
gates.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,alpha=0.05,me.est="Keff")

set.seed(1234)
gates.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,me.est="LiJi")

set.seed(1234)
gates.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,me.est="Galwey")

# 3. tTS and tProd
set.seed(1234)
tTS.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,tau=0.5,n.sim=10000)

set.seed(1234)
tProd.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,tau=0.05,n.sim=1000)
