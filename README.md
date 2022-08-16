# CAMO
Github repository for a molecular congruence analysis framework for evaluating model organisms (CAMO)


## Install CAMO package
CAMO package in be installed from Github by running the following code in R console.

```{R}
library(devtools)
install_github("https://github.com/CAMO-R/Rpackage") 
```
Alternatively, the complied source package CAMO_1.0.tar.gz can be downloaded at https://github.com/CAMO-R/other/tree/main/Package_dependences and installed locally by running the following code in R console.

```{R}
install.packages("~/CAMO_1.0.tar.gz", repos = NULL, type = "source")
```


Please make sure all dependency packages are installed. The code for installing dependency pakcages as follows:
```R
## from CRAN
CRAN.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
}
CRAN.packages(c("ROCR", "quantmod", "xts","RcppArmadillo", "Rcpp", "MASS", "parallel", "devtools", "methods", "igraph", "gridExtra", "grid", "ggplot2", "gplots", "reticulate"))

## from Bioconductor
Bioconductor.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        BiocManager::install(new.pkg, dependencies = TRUE)
}
Bioconductor.packages(c("DESeq2", "limma", "ConsensusClusterPlus", "pathview", "KEGGgraph", "KEGGREST", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dm.eg.db", "reactome.db"))
```

Some of the dependences that are not on CRAN and Bioconductor can be downloaded at https://github.com/CAMO-R/other/tree/main/Package_dependences.

MacOS users may encounter issues realted to C++ compiler in Big Sur. This tutorial might be helpful: http://yiqingxu.org/public/BigSurError.pdf


## Full tutorial (to be updated)


http://htmlpreview.github.io/?https://github.com/CAMO-R/Rpackage/blob/main/vignettes/CAMO.html

