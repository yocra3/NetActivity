
# NetActivity

<!-- badges: start -->
<!-- badges: end -->

The goal of NetActivity is to compute gene set scores. Gene set scores are computed based on the weights obtained from training a sparsely-connected autoencoder.

## Installation

You can install the current release version of NetActivity from [Bioconductor](https://bioconductor.org/) with:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("NetActivity")
```

You can install the development version of NetActivity from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yocra3/NetActivity")
devtools::install_github("yocra3/NetActivityData")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(NetActivity)
library(airway)
data(airway)
ddsSE <- DESeq2::DESeqDataSet(airway, design = ~ cell + dex)
vst <- DESeq2::varianceStabilizingTransformation(ddsSE)
out <- prepareSummarizedExperiment(vst, "gtex_gokegg")
scores <- computeGeneSetScores(out, "gtex_gokegg")
```

