
<!-- README.md is generated from README.Rmd. Please edit that file -->
dcanr: Differential co-expression/association network analysis
==============================================================

Application and evaluation of methods to infer differential co-expression/associations. An evaluation framework for testing differential co-expression methods in R. Simulated expression data with differential associations is used to assess inference methods. The package provides methods to evaluate inference methods and visualise predictions. Inference methods without implementations have been added to the package along with interfaces to existing implementations to enable ease of application and evaluations.

Installation
------------

You can install dcanr from github with:

``` r
# install.packages("devtools")
devtools::install_github("DavisLaboratory/dcanr")
```

Example
-------

This example shows how a differential network can be derived. Simulated data within the package is used.

``` r
library(dcanr)

#load simulated data
data(sim102)
#get expression data and conditions for 'UME6' knock-down
simdata <- getSimData(sim102, cond.name = 'UME6', full = FALSE)
emat <- simdata$emat
ume6_kd <- simdata$condition
#apply the z-score method with Spearman correlations
z_scores <- dcScore(emat, ume6_kd, cor.method = 'spearman')
#perform a statistical test: the z-test is selected automatically
raw_p <- dcTest(z_scores, emat, ume6_kd)
#adjust p-values (raw p-values from dcTest should NOT be modified)
adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')
#get the differential network
dcnet <- dcNetwork(z_scores, adj_p)
#> Warning in dcNetwork(z_scores, adj_p): default thresholds being selected
plot(dcnet, vertex.label = '')
```

![](README-example-1.png)
