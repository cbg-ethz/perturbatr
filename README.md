# perturbatr <img src="https://rawgit.com/cbg-ethz/perturbatr/master/inst/figure/sticker.png" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/cbg-ethz/perturbatr.svg?branch=master)](https://travis-ci.org/cbg-ethz/perturbatr)
[![Build app](https://ci.appveyor.com/api/projects/status/a28cs08ug9qng8hn?svg=true)](https://ci.appveyor.com/project/dirmeier/perturbatr)
[![codecov](https://codecov.io/gh/cbg-ethz/perturbatr/branch/master/graph/badge.svg)](https://codecov.io/gh/cbg-ethz/perturbatr)
[![bioc](https://bioconductor.org/shields/years-in-bioc/perturbatr.svg)](https://bioconductor.org/packages/release/bioc/html/perturbatr.html)

Analysis of high-throughput gene perturbation screens in R.

## Introduction

`perturbatr` does stage-wise analysis of large-scale genetic
perturbation screens for integrated data sets consisting of multiple screens.
For multiple integrated perturbation screens a hierarchical model that
considers the variance between different biological conditions is fitted.
That means that we first estimate relative effect sizes for all genes.
The resulting hit lists is then further extended using a network
propagation algorithm to correct for false negatives. and positives.

```{r}
data(rnaiscreen)
graph <- readRDS(
  system.file("extdata", "graph_file.tsv", package = "perturbatr"))

frm   <- Readout ~ Condition +
                   (1|GeneSymbol) + (1|Condition:GeneSymbol) +
                   (1|ScreenType) + (1|Condition:ScreenType)
ft    <- hm(rnaiscreen, formula = frm)
diffu <- diffuse(ft, graph=graph, r=0.3)

plot(diffu)
```

## Installation

You can install and use `perturbatr` either as an `R` library from [Bioconductor](https://doi.org/doi:10.18129/B9.bioc.perturbatr),
or by downloading the [tarball](https://github.com/cbg-ethz/perturbatr/releases).

If you want to use the **recommended** way using Bioconductor just call:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("perturbatr")
  
library(perturbatr)
```

from the R-console.

Installing the package using the downloaded tarball works like this:

```bash
  R CMD install <perturbatr.tar.gz>
```

where `perturbatr.tar.gz` is the downloaded tarball.

## Documentation

Load the package using `library(perturbatr)`. We provide a vignette for the package that can be called using: `vignette("perturbatr")`.

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@bsse.ethz.ch">simon.dirmeier@bsse.ethz.ch</a>
