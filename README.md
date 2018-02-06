# perturbatr <img src="https://rawgit.com/cbg-ethz/perturbatr/master/inst/figure/sticker.png" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/cbg-ethz/perturbatr.svg?branch=master)](https://travis-ci.org/cbg-ethz/perturbatr)
[![codecov](https://codecov.io/gh/cbg-ethz/perturbatr/branch/master/graph/badge.svg)](https://codecov.io/gh/cbg-ethz/perturbatr)
[![bioc](https://bioconductor.org/shields/years-in-bioc/perturbatr.svg)](https://bioconductor.org/packages/release/bioc/html/perturbR.html)

Analysis of high-throughput gene perturbation screens in R.

## Introduction

`perturbatr` does stage-wise analysis of large-scale genetic perturbation screens for integrated datasets consisting of several screens. The package provides various tools for normalisation, plotting and analysis. For single perturbation screens classical analyses using hypothesis testing are implemented. For multiple integrated perturbation screens we developed a hierarchical model that considers the variance between the different biological settings the data are derived from. Our model first estimates an overall relative size for all genes (hits). The resulting hit lists is then be further extended using a network propagation algorithm to correct for false negatives and positives.

## Installation

Download the latest `perturbatr` release and install the package using:

```bash
  R CMD install <perturbatr.tar.gz>
```
where `perturbatr.tar.gz` is the downloaded tarball.

## Documentation

Load the package using `library(perturbatr)`. We provide a vignette for the package that can be called using: `vignette("perturbatr")`.

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@bsse.ethz.ch">simon.dirmeier@bsse.ethz.ch</a>
