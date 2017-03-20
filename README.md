<h1 align="center"> knockout </h1>

[![Project Status](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/dirmeier/knockout.svg?branch=master)](https://travis-ci.org/dirmeier/knockout)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dirmeier/knockout?branch=master&svg=true)](https://ci.appveyor.com/project/dirmeier/knockout)
[![codecov](https://codecov.io/gh/dirmeier/knockout/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/knockout)

Analysis of high-throughput gene perturbation screens in R.

## Introduction

`knockout` does analysis of large-scale gene knockout/knockdown screens.
For this several preprocessing and data normalization techniques such as *median-polish* or *quantile-quantile normalization* are implemented. 
On the normalized data-set essential hits can be prioritized using state-of-the-art analysis tools, such as *gespeR* or *pmm*. 
The resulting hit lists can be further extended using network diffusion algorithms, such as *Markov random walks with restarts* or the well-known *heat equation*.
Eventually hits can be analyses using *GSEA*, etc.

## Installation
 
Install `knockout` using:
```{r}
devtools::install_github("dirmeier/knockout") 
```
from the R-console.

## Usage

Load the package using `library(knockout)`. We (maybe if I already was in the mood) provide a vignette for the package that can be called using: `vignette("knockout")`.
Basically that is all you have to know.

## References



## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@bsse.ethz.ch">simon.dirmeier@bsse.ethz.ch</a>
