# knockout <img src="https://cdn.rawgit.com/dirmeier/knockout/12ca810f/inst/figure/sticker.svg" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/cbg-ethz/knockout.svg?branch=master)](https://travis-ci.org/cbg-ethz/knockout)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/cbg-ethz/knockout?branch=master&svg=true)](https://ci.appveyor.com/project/cbg-ethz/knockout)
[![codecov](https://codecov.io/gh/cbg-ethz/knockout/branch/master/graph/badge.svg)](https://codecov.io/gh/cbg-ethz/knockout)

Analysis of high-throughput gene perturbation screens in R.

## Introduction

`knockout` does analysis of large-scale RNAi interference screens for pan-pathogenic datasets.
The package provides various tools for normalisation, plotting and analysis. For single pathogen
screens classical analyses using hypothesis testing are implemented. For pan-pathogenic 
screens we developed a random effects model that exploits the different biological settings
the data are derived from. The resulting hit lists can be further extended using network diffusion 
algorithms, such as *Markov random walks with restarts*.

## Installation
 
Download the latest `knockout` release and install the package using:

```bash
  R CMD install <netreg.tar.gz>
```
wher <netreg.tar.gz> is the downloaded tarball.

## Usage

Load the package using `library(knockout)`. We provide a vignette for the package that can be called using: `vignette("knockout")`. Basically that is all you have to know.

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@bsse.ethz.ch">simon.dirmeier@bsse.ethz.ch</a>
