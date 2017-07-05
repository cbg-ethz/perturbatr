<h1 align="center"> knockout </h1>

[![Project Status](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/dirmeier/knockout.svg?branch=master)](https://travis-ci.org/dirmeier/knockout)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dirmeier/knockout?branch=master&svg=true)](https://ci.appveyor.com/project/dirmeier/knockout)
[![codecov](https://codecov.io/gh/dirmeier/knockout/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/knockout)

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
