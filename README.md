# knockdown <img src="https://rawgit.com/cbg-ethz/knockdown/master/inst/figure/sticker.png" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/cbg-ethz/knockdown.svg?branch=master)](https://travis-ci.org/cbg-ethz/knockdown)
[![codecov](https://codecov.io/gh/cbg-ethz/knockdown/branch/master/graph/badge.svg)](https://codecov.io/gh/cbg-ethz/knockdown)

Analysis of high-throughput gene perturbation screens in R.

## Introduction

`knockdown` does analysis of large-scale RNAi interference screens for pan-pathogenic datasets.
The package provides various tools for normalisation, plotting and analysis. For single pathogen
screens classical analyses using hypothesis testing are implemented. For pan-pathogenic
screens we developed a random effects model that exploits the different biological settings
the data are derived from. The resulting hit lists can be further extended using network diffusion
algorithms, such as *Markov random walks with restarts*.

## Installation

Download the latest `knockdown` release and install the package using:

```bash
  R CMD install <knockdown.tar.gz>
```
where `knockdown.tar.gz` is the downloaded tarball.

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@bsse.ethz.ch">simon.dirmeier@bsse.ethz.ch</a>
