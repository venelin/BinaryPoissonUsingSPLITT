
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis build status](https://travis-ci.org/venelin/BinaryPoissonUsingSPLITT.svg?branch=master)](https://travis-ci.org/venelin/BinaryPoissonUsingSPLITT) [![Coverage status](https://codecov.io/gh/venelin/BinaryPoissonUsingSPLITT/branch/master/graph/badge.svg)](https://codecov.io/github/venelin/BinaryPoissonUsingSPLITT?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/BinaryPoissonUsingSPLITT?color=blue)](http://cran.r-project.org/web/packages/BinaryPoissonUsingSPLITT) [![Downloads](http://cranlogs.r-pkg.org/badges/BinaryPoissonUsingSPLITT?color=blue)](http://cran.rstudio.com/package=BinaryPoissonUsingSPLITT)

BinaryPoissonUsingSPLITT
========================

Likelihood calculation of the binary Poisson model using the SPLITT library. This package implementats the example of a discrete model from the article 'Parallel Likelihood Calculation for Phylogenetic Comparative Models: the 'SPLITT' C++ Library'. A preprint of the article is available from <https://doi.org/10.1101/235739>. The function 'MiniBenchmark' allows to compare the calculation times for different tree sizes. See [this guide](https://venelin.github.io/SPLITT/articles/SPLITTRcppModules.html) for a tutorial.

Installation
------------

You can install the released version of BinaryPoissonUsingSPLITT from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("BinaryPoissonUsingSPLITT")
```

Examples
--------

-   Calculating the likelihood of the PMM for a given tree, data and model parameters.

``` r
library(ape)
library(BinaryPoissonUsingSPLITT)
#> Loading required package: Rcpp

N <- 100

x <- 1
q01 <- 0.02
q10 <- 0.8

tree <- rtree(N)
x <- sample(c(0, 1), size = N, replace = TRUE)

cat("logLikelihood using R:", BinaryPoissonModelLogLik(x, tree, 1, 0.02, 0.8), "\n")
#> logLikelihood using R: 1.537735e-42
cat("logLikelihood using R:", BinaryPoissonModelLogLikCpp(x, tree, 1, 0.02, 0.8), "\n")
#> logLikelihood using R: 1.537735e-42
```

-   Performing a benchmark to measure the likelihood calculation times using different parallelization strategies:

``` r
# N specifies the size of the phylogenetic tree. 
# Ntests specifies the number of executions in one time measurement 
# (the more the better, but also slower).
MiniBenchmark(N = 100, Ntests = 100)
#> Performing a mini-benchmark of the PMM log-likelihood calculation with 
#>       a tree of size N= 100 ;
#> Calling each likelihood calculation Ntests= 100  times ...
#> CPU:  Intel(R) Core(TM) i7-4850HQ CPU @ 2.30GHz 
#> OpenMP version:  201107 
#> Number of threads: 8 
#> Measuring calculation times...
#>    model                                            mode time.ms
#> 1    PMM                                      R (serial)    2.90
#> 2    PMM                                      C++ (AUTO)    0.04
#> 3    PMM              C++ (SINGLE_THREAD_LOOP_POSTORDER)    0.03
#> 4    PMM                 C++ (SINGLE_THREAD_LOOP_PRUNES)    0.04
#> 5    PMM                 C++ (SINGLE_THREAD_LOOP_VISITS)    0.05
#> 6    PMM                  C++ (MULTI_THREAD_LOOP_PRUNES)    0.08
#> 7    PMM                  C++ (MULTI_THREAD_LOOP_VISITS)    0.09
#> 8    PMM C++ (MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES)    0.13
#> 9    PMM                  C++ (MULTI_THREAD_VISIT_QUEUE)    1.67
#> 10   PMM     C++ (MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION)    0.17
#> 11   PMM                        C++ (HYBRID_LOOP_PRUNES)    0.33
#> 12   PMM                        C++ (HYBRID_LOOP_VISITS)    0.36
#> 13   PMM       C++ (HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES)    0.27
MiniBenchmark(N = 1000, Ntests = 10)
#> Performing a mini-benchmark of the PMM log-likelihood calculation with 
#>       a tree of size N= 1000 ;
#> Calling each likelihood calculation Ntests= 10  times ...
#> CPU:  Intel(R) Core(TM) i7-4850HQ CPU @ 2.30GHz 
#> OpenMP version:  201107 
#> Number of threads: 8 
#> Measuring calculation times...
#>    model                                            mode time.ms
#> 1    PMM                                      R (serial)     3.0
#> 2    PMM                                      C++ (AUTO)     0.1
#> 3    PMM              C++ (SINGLE_THREAD_LOOP_POSTORDER)     0.1
#> 4    PMM                 C++ (SINGLE_THREAD_LOOP_PRUNES)     0.0
#> 5    PMM                 C++ (SINGLE_THREAD_LOOP_VISITS)     0.1
#> 6    PMM                  C++ (MULTI_THREAD_LOOP_PRUNES)     0.5
#> 7    PMM                  C++ (MULTI_THREAD_LOOP_VISITS)     0.2
#> 8    PMM C++ (MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES)     0.6
#> 9    PMM                  C++ (MULTI_THREAD_VISIT_QUEUE)     1.8
#> 10   PMM     C++ (MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION)     0.3
#> 11   PMM                        C++ (HYBRID_LOOP_PRUNES)     0.3
#> 12   PMM                        C++ (HYBRID_LOOP_VISITS)     0.3
#> 13   PMM       C++ (HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES)     0.3
MiniBenchmark(N = 10000, Ntests = 10)
#> Performing a mini-benchmark of the PMM log-likelihood calculation with 
#>       a tree of size N= 10000 ;
#> Calling each likelihood calculation Ntests= 10  times ...
#> CPU:  Intel(R) Core(TM) i7-4850HQ CPU @ 2.30GHz 
#> OpenMP version:  201107 
#> Number of threads: 8 
#> Measuring calculation times...
#>    model                                            mode time.ms
#> 1    PMM                                      R (serial)     3.0
#> 2    PMM                                      C++ (AUTO)     0.1
#> 3    PMM              C++ (SINGLE_THREAD_LOOP_POSTORDER)     0.0
#> 4    PMM                 C++ (SINGLE_THREAD_LOOP_PRUNES)     0.0
#> 5    PMM                 C++ (SINGLE_THREAD_LOOP_VISITS)     0.1
#> 6    PMM                  C++ (MULTI_THREAD_LOOP_PRUNES)     0.4
#> 7    PMM                  C++ (MULTI_THREAD_LOOP_VISITS)     0.2
#> 8    PMM C++ (MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES)     0.1
#> 9    PMM                  C++ (MULTI_THREAD_VISIT_QUEUE)     1.7
#> 10   PMM     C++ (MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION)     0.0
#> 11   PMM                        C++ (HYBRID_LOOP_PRUNES)     0.1
#> 12   PMM                        C++ (HYBRID_LOOP_VISITS)     0.4
#> 13   PMM       C++ (HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES)     0.2
```
