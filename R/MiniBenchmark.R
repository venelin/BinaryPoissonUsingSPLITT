# MiniBenchmark.R
# BinaryPoissonUsingSPLITT
# 
# Copyright 2017 Venelin Mitov
# 
# This file is part of SPLITT: a generic C++ library for Serial and Parallel
# Lineage Traversal of Trees.
# 
# SPLITT is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# 
# SPLITT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with SPLITT.  If not, see
# <http://www.gnu.org/licenses/>.
# 
# @author Venelin Mitov

#' Perform a quick test of the PMM example likelihood calculation
#' @description This function runs a small benchmark to evaluate the SPLITT-package 
#' installation on a given computer. 
#' @param N number of tips in the test phylogenetic tree, default 10000
#' @param Ntests number of calculations within a call to sys.time (the resulting
#' times are averages Ntests calls). Default: 10.
#' @return a data.frame containing the response times for the PMM log-likelihood
#' calculation.
MiniBenchmark <- function(N = 10000, Ntests = 10) {
  
  # needed to pass the check
  speedupX <- NULL
  time.ms <- NULL
  
  cat("Performing a mini-benchmark of the PMM log-likelihood calculation with 
      a tree of size N=", 
      N, ";\nCalling each likelihood calculation Ntests=", 
      Ntests, " times ...\n")
  
  set.seed(10)
  N <- 100
  
  
  x0 <- 1
  q01 <- 0.02
  q10 <- 0.8

  tree <- rtree(N)
  x <- sample(c(0, 1), size = N, replace = TRUE)

    
  if(R.version[['os']]=='linux-gnu') {
    # this only works on linux
    cpuInfo <- system("cat /proc/cpuinfo | grep 'model name' | uniq", intern = TRUE)
  } else if(Sys.info()["sysname"] == "Darwin") {
    # this only works on mac OS x
    cpuInfo <- system("sysctl -a -n machdep.cpu.brand_string", intern = TRUE)
  } else {
    cpuInfo <- "(I don't know how to figure it out on your OS, sorry.)"
  }
  
  cat("CPU: ", cpuInfo, "\n")
  
  ord <- reorder(tree, order = "postorder", index.only = TRUE)
  
  cppPMMObject <- NewBinaryPoissonModelCppObject(x, tree)
  
  if(cppPMMObject$algorithm$VersionOPENMP == 0) {
    cat("It seems that OpenMP was disabled during C++ compilation. For parallel
        tree traversal OpenMP should be available on your system and supported
        by your C++ compiler. Only serial tree traversal is possible with 
        the current SPLITT installation. Please, read the user guide for further 
        instructions on how to build the SPLITT library.")
  } else {
    cat("OpenMP version: ", cppPMMObject$algorithm$VersionOPENMP, "\n")    
  }
  
  cat("Number of threads:",  cppPMMObject$algorithm$NumOmpThreads, "\n")
  
  cat("Measuring calculation times...\n")
  
  # warm-up for mode AUTO
  for(t in 1:1000) BinaryPoissonModelLogLikCpp(x, tree, x0, q01, q10, cppPMMObject, 0)
  
  tR <- system.time(for(t in seq_len(Ntests/10)) BinaryPoissonModelLogLik(x, tree, x0, q01, q10, ord = ord))[3] / (Ntests/10)
  unname(tR)
  
  measureTimeBinaryPoissonModelCpp <- function(mode) {
    unname(
      system.time(
        for(t in seq_len(Ntests)) 
          BinaryPoissonModelLogLikCpp(x, tree, x0, q01, q10, cppPMMObject, mode)
      )[3] / Ntests*1000
    )
  }
  
  resultsBinaryPoissonModel <- rbind(
    data.frame(model = "PMM", 
               mode = "R (serial)", 
               time.ms = tR*1000),
    data.frame(model = "PMM", 
               mode = "C++ (AUTO)", 
               time.ms = measureTimeBinaryPoissonModelCpp(0)),
    data.frame(model = "PMM", 
               mode = "C++ (SINGLE_THREAD_LOOP_POSTORDER)", 
               time.ms = measureTimeBinaryPoissonModelCpp(10)),
    data.frame(model = "PMM", 
               mode = "C++ (SINGLE_THREAD_LOOP_PRUNES)", 
               time.ms = measureTimeBinaryPoissonModelCpp(11)),
    data.frame(model = "PMM", 
               mode = "C++ (SINGLE_THREAD_LOOP_VISITS)", 
               time.ms = measureTimeBinaryPoissonModelCpp(12)),
    data.frame(model = "PMM", 
               mode = "C++ (MULTI_THREAD_LOOP_PRUNES)", 
               time.ms = measureTimeBinaryPoissonModelCpp(21)),
    data.frame(model = "PMM", 
               mode = "C++ (MULTI_THREAD_LOOP_VISITS)", 
               time.ms = measureTimeBinaryPoissonModelCpp(22)),
    data.frame(model = "PMM", 
               mode = "C++ (MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES)", 
               time.ms = measureTimeBinaryPoissonModelCpp(23)),
    data.frame(model = "PMM", 
               mode = "C++ (MULTI_THREAD_VISIT_QUEUE)", 
               time.ms = measureTimeBinaryPoissonModelCpp(24)),
    data.frame(model = "PMM", 
               mode = "C++ (MULTI_THREAD_LOOP_PRUNES_NO_EXCEPTION)", 
               time.ms = measureTimeBinaryPoissonModelCpp(25)),
    data.frame(model = "PMM", 
               mode = "C++ (HYBRID_LOOP_PRUNES)", 
               time.ms = measureTimeBinaryPoissonModelCpp(31)),
    data.frame(model = "PMM", 
               mode = "C++ (HYBRID_LOOP_VISITS)", 
               time.ms = measureTimeBinaryPoissonModelCpp(32)),
    data.frame(model = "PMM", 
               mode = "C++ (HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES)", 
               time.ms = measureTimeBinaryPoissonModelCpp(33))
  )
  
  rownames(resultsBinaryPoissonModel) <- NULL
  
  resultsBinaryPoissonModel
}
