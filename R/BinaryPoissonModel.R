# BinaryPoissonModel.R
# BinaryPoissonUsingSPLITT
# 
# Copyright 2018 Venelin Mitov
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


#' Calculate the likelihood of a binary trait Poisson substituion model, 
#' given a tree and binary trait data at the tips. 
#' 
#' @description The binary trait Poisson substitution model assumes that a 
#' trait (character) with two possible states (0 and 1) changes it's state at 
#' random with constant rates of substitution from 0 to 1 (q01) and from
#' 1 to 0 (q10). The model parameters are the root-state x0, and the two rate
#' parameters q01 and q10. The likelihood for given model parameters <x0, q01, q10>
#' is defined as the probability of observing a given vector of states at the 
#' tips of a fixed phylogenetic tree. 
#' 
#' @param x an integer vector of size N, where N is the number of tips in tree.
#' @param tree a phylo object
#' @param x0,q01,q10 model parameters.
#' @param ord indices of the rows in tree$edge in pruning order. Defaults to 
#' \code{reorder(tree, order = "postorder", index.only = TRUE))}
#' 
#' @return the likelihood value.
BinaryPoissonModelLogLik <- function(
  x, tree, x0, q01, q10, 
  ord = reorder(tree, order = "postorder", index.only = TRUE)) {
  
  # number of tips in the tree
  N <- length(tree$tip.label)
  # total number of nodes in the tree (tips, internal nodes and root node)
  M <- nrow(tree$edge) + 1L
  # state variables for each node
  L0 <- L1 <- rep(1.0, M)
  
  # rate matrix
  Q <- rbind(c(-q01, q01),
             c(q10, -q10))
  
  eigQ <- eigen(Q)
  lambdaQ <- eigQ$values
  vecQ <- eigQ$vectors
  vecQ_1 <- solve(vecQ)
  
  for(o in ord) {
    # daughter node
    i <- tree$edge[o, 2]
    # parent node
    j <- tree$edge[o, 1]
    # branch length
    t <- tree$edge.length[o]
    
    Pt <- vecQ %*% diag(exp(lambdaQ * t)) %*% vecQ_1
    
    if(i <= N) {
      # initialize a tip node
      L0[i] <- as.double(x[i] == 0)
      L1[i] <- as.double(x[i] == 1)
    }
    
    L0[j] <- L0[j] * (Pt[1, 1]*L0[i] + Pt[1, 2]*L1[i])
    L1[j] <- L1[j] * (Pt[2, 1]*L0[i] + Pt[2, 2]*L1[i])
  }
  
  # for phylo objects, N+1 denotes the root node
  if(x0 == 0) {
    L0[N+1]
  } else {
    L1[N+1]
  }
}

#' Calculate the likelihood of a Poisson binary trait substituion model, 
#' given a tree and binary trait data at the tips.
#'  
#' @inheritParams BinaryPoissonModelLogLik
#' @param cppObject a previously created object returned by \code{\link{NewBinaryPoissonModelCppObject}}
#' @param mode an integer denoting the mode for traversing the tree, i.e. serial vs parallel.
#' 
#' @return the likelihood value.
BinaryPoissonModelLogLikCpp <- function(x, tree, x0, q01, q10, 
                         cppObject = NewBinaryPoissonModelCppObject(x, tree),
                         mode = getOption("SPLITT.postorder.mode", 0)) {
  L0L1 <- cppObject$TraverseTree(c(q01, q10), mode)
  
  if(x0 == 0) {
    L0L1[1]
  } else {
    L0L1[2]
  }
}


#' Create an instance of the Rcpp module for a given tree and trait data
#'
#' @inheritParams BinaryPoissonModelLogLik
#' @return an object to be passed as argument of the \link{BinaryPoissonModelLogLikCpp} function.
#' @seealso \link{BinaryPoissonModelLogLikCpp}
NewBinaryPoissonModelCppObject <- function(x, tree) {
  BinaryPoissonUsingSPLITT__TraversalTaskBinaryPoissonModel$new(tree, x[1:length(tree$tip.label)])
}
