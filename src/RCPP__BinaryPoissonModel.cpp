/**
  *  RCPP__BinaryPoissonModel.cpp
  *  SPLITT
  *
  * Copyright 2017 Venelin Mitov
  *
  * This file is part of SPLITT: a generic C++ library for Serial and Parallel
  * Lineage Traversal of Trees.
  *
  * SPLITT is free software: you can redistribute it and/or modify
  * it under the terms of the GNU Lesser General Public License as
  * published by the Free Software Foundation, either version 3 of
  * the License, or (at your option) any later version.
  *
  * SPLITT is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU Lesser General Public License for more details.
  *
  * You should have received a copy of the GNU Lesser General Public
  * License along with SPLITT.  If not, see
  * <http://www.gnu.org/licenses/>.
  *
  * @author Venelin Mitov
  */

#include <Rcpp.h>
#include "./BinaryPoissonModel.h"
    
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends("RcppArmadillo")]]

using namespace SPLITT;
using namespace BinaryPoissonUsingSPLITT;

typedef TraversalTask< BinaryPoissonModel<OrderedTree<uint, double>> > TraversalTaskBinaryPoissonModel;


TraversalTaskBinaryPoissonModel* CreateTraversalTaskBinaryPoissonModel(
    Rcpp::List const& tree, uvec const& values) {
  
  Rcpp::IntegerMatrix branches = tree["edge"];
  uvec parents(branches.column(0).begin(), branches.column(0).end());
  uvec daughters(branches.column(1).begin(), branches.column(1).end());
  vec t = Rcpp::as<vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  uvec tip_names = Seq(uint(1), num_tips);
  
  typename TraversalTaskBinaryPoissonModel::DataType data(tip_names, values);
  
  return new TraversalTaskBinaryPoissonModel(parents, daughters, t, data);
}


// This will enable returning a copy of the `TraversalAlgorithm`-object stored in
// a `TraversalTaskBinaryPoissonModel` object to a R. This will be used in the MiniBenchmark
// R-function to check things like the OpenMP version used during compilation and
// the number of OpenMP threads at runtime. 
RCPP_EXPOSED_CLASS_NODECL(TraversalTaskBinaryPoissonModel::AlgorithmType)
  
RCPP_MODULE(BinaryPoissonUsingSPLITT__TraversalTaskBinaryPoissonModel) {
  
  // Expose the properties VersionOPENMP and NumOmpThreads from the base 
  // TraversalAlgorithm class
  Rcpp::class_<TraversalTaskBinaryPoissonModel::AlgorithmType::ParentType> (
      "BinaryPoissonUsingSPLITT__BinaryPoissonModel__TraversalAlgorithm"
    )
  .property( "VersionOPENMP",
             &TraversalTaskBinaryPoissonModel::AlgorithmType::ParentType::VersionOPENMP )
  .property( "NumOmpThreads",
             &TraversalTaskBinaryPoissonModel::AlgorithmType::ParentType::NumOmpThreads )
  ;

  // Expose the TraversalTaskBinaryPoissonModel::AlgorithmType specifying that it derives 
  // from the base TraversalAlgorithm class
  Rcpp::class_<TraversalTaskBinaryPoissonModel::AlgorithmType> (
      "BinaryPoissonUsingSPLITT__BinaryPoissonModel__AlgorithmType"
    )
  .derives<TraversalTaskBinaryPoissonModel::AlgorithmType::ParentType>(
      "BinaryPoissonUsingSPLITT__BinaryPoissonModel__TraversalAlgorithm"
    )
  ;
  
  // Finally, expose the TraversalTaskBinaryPoissonModel class - this is the main class in 
  // the module, which will be instantiated from R using the factory function
  // we've just written.
  Rcpp::class_<TraversalTaskBinaryPoissonModel>( "BinaryPoissonUsingSPLITT__TraversalTaskBinaryPoissonModel" )
  // The <argument-type-list> MUST MATCH the arguments of the factory function 
  // defined above.
  .factory<Rcpp::List const&, uvec const&>( &CreateTraversalTaskBinaryPoissonModel )
  // Expose the method that we will use to execute the TraversalTask
  .method( "TraverseTree", &TraversalTaskBinaryPoissonModel::TraverseTree )
  // Expose the algorithm property
  .property( "algorithm", &TraversalTaskBinaryPoissonModel::algorithm )
  ;
}

