/**
 *  BinaryPoissonModel.h
 *  SPLITT
 *
 * Copyright 2018 Venelin Mitov
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
#ifndef BinaryPoissonModel_H_
#define BinaryPoissonModel_H_

#include "./SPLITT.h"
#include "./DiscreteTraitData.h"
#include <iostream>
#include <armadillo>
#include <cmath>

namespace BinaryPoissonUsingSPLITT {

using namespace SPLITT;

template<class Tree>
class BinaryPoissonModel: public TraversalSpecification<Tree> {

public:
  typedef BinaryPoissonModel<Tree> MyType;
  typedef TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef vec ParameterType;
  typedef DiscreteTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;

  // instantaneous rate parameters of transition from 0 to 1 and vice-versa
  double q01, q10;
  
  // rate matrix
  arma::mat22 Q; 
  arma::vec2 lambdaQ;
  arma::mat22 vecQ;
  arma::mat22 vecQ_1;
  
  // binary-trait vector at the tips of the tree.
  uvec x;
  
  // node-states 
  vec L0, L1;
  
  BinaryPoissonModel(TreeType const& tree, DataType const& input_data):
    BaseType(tree) {

    if(input_data.x_.size() != this->ref_tree_.num_tips()) {
      std::ostringstream oss;
      oss<<"The vector x must be the same length as the number of tips ("<<
        this->ref_tree_.num_tips()<<"), but were"<<input_data.x_.size()<<".";
      throw std::invalid_argument(oss.str());
    } else {
      
      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->x = At(input_data.x_, ordNodes);
      this->L0 = vec(this->ref_tree_.num_nodes());
      this->L1 = vec(this->ref_tree_.num_nodes());
      this->Q = arma::mat22();
    }
  };

  StateType StateAtRoot() const {
    vec res(2);
    res[0] = L0[this->ref_tree_.num_nodes() - 1];
    res[1] = L1[this->ref_tree_.num_nodes() - 1];
    return res;
  };
  
  void SetParameter(ParameterType const& par) {
    if(par.size() != 2) {
      throw std::invalid_argument(
          "The par vector should be of length 2 with \
      elements corresponding to q01 and q10.");
    }
    if(par[0] <= 0 || par[1] <= 0) {
      throw std::logic_error("The parameters q01 and q10 should be positive.");
    }
    this->q01 = par[0];
    this->q10 = par[1];
    
    this->Q(0,0) = -this->q01;
    this->Q(0,1) = this->q01;
    this->Q(1,0) = this->q10;
    this->Q(1,1) = -this->q10;
    
    using namespace arma;
    cx_vec2 cx_lambdaQ;
    cx_lambdaQ.fill(std::complex<double>(0.0, 0.0));
    cx_mat22 cx_vecQ;
    
    eig_gen(cx_lambdaQ, cx_vecQ, this->Q);
    
    this->lambdaQ = real(cx_lambdaQ);
    this->vecQ = real(cx_vecQ);
    this->vecQ_1 = inv(this->vecQ);
    
  }
  
  inline void InitNode(uint i) {
    if(i < this->ref_tree_.num_tips()) {
      L0[i] = x[i] == 0 ? 1.0 : 0.0;  
      L1[i] = x[i] == 1 ? 1.0 : 0.0;
    } else {
      L0[i] = L1[i] = 1.0;
    }
  }
  
  inline void VisitNode(uint i) {
    
  }
  
  inline void PruneNode(uint i, uint j) {
    double t = this->ref_tree_.LengthOfBranch(i);
    
    using namespace arma;
    // transition probability matrix
    mat22 Pt = this->vecQ * diagmat(exp(this->lambdaQ * t)) * this->vecQ_1;
    
    L0[j] = L0[j] * (Pt(0, 0)*L0[i] + Pt(0, 1)*L1[i]);
    L1[j] = L1[j] * (Pt(1, 0)*L0[i] + Pt(1, 1)*L1[i]);
  }
};
}

#endif //BinaryPoissonModel_H_
