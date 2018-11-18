/**
 *  DiscreteTraitData.h
 *  BinaryPoissonUsingSPLITT
 *
 * Copyright 2018 Venelin Mitov
 *
 * This file is part of BinaryPoissonUsingSPLITT.
 *
 * BinaryPoissonUsingSPLITT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * BinaryPoissonUsingSPLITT is distributed in the hope that it will be useful,
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
#ifndef DiscreteTraitData_H_
#define DiscreteTraitData_H_

#include "./SPLITT.h"

namespace BinaryPoissonUsingSPLITT {

using namespace SPLITT;

template<class NameType>
struct DiscreteTraitData {
  // use const references to avoid copying of long vectors
  std::vector<NameType> const& names_;
  uvec const& x_;
  DiscreteTraitData(
    std::vector<NameType> const& names,
    uvec const& x): names_(names), x_(x) {}
};
}
#endif //DiscreteTraitData_H_
