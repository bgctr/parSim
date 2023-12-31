/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef LATTICE_PHYS_BOUNDARY_FORCE_3D_HH
#define LATTICE_PHYS_BOUNDARY_FORCE_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysBoundaryForce3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticePhysBoundaryForce3D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce3D(
  SuperLattice<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "physBoundaryForce";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysBoundaryForce3D<T, DESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        _indicatorF->getBlockIndicatorF(iC),
        this->_converter));
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysBoundaryForce3D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce3D(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  SuperGeometry<T,3>& superGeometry, const int material,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysBoundaryForce3D(sLattice,
                                    superGeometry.getMaterialIndicator(material),
                                    converter)
{ }

template <typename T, typename DESCRIPTOR>
BlockLatticePhysBoundaryForce3D<T,DESCRIPTOR>::BlockLatticePhysBoundaryForce3D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockIndicatorF3D<T>&                  indicatorF,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice, converter, 3),
    _indicatorF(indicatorF),
    _blockGeometry(indicatorF.getBlockGeometry())
{
  this->getName() = "physBoundaryForce";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysBoundaryForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T();
  }

  if (_indicatorF(input)) {
    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
      // Get direction
      const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
      // Get next cell located in the current direction
      // Check if the next cell is a fluid node
      if (_blockGeometry.get({input[0] + c[0], input[1] + c[1], input[2] + c[2]}) == 1) {
        // Get f_q of next fluid cell where l = opposite(q)
        T f = this->_blockLattice.get(input[0] + c[0], input[1] + c[1], input[2] + c[2])[iPop];
        // Get f_l of the boundary cell
        // Add f_q and f_opp
        f += this->_blockLattice.get(input)[descriptors::opposite<DESCRIPTOR>(iPop)];
        // Update force
        for (int i = 0; i < this->getTargetDim(); ++i) {
          output[i] -= c[i] * f;
        }
      }
    }
    for (int i = 0; i < this->getTargetDim(); ++i) {
      output[i] = this->_converter.getPhysForce(output[i]);
    }
  }
  return true;
}

}
#endif
