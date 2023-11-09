/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_CUBOID_2D_HH
#define LATTICE_CUBOID_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticeCuboid2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T,typename DESCRIPTOR>
SuperLatticeCuboid2D<T,DESCRIPTOR>::SuperLatticeCuboid2D(
  SuperLattice<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "cuboid";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeCuboid2D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC),
                                this->_sLattice.getLoadBalancer().glob(iC)) );
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticeCuboid2D<T,DESCRIPTOR>::BlockLatticeCuboid2D
(BlockLattice<T,DESCRIPTOR>& blockLattice, const int iC)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _iC(iC)
{
  this->getName() = "cuboid";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeCuboid2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = _iC + 1;
  return false;
}

}
#endif
