/*  This file is part of the OpenLB library
 *
 *  Copyright (C) Adrian Kummerlaender
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

#ifndef BLOCK_MAX_3D_H
#define BLOCK_MAX_3D_H

#include "blockBaseF3D.h"
#include "geometry/cuboid3D.h"
#include "indicator/blockIndicatorBaseF3D.h"

namespace olb {


/// BlockMax3D returns the max in each component of f on a indicated subset
template <typename T, typename W = T>
class BlockMax3D final : public BlockF3D<W> {
private:
  BlockF3D<W>&          _f;
  BlockIndicatorF3D<T>& _indicatorF;
  Cuboid3D<T>&          _cuboid;
public:
  BlockMax3D(BlockF3D<W>&          f,
             BlockIndicatorF3D<T>& indicatorF,
             Cuboid3D<T>&          cuboid);
  bool operator() (W output[], const int input[]) override;
};


}

#endif
