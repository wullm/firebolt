/*******************************************************************************
 * This file is part of DEXM.
 * Copyright (c) 2020 Willem Elbers (whe@willemelbers.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef EVALUATE_H
#define EVALUATE_H

#include "grids.h"

double gridCIC(const double *box, int N, double boxlen, double x, double y, double z);
double evalDensity(const struct grids *grs, const struct multipoles *m,
                   double x, double y, double z, double qx, double qy,
                   double qz);
double evalDensityBin(const struct grids *grs, double x, double y, double z,
                      double nx, double ny, double nz, int index_q);
#endif
