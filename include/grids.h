/*******************************************************************************
 * This file is part of Firebolt.
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

#ifndef GRIDS_H
#define GRIDS_H

#include <fftw3.h>
#include "input.h"
#include "multipoles.h"

struct grids {
  double boxlen;
  int N; //the grid is N*N*N
  int q_size;
  int l_size;
  double *grids; //all the grids are stored in sequence
};

int initGrids(const struct params *pars, const struct multipoles *m,
              struct grids *grs);
int cleanGrids(struct grids *grs);

int generateGrids(const struct params *pars, const struct units *us,
                  const struct cosmology *cosmo,
                  const struct multipoles *m,
                  const fftw_complex *grf,
                  struct grids *grs);

#endif
