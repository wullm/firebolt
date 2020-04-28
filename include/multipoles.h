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

#ifndef MULTIPOLES_H
#define MULTIPOLES_H

#include "input.h"
#include "perturb_data.h"

struct multipoles {
  int k_size;
  int q_size;
  int l_size;
  double *Psi;
  double *k;
  double *q;
};


/* Fermi-Dirac distribution function */
double f0(double q);
double compute_dlnf0_dlnq(double q, double h);

/* Multipole functions */
int initMultipoles(struct multipoles *m, int k_size, int q_size, int l_size,
                   double q_min, double q_max, double k_min, double k_max);
int evolveMultipoles(struct multipoles *m, const struct perturb_data *ptdat,
                     double tau_ini, double tau_fin, double tol, double mass,
                     double c_vel, short verbose);
int convertMultipoleBasis_L2m(struct multipoles *mL, struct multipoles *mm, int l_max);
int cleanMultipoles(struct multipoles *m);



#endif
