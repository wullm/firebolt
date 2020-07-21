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

#ifndef PERTURB_DATA_H
#define PERTURB_DATA_H

#include "input.h"

struct perturb_data {
  int k_size;
  int tau_size;
  int n_functions;
  double *delta;
  double *dydt; //time derivatives of delta
  double *k;
  double *log_tau;
  double *redshift;
  char **titles;
};

int readPerturb(struct params *pars, struct units *us, struct perturb_data *pt);
int computeDerivs(struct perturb_data *pt);
int cleanPerturb(struct perturb_data *pt);
double unitConversionFactor(char *title, double unit_length_factor,
                            double unit_time_factor);

#endif
