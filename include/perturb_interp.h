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

#ifndef PERTURB_INTERP_H
#define PERTURB_INTERP_H

#include "perturb_data.h"

int rend_interp_init(const struct perturb_data *pt);
int rend_interp_switch_source(const struct perturb_data *pt, int index_src, int spline);
int rend_interp_free(const struct perturb_data *pt);
double rend_interp(double k, double tau, int spline);
// double rend_dydt_interp(double k, double tau, int spline);

#endif
