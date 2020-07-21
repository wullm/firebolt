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

#ifndef FB_EVOLVE_H
#define FB_EVOLVE_H

int evolve_gsl(double **Psi, double q, double k, int l_max, double tau_ini,
               double tau_final, double mass, double c_vel, double dlnf0_dlnq,
               double (*redshift_func)(double log_tau),
               double (*h_prime_func)(double k, double log_tau),
               double (*eta_prime_func)(double k, double log_tau),
               double tolerance);

#endif
