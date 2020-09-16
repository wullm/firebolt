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
int evolveMultipoles(struct multipoles *m, double tau_ini, double tau_fin,
                     double tol, double mass, double c_vel,
                     double (*redshift_func)(double log_tau),
                     double (*h_prime_func)(double k, double log_tau),
                     double (*eta_prime_func)(double k, double log_tau),
                     short verbose);
int convertMultipoleBasis_L2m(struct multipoles *mL, struct multipoles *mm, int l_max);
int convertMultipoleGauge_Nb(struct multipoles *mL,
                             double log_tau, double a, double mass,
                             double c_vel,
                             double (*delta_shift_func)(double k, double log_tau),
                             double (*theta_shift_func)(double k, double log_tau));
int resetMultipoles(struct multipoles *m);
int cleanMultipoles(struct multipoles *m);

int computeDensity(struct multipoles *m, double *density, double mass, double a,
                   double rho_nu, double factor_ncdm);


#endif
