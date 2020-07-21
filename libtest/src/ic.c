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

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "../include/ic.h"
#include "../../include/multipoles.h"




/* Generate initial conditions for momentum q and wavenumber k at conformal
 * time tau */
 void generate_ics(const struct background *bg, const struct background_title_ids *bti,
                   double q, double k, double tau, double **Psi, int l_max) {

    /* Only Psi_0 through Psi_3 will have non-zero values at early times */
    /* Allocate and initialize all Psi's at 0. */
    // *Psi = calloc(l_max + 1, sizeof(double));

    /* We assume adiabatic initial conditions. First, we need to calculate
     * some auxilliary quantities. */

    /* Initial adiabatic curvature perturbation (default 1) */
    double ini_curv = 1.0;

    /* Find a good row */
    int therow = 0;
    // for (int i=0; i<bg->nrow; i++) {
    //     if (bg->functions[1][i] > tau) {
    //         therow = i;
    //         break;
    //     }
    // }

    /* Hubble constant today and at tau */
    double H0 = bg->functions[bti->id_H][bg->nrow - 1];
    double H = bg->functions[bti->id_H][therow];

    /* Scale factor and redshift */
    double z = bg->z[therow];
    double a = 1.0/(1+z);

    /* Densities */
    double rho_crit = bg->functions[bti->id_rho_crit][therow];
    // double rho_tot = bg->functions[bti->id_rho_tot][therow];

    /* Densities today */
    double rho_crit0 = bg->functions[bti->id_rho_crit][bg->nrow - 1];
    double rho_tot0 = bg->functions[bti->id_rho_tot][bg->nrow - 1];

    /* Curvature parameter */
    double Omega_k0 = (rho_crit0 - rho_tot0)/rho_crit0;
    double K = -a*a*Omega_k0*H0*H0;

    /* Calculate total matter and radiation densities */
    double rho_m = 0, rho_r = 0;

    /* Photons */
    rho_r += bg->functions[bti->id_rho_g][therow];
    /* Baryons */
    rho_m += bg->functions[bti->id_rho_b][therow];
    /* CDM */
    if (bti->has_cdm) rho_m += bg->functions[bti->id_rho_cdm][therow];
    /* Decaying dark matter, dark radiation */
    if (bti->has_dcdm) rho_m += bg->functions[bti->id_rho_dcdm][therow];
    if (bti->has_dr) rho_r += bg->functions[bti->id_rho_dr][therow];
    /* Scalar field (rho_r += 3P, rho_m += rho - 3P)*/
    if (bti->has_scf) rho_r += 3 * bg->functions[bti->id_p_scf][therow];
    if (bti->has_scf) rho_m += bg->functions[bti->id_rho_scf][therow] - 3 * bg->functions[bti->id_p_scf][therow];
    /* Ultra-relativistic species */
    if (bti->has_ur) rho_r += bg->functions[bti->id_rho_ur][therow];
    /* Interacting dark matter and radiation */
    if (bti->has_idm_dr) rho_m += bg->functions[bti->id_rho_idm_dr][therow];
    if (bti->has_idr) rho_r += bg->functions[bti->id_rho_idr][therow];
    /* Massive neutrinos */
    for (int i=0; i<bti->n_ncdm; i++) {
        /* Radiation contribution is 3P, density contribution is rho - 3P. */
        rho_r += 3 * bg->functions[bti->id_p_ncdm[i]][therow];
        rho_m += bg->functions[bti->id_rho_ncdm[i]][therow] - 3 * bg->functions[bti->id_p_ncdm[i]][therow];
    }

    /* Matter and radiation density parameters */
    double Omega_m = rho_m / rho_crit;
    double Omega_r = rho_r / rho_crit;

    /* Parameter which enters into the exact solution for matter+radiation
     * cosmology (without Lambda). */
    double omega = a*Omega_m*H/sqrt(Omega_r);


    /* Parameter that is equal to 1 in a flat universe */
    double ss = 1.0 - 3*K/(k*k);

    /* Contribution of neutrinos and relativistic collisionless relics */
    double rho_nu = 0;
    if (bti->has_dr) rho_nu += bg->functions[bti->id_rho_dr][therow];
    if (bti->has_ur) rho_nu += bg->functions[bti->id_rho_ur][therow];
    if (bti->has_idr) rho_nu += bg->functions[bti->id_rho_idr][therow];
    /* Massive neutrinos */
    for (int i=0; i<bti->n_ncdm; i++) {
        /* Radiation contribution is 3P, density contribution is rho - 3P. */
        rho_nu += bg->functions[bti->id_rho_ncdm[i]][therow];
    }

    /* Fraction of relics in radiation density */
    double fnu = rho_nu/rho_r;

    // printf("Omega_m, Omega_r = %f %f\n", Omega_m, Omega_r);
    // printf("We found %d ncdm species\n", n_ncdm);
    // printf("z = %e, \t a = %e, \t H = %e\n\n", z, a, H0);
    // printf("omega = %e, s^2 = %e, fnu = %e\n", omega, ss, fnu);

    /* More auxilliary quantities */
    double ktau = k*tau;

    /* Now, calculate the initial conditions for ultra-relativistic species */
    double delta_ur = - ktau*ktau/3 * (1.0 - omega*tau/5) * ini_curv * ss;
    double theta_ur = - k * ktau*ktau*ktau/36 / (4*fnu + 15) * (4*fnu + 11 + 12*ss
                      - 3*(8*fnu*fnu + 50*fnu + 275)/(2*fnu + 15)/20 * tau * omega)
                      * ini_curv * ss;
    double sigma_ur = ktau*ktau / (45 + 12*fnu) * (3*ss - 1)
                      * (1 + 0.25 * (4*fnu - 5)/(2*fnu+15) * tau * omega) * ini_curv;
    double l3_ur = ktau * ktau * ktau * 2. / (12*fnu + 45) / 7.;

    // printf("%e %e %e %e\n", delta_ur, theta_ur, sigma_ur, l3_ur);

    /* Now calculate the initial conditions for the neutrino species */
    double M_nu = 0.2;
    double eps = hypot(q, a * M_nu);

    /* Derivative of the distribution function */
    double h = 0.01;
    double dlnf0_dlnq = compute_dlnf0_dlnq(q, h);

    /* Calculate the Psi's */
    (*Psi)[0] = - 0.25 * delta_ur * dlnf0_dlnq;
    (*Psi)[1] = - eps / (3*q*k) * theta_ur * dlnf0_dlnq;
    (*Psi)[2] = - 0.5 * sigma_ur * dlnf0_dlnq;
    (*Psi)[3] = - 0.25 * l3_ur * dlnf0_dlnq;
}
