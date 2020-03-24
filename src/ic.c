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

/* Fermi-Dirac distribution function */
double f0(double q) {
    double ksi = 0; //potential
    return 1.0/pow(2*M_PI,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));
}

// double dlnf0_dlnq(double q) {
//     return 1.0/pow(2*_PI_,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));
// }


/* Generate initial conditions for momentum q and wavenumber k at conformal
 * time tau */
void generate_ics(const struct background *bg, double q, double k, double tau,
                  double **Psi, int l_max) {

    /* Only Psi_0 through Psi_3 will have non-zero values at early times */
    /* Allocate and initialize all Psi's at 0. */
    *Psi = calloc(l_max + 1, sizeof(double));

    /* We assume adiabatic initial conditions. First, we need to calculate
     * some auxilliary quantities. */

    /* Initial adiabatic curvature perturbation (default 1) */
    double ini_curv = 1.0;

    /* Find a good row */
    int therow;
    for (int i=0; i<bg->nrow; i++) {
        if (bg->functions[1][i] > tau) {
            therow = i;
            break;
        }
    }

    /* Scan for the appropriate column titles */
    int id_rho_crit, id_rho_tot, id_H;
    int id_rho_g, id_rho_b, id_rho_cdm, id_rho_ur;
    int id_rho_idr, id_rho_idm_dr, id_rho_scf, id_rho_dr, id_rho_dcdm;
    int id_p_scf;
    char has_cdm = 0, has_ur = 0, has_idr = 0, has_idm_dr = 0, has_scf = 0;
    char has_dr = 0, has_dcdm = 0;
    /* Count number ncdm species among the columns */
    int n_ncdm = 0;
    for (int i=0; i<bg->ncol; i++) {
        char word[50];
        sscanf(bg->titles[i], "%s", word);

        if (strcmp(word, "H") == 0) {
            id_H = i;
        } else if (strcmp(word, "(.)rho_crit") == 0) {
            id_rho_crit = i;
        } else if (strcmp(word, "(.)rho_tot") == 0) {
            id_rho_tot = i;
        } else if (strcmp(word, "(.)rho_g") == 0) {
            id_rho_g = i;
        } else if (strcmp(word, "(.)rho_b") == 0) {
            id_rho_b = i;
        } else if (strcmp(word, "(.)rho_cdm") == 0) {
            has_cdm = 1;
            id_rho_cdm = i;
        } else if (strcmp(word, "(.)rho_ur") == 0) {
            has_ur = 1;
            id_rho_ur = i;
        } else if (strcmp(word, "(.)rho_idr") == 0) {
            has_idr = 1;
            id_rho_idr = i;
        } else if (strcmp(word, "(.)rho_idm_dr") == 0) {
            has_idm_dr = 1;
            id_rho_idm_dr = i;
        } else if (strcmp(word, "(.)rho_scf") == 0) {
            has_scf = 1;
            id_rho_scf = i;
        } else if (strcmp(word, "(.)rho_dr") == 0) {
            has_dr = 1;
            id_rho_dr = i;
        } else if (strcmp(word, "(.)rho_dcdm") == 0) {
            has_dcdm = 1;
            id_rho_dcdm = i;
        } else if (strcmp(word, "(.)p_scf") == 0) {
            has_scf = 1;
            id_p_scf = i;
        }

        /* Look for ncdm density columns */
        int ncdm_num;
        if (sscanf(bg->titles[i], "(.)rho_ncdm[%d]", &ncdm_num) > 0) {
            n_ncdm++;
        }
    }

    /* Scan for ncdm columns */
    int *id_rho_ncdm = malloc(n_ncdm * sizeof(int));
    int *id_p_ncdm = malloc(n_ncdm * sizeof(int));
    for (int i=0; i<bg->ncol; i++) {
        char word[50];
        sscanf(bg->titles[i], "%s", word);

        for (int j = 0; j<n_ncdm; j++) {
            char search_rho[50];
            char search_p[50];
            sprintf(search_rho, "(.)rho_ncdm[%d]", j);
            sprintf(search_p, "(.)p_ncdm[%d]", j);

            if (strcmp(word, search_rho) == 0) {
                id_rho_ncdm[j] = i;
            } else if (strcmp(word, search_p) == 0) {
                id_p_ncdm[j] = i;
            }
        }
    }

    /* Hubble constant today and at tau */
    double H0 = bg->functions[id_H][bg->nrow - 1];
    double H = bg->functions[id_H][therow];

    /* Scale factor and redshift */
    double z = bg->z[therow];
    double a = 1.0/(1+z);

    /* Densities */
    double rho_crit = bg->functions[id_rho_crit][therow];
    double rho_tot = bg->functions[id_rho_tot][therow];

    /* Curvature parameter */
    double K = a*a*(rho_tot - rho_crit);

    /* Calculate total matter and radiation densities */
    double rho_m = 0, rho_r = 0;

    /* Photons */
    rho_r += bg->functions[id_rho_g-1][therow];
    /* Baryons */
    rho_m += bg->functions[id_rho_b-1][therow];
    /* CDM */
    if (has_cdm) rho_m += bg->functions[id_rho_cdm-1][therow];
    /* Decaying dark matter, dark radiation */
    if (has_dcdm) rho_m += bg->functions[id_rho_dcdm-1][therow];
    if (has_dr) rho_r += bg->functions[id_rho_dr-1][therow];
    /* Scalar field (rho_r += 3P, rho_m += rho - 3P)*/
    if (has_scf) rho_r += 3 * bg->functions[id_p_scf-1][therow];
    if (has_scf) rho_m += bg->functions[id_rho_scf-1][therow] - 3 * bg->functions[id_p_scf-1][therow];
    /* Ultra-relativistic species */
    if (has_ur) rho_r += bg->functions[id_rho_ur-1][therow];
    /* Interacting dark matter and radiation */
    if (has_idm_dr) rho_m += bg->functions[id_rho_idm_dr-1][therow];
    if (has_idr) rho_r += bg->functions[id_rho_idr-1][therow];
    /* Massive neutrinos */
    for (int i=0; i<n_ncdm; i++) {
        /* Radiation contribution is 3P, density contribution is rho - 3P. */
        rho_r += 3 * bg->functions[id_p_ncdm[i]-1][therow];
        rho_m += bg->functions[id_rho_ncdm[i]-1][therow] - 3 * bg->functions[id_p_ncdm[i]-1][therow];
    }

    /* Matter and radiation density parameters */
    double Omega_m = rho_m / rho_crit;
    double Omega_r = rho_r / rho_crit;

    printf("Omega_m, Omega_r = %f %f\n", Omega_m, Omega_r);
    printf("We found %d ncdm species\n", n_ncdm);

    /* Parameter which enters into the exact solution for matter+radiation
     * cosmology (without Lambda). */
    double omega = a*Omega_m*H/sqrt(Omega_r);

    /* Parameter that is equal to 1 in a flat universe */
    double ss = 1.0 - 3*K/(k*k);

    /* Contribution of neutrinos and relativistic collisionless relics */
    double rho_nu = 0;
    if (has_dr) rho_nu += bg->functions[id_rho_dr-1][therow];
    if (has_ur) rho_nu += bg->functions[id_rho_ur-1][therow];
    if (has_idr) rho_nu += bg->functions[id_rho_idr-1][therow];
    /* Massive neutrinos */
    for (int i=0; i<n_ncdm; i++) {
        /* Radiation contribution is 3P, density contribution is rho - 3P. */
        rho_nu += bg->functions[id_rho_ncdm[i]-1][therow];
    }

    /* Fraction of relics in radiation density */
    double fnu = rho_nu/rho_r;

    printf("omega = %f, s^2 = %f, fnu = %f\n", omega, ss, fnu);

    /* More auxilliary quantities */
    double ktau = k*tau;

    /* Now, calculate the initial conditions for ultra-relativistic species */
    double delta_ur = - ktau*ktau/3 * (1.0 - omega*tau/5) * ini_curv * ss;
    double theta_ur = - k * ktau*ktau*ktau/36 / (4*fnu + 15) * (4*fnu + 11 + 12*ss
                      - 3*(8*fnu*fnu + 50*fnu + 275)/(2*fnu + 15)/20 * tau * omega)
                      * ini_curv * ss;
    double sigma_ur = ktau*ktau / (45 + 12*fnu) * (3*ss - 1)
                      * (1 + 0.25 * (4*fnu - 5)/(2*fnu+15) * tau * omega) * ini_curv;
    double l3_ur = ktau * ktau * 2. / (12*fnu + 45) / 7.;

    /* Now calculate the initial conditions for the neutrino species */
    double M_nu = 0.2;
    double eps = hypot(q, a * M_nu);

    /* Derivative of the distribution function */
    double h = 0.01;
    double df0_dq = 0, dlnf0_dlnq;

    df0_dq += (1./12.) * f0(q - 2*h);
    df0_dq -= (8./12.) * f0(q - 1*h);
    df0_dq += (8./12.) * f0(q + 1*h);
    df0_dq -= (1./12.) * f0(q + 2*h);
    df0_dq /= h;

    double f0_eval = f0(q);
    if (fabs(f0_eval) > 0) {
        dlnf0_dlnq = q/f0_eval * df0_dq;
    } else {
        dlnf0_dlnq = -q;
    }

    /* Calculate the Psi's */
    (*Psi)[0] = - 0.25 * delta_ur * dlnf0_dlnq;
    (*Psi)[1] = - eps / (3*q*k) * theta_ur * dlnf0_dlnq;
    (*Psi)[2] = - 0.5 * sigma_ur * dlnf0_dlnq;
    (*Psi)[3] = - 0.25 * l3_ur * dlnf0_dlnq;

    // **Psi + 1 = -1.;

    // printf("%e %e %e\n", f0_eval, df0_dq, dlnf0_dlnq);

    // printf("%f %f %f %f\n", delta_ur, theta_ur, sigma_ur, l3_ur);
}
