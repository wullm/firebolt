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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/firebolt.h"

const char *fname;

int main(int argc, char *argv[]) {
    if (argc == 1) {
        printf("No parameter file specified.\n");
        return 0;
    }

    /* Read options */
    const char *fname = argv[1];
    printf("The parameter file is %s\n", fname);

    /* Firebolt structures */
    struct params pars;
    struct units us;
    struct cosmology cosmo;
    struct background bg;
    struct background_title_ids bti;
    struct perturb_data ptdat;

    /* Read parameters & data */
    readParams(&pars, fname);
    readUnits(&us, fname);
    readCosmology(&cosmo, fname);
    readPerturb(&pars, &us, &ptdat);
    readBackground(&pars, &us, &cosmo, &bg);
    parseBackgroundTitles(&bg, &bti);

    /* Initialize interpolation splines */
    rend_interp_init(&ptdat);
    bg_interp_init(&bg);

    /* The system to solve */
    // double q = 1;
    double k = 0.5;
    double tau_ini = 0.05;
    // double final_tau = bg.functions[1][bg.nrow - 1];
    double final_tau = 2150; //approx z=40
    double M = 1189.133740; //mass in units of neutrino temperature today (0.2 eV)
    printf("[k, tau_ini, z] = [%f, %f, %e]\n", k, tau_ini, bg_z_at_log_tau(log(tau_ini)));

    /* Size of the problem */
    int l_max = pars.MaxMultipole;

    // rend_interp_switch_source(&ptdat, 4, 0);
    // rend_interp_switch_source(&ptdat, 5, 1);
    //
    // for (int i=0; i<100; i++) {
    //     double y = exp(log(1e-5) + (log(10) - log(1e-5)) * (i+1) * 0.01);
    //     double h_prime = rend_interp(y, log(80.), 0);
    //     double eta_prime = rend_interp(y, log(80.), 1);
    //     // printf("%f %e %e\n", y, h_prime/6, h_prime/15 + 2*eta_prime/5);
    //     printf("%f %e %e\n", y, h_prime, eta_prime);
    // }
    //
    // printf("So far so good\t z=%f\n", bg_z_at_log_tau(log(80.)));

    double q_max = pars.MaxMomentum;
    int q_steps = pars.NumberMomentumBins;

    printf("[l_max, q_steps, q_max] = [%d, %d, %f]\n", l_max, q_steps, q_max);
    printf("\n");

    double z = bg_z_at_log_tau(log(final_tau));
    double a = 1./(1+z);

    double degeneracy = 1.;
    double factor_ncdm = degeneracy * 4 * M_PI * pow(cosmo.T_nu0 * us.kBoltzmann, 4) / pow(us.hPlanck/(2*M_PI), 3) / pow(us.SpeedOfLight, 7);
    double factor = factor_ncdm / (a*a*a*a);

    rend_interp_switch_source(&ptdat, 0, 0);
    rend_interp_switch_source(&ptdat, 1, 1);

    double rho_delta_nu = 0;
    double rho_plus_p_theta_nu = 0;
    double rho_plus_p_shear_nu = 0;

    for (int j=0; j<q_steps; j++) {
        double *Psi;

        double dq = q_max/q_steps;
        double q = (j+0.5) * dq;
        double eps = hypot(q, a*M);

        /* Trapezoid rule */
        if (j == q_steps - 1) {
            dq *= 0.5;
        }

        /* Generate initial conditions */
        generate_ics(&bg, &bti, q, k, tau_ini, &Psi, l_max);

        /* Print the initial conditions */
        // printf("\n\n");
        // printf("%f %e %e %e %e %e\n", tau, Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

        /* Derivative of the distribution function (5-point stencil) */
        double y = 0.0001;
        double dlnf0_dlnq = compute_dlnf0_dlnq(q, y);
        double f0_eval = f0(q);

        evolve_gsl(&Psi, &ptdat, &bg, q, k, l_max, tau_ini, final_tau, M, dlnf0_dlnq);
        // printf("%f %e %e %e %e %e\n", q, Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
        printf("%f %e %e %e %e %e\n", q, q*q*eps*Psi[0]*f0_eval*dq, q*q*q*Psi[1]*f0_eval*dq, q*q*q*q/eps*Psi[2]*f0_eval*dq, Psi[3], Psi[4]);

        rho_delta_nu += q*q*eps*Psi[0]*f0_eval*dq;
        rho_plus_p_theta_nu += q*q*q*Psi[1]*f0_eval*dq;
        rho_plus_p_shear_nu += q*q*q*q/eps*Psi[2]*f0_eval*dq;

        /* Free the integrated variables */
        free(Psi);
    }

    rho_delta_nu *= factor;
    rho_plus_p_theta_nu *= k*factor;
    rho_plus_p_shear_nu *= 2*factor/3;

    bg_interp_switch_func(&bg, bti.id_rho_ncdm[1]);
    double rho_nu = bg_func_at_log_tau(log(final_tau));
    bg_interp_switch_func(&bg, bti.id_p_ncdm[1]);
    double p_nu = bg_func_at_log_tau(log(final_tau));
    double w_nu = p_nu/rho_nu;

    double delta_nu = rho_delta_nu/rho_nu;
    double theta_nu = rho_plus_p_theta_nu/(rho_nu + p_nu);
    double shear_nu = rho_plus_p_shear_nu/(rho_nu + p_nu);


    /* Transformations to N-body gauge */
    rend_interp_switch_source(&ptdat, 0, 0); // h'
    rend_interp_switch_source(&ptdat, 1, 1); // eta'
    double h_prime = rend_interp(k, log(final_tau), 0);
    double eta_prime = rend_interp(k, log(final_tau), 1);
    double alpha = (h_prime + 6*eta_prime)/(2*k*k);

    rend_interp_switch_source(&ptdat, 2, 0); // H_T_Nb_prime
    rend_interp_switch_source(&ptdat, 3, 1); // t_tot

    double H_T_Nb_prime = rend_interp(k, log(final_tau), 0);
    double theta_shift = H_T_Nb_prime + alpha*k*k;
    double theta_tot = rend_interp(k, log(final_tau), 1)  - theta_shift;
    /* End transformation variables */

    /* Little h correction for theta's */
    theta_shift *= (cosmo.h * cosmo.h);
    theta_tot *= (cosmo.h * cosmo.h);

    /* Determine a'/a = H * a */
    bg_interp_switch_func(&bg, bti.id_H);
    double H = bg_func_at_log_tau(log(final_tau));
    double H_conf = H*a;

    printf("\n\n");
    printf("alpha \t\t= %e\n", alpha);
    printf("alpha*k^2 \t= %e\n", alpha*k*k);
    printf("theta_tot \t= %e (%e)\n", theta_tot + theta_shift, theta_tot);
    printf("H_T_Nb_prime \t= %e\n", H_T_Nb_prime);
    printf("theta_shift \t= %e\n", theta_shift);
    printf("delta_shift \t= %e\n", 3*(1+w_nu)*H_conf*theta_tot/k/k);
    printf("w \t\t= %f\n", w_nu);
    printf("H \t\t= %f\n\n\n", H_conf);

    /* Synchronous to N-body gauge transformation */
    delta_nu += 3*(1+w_nu)*H_conf*theta_tot/k/k;
    theta_nu += theta_shift;

    /* Compare with CLASS results */
    rend_interp_switch_source(&ptdat, 4, 0);
    rend_interp_switch_source(&ptdat, 5, 1);

    double class_delta_nu = rend_interp(k, log(final_tau), 0);
    double class_theta_nu = rend_interp(k, log(final_tau), 1);

    /* Little h correction for theta's */
    class_theta_nu *= cosmo.h * cosmo.h;

    /* Also determine the shear */
    rend_interp_switch_source(&ptdat, 6, 1);
    double class_shear_nu = rend_interp(k, log(final_tau), 1);

    printf("%f %e %e %e\n", k, delta_nu, theta_nu, shear_nu);
    printf("%f %e %e %e\n\n", k, class_delta_nu, class_theta_nu, class_shear_nu);
    printf("%f %f %f %f\n", k, class_delta_nu/delta_nu, class_theta_nu/theta_nu, class_shear_nu/shear_nu);

    /* Release the interpolation splines */
    bg_interp_free(&bg);
    rend_interp_free(&ptdat);

    /* Clean up the remaining structures */
    cleanPerturb(&ptdat);
    cleanBackgroundTitles(&bti);
    cleanBackground(&bg);
    cleanParams(&pars);
}
