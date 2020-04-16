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
#include <sys/time.h>

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

    /* Timer */
    struct timeval stop, start;
    gettimeofday(&start, NULL);

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

    /* Comppute derivatives */
    // computeDerivs(&ptdat);

    /* Initialize interpolation splines */
    rend_interp_init(&ptdat);
    bg_interp_init(&bg);

    /* The system to solve */
    double k = pars.kSingle;
    double tau_ini = exp(ptdat.log_tau[0]);
    double tau_fin = pars.tauFinalSingle;
    double M = cosmo.M_nu * us.ElectronVolt / (cosmo.T_nu0 * us.kBoltzmann);

    /* Size of the problem */
    int l_max = pars.MaxMultipole;
    double q_max = pars.MaxMomentum;
    int q_steps = pars.NumberMomentumBins;
    double tol = pars.Tolerance;

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

    printf("\n");
    printf("[k] = [%f]\n", k);
    printf("[tau_ini, z_ini] = [%f, %.2e]\n", tau_ini, bg_z_at_log_tau(log(tau_ini)));
    printf("[tau_fin, z_fin] = [%f, %.2e]\n", tau_fin, bg_z_at_log_tau(log(tau_fin)));
    printf("[l_max, q_steps, q_max, tol] = [%d, %d, %.1f, %.3e]\n", l_max, q_steps, q_max, tol);
    printf("\n");

    double degeneracy = cosmo.Degeneracy;
    double factor_ncdm = degeneracy * 4 * M_PI * pow(cosmo.T_nu0 * us.kBoltzmann, 4) / pow(us.hPlanck/(2*M_PI), 3) / pow(us.SpeedOfLight, 7);

    /* Find initial conditions */
    double ini_delta_nu;
    double ini_theta_nu;
    double ini_shear_nu;
    double ini_l3_nu;

    {
        double log_tau = log(tau_ini);

        /* Redshift */
        // double z = bg_z_at_log_tau(log_tau);
        // double a = 1./(1+z);

        /* Obtain the background density and pressure */
        // bg_interp_switch_func(&bg, bti.id_rho_ncdm[0]);
        // double rho_nu = bg_func_at_log_tau(log_tau);
        // bg_interp_switch_func(&bg, bti.id_p_ncdm[0]);
        // double p_nu = bg_func_at_log_tau(log_tau);
        // double w_nu = p_nu/rho_nu; //equation of state

        // /* Determine a'/a = H * a */
        // bg_interp_switch_func(&bg, bti.id_H);
        // double H = bg_func_at_log_tau(log_tau);
        // double H_conf = H*a; // = a'/a

        // /* Transformations to N-body gauge */
        // rend_interp_switch_source(&ptdat, 0, 0); // h'
        // rend_interp_switch_source(&ptdat, 1, 1); // eta'
        // double h_prime = rend_interp(k, log_tau, 0);
        // double eta_prime = rend_interp(k, log_tau, 1);
        // double alpha = (h_prime + 6*eta_prime)/(2*k*k);
        //
        // rend_interp_switch_source(&ptdat, 2, 0); // H_T_Nb_prime
        // rend_interp_switch_source(&ptdat, 3, 1); // t_tot
        //
        // double H_T_Nb_prime = rend_interp(k, log_tau, 0);
        // double theta_shift = H_T_Nb_prime + alpha*k*k;
        // double theta_tot = rend_interp(k, log_tau, 1) - theta_shift;
        //
        // /* Little h correction for theta's */
        // theta_shift *= (cosmo.h * cosmo.h);
        // theta_tot *= (cosmo.h * cosmo.h);

        /* Synchronous to N-body gauge transformation */
        // delta_nu += 3*(1+w_nu)*H_conf*theta_tot/k/k;
        // theta_nu += theta_shift;

        /* Compare with CLASS results */
        rend_interp_switch_source(&ptdat, 4, 0);
        rend_interp_switch_source(&ptdat, 5, 1);

        double class_delta_nu = rend_interp(k, log_tau, 0);
        double class_theta_nu = rend_interp(k, log_tau, 1);

        /* Little h correction for theta's */
        // class_theta_nu *= cosmo.h * cosmo.h;

        // class_delta_nu -= 3*(1+w_nu)*H_conf*theta_tot/k/k;
        // class_theta_nu -= theta_shift;

        /* Also determine the shear and l3 */
        rend_interp_switch_source(&ptdat, 6, 0); //shear
        double class_shear_nu = rend_interp(k, log_tau, 0);
        rend_interp_switch_source(&ptdat, 8, 1); //l3
        double class_l3_nu = rend_interp(k, log_tau, 1);

        ini_delta_nu = class_delta_nu;
        ini_theta_nu = class_theta_nu;
        ini_shear_nu = class_shear_nu;
        ini_l3_nu = class_l3_nu;
    }

    /* Switch back to h' and eta' */
    rend_interp_switch_source(&ptdat, 0, 0); // h'
    rend_interp_switch_source(&ptdat, 1, 1); // eta'

    printf("Initial conditions done.\n");

    printf("The initial time is %f\n", exp(ptdat.log_tau[0]));

    /* The variables that are integrated (just for this k) */
    double rho_delta_nu = 0;
    double rho_plus_p_theta_nu = 0;
    double rho_plus_p_shear_nu = 0;
    double rho_l3_nu = 0;
    double delta_p_nu = 0;

    double q_min = 0.01;

    /* For each momentum bin */
    for (int j=0; j<q_steps; j++) {
        double dlogq = (log(q_max) - log(q_min))/q_steps;
        double q = q_min * exp((j+0.5) * dlogq);
        double dq = exp((j+0.5) * dlogq - 1);
        // double q = (j+0.5) * dq;

        /* Derivative of the distribution function (5-point stencil) */
        double y = 0.0001;
        double dlnf0_dlnq = compute_dlnf0_dlnq(q, y);
        double f0_eval = f0(q);

        /* The neutrino multipoles */
        double *Psi = calloc(l_max+1,sizeof(double));

        /* Generate initial conditions */
        // generate_ics(&bg, &bti, q, k, exp(ptdat.log_tau[0]), &Psi, l_max);

        /* At the initial time */
        double z0 = bg_z_at_log_tau(ptdat.log_tau[0]);
        double a0 = 1./(1+z0);
        double eps0 = hypot(q, a0*M);

        Psi[0] = - 0.25 * ini_delta_nu * dlnf0_dlnq;
        Psi[1] = - eps0 / (3*q*k) * ini_theta_nu * dlnf0_dlnq;
        Psi[2] = - 0.5 * ini_shear_nu * dlnf0_dlnq;
        Psi[3] = - 0.25 * ini_l3_nu * dlnf0_dlnq;

        double log_tau_fin = log(tau_fin);

        double z = bg_z_at_log_tau(log_tau_fin);
        double a = 1./(1+z);
        double eps = hypot(q, a*M);

        evolve_gsl(&Psi, &ptdat, &bg, q, k, l_max, tau_ini, tau_fin, M, dlnf0_dlnq, tol);

        /* Do the momentum integrals */
        rho_delta_nu += q*q*eps*Psi[0]*f0_eval*dq;
        rho_plus_p_theta_nu += q*q*q*Psi[1]*f0_eval*dq;
        rho_plus_p_shear_nu += q*q*q*q/eps*Psi[2]*f0_eval*dq;
        rho_l3_nu += q*q*eps*Psi[3]*f0_eval*dq;
        delta_p_nu += q*q*q*q/eps*Psi[0]*f0_eval*dq;

        printf("%f %e %e %e %e\n", q, Psi[0], Psi[1], Psi[2], Psi[3]);
        // printf("%f %e %e %e %e\n", q, Psi[0]*f0_eval, Psi[1]*f0_eval, Psi[2]*f0_eval, Psi[3]*f0_eval);
        // printf("%f %e %e %e %e\n", q, q*q*eps*Psi[0]*f0_eval*dq, q*q*q*Psi[1]*f0_eval*dq, q*q*q*q/eps*Psi[2]*f0_eval*dq, q*q*eps*Psi[3]*f0_eval*dq);

        free(Psi);
    }

    printf("Done with integrating. Processing the moments.\n");

    /* Post-process the integrated moments */
    double log_tau_fin = log(tau_fin);
    double z = bg_z_at_log_tau(log_tau_fin);
    double a = 1./(1+z);
    double factor = factor_ncdm / (a*a*a*a);

    /* Obtain the background density and pressure */
    bg_interp_switch_func(&bg, bti.id_rho_ncdm[0]);
    double rho_nu = bg_func_at_log_tau(log_tau_fin);
    bg_interp_switch_func(&bg, bti.id_p_ncdm[0]);
    double p_nu = bg_func_at_log_tau(log_tau_fin);
    double w_nu = p_nu/rho_nu; //equation of state

    printf("\n");
    printf("w = %f\n",w_nu);
    printf("\n");

    /* Determine a'/a = H * a */
    // bg_interp_switch_func(&bg, bti.id_H);
    // double H = bg_func_at_log_tau(log_tau_fin);
    // double H_conf = H*a; // = a'/a

    rho_delta_nu *= factor;
    rho_plus_p_theta_nu *= k*factor;
    rho_plus_p_shear_nu *= 2*factor/3;
    rho_l3_nu *= factor;
    delta_p_nu *= factor/3;

    double delta_nu = rho_delta_nu/rho_nu;
    double theta_nu = rho_plus_p_theta_nu/(rho_nu + p_nu);
    double shear_nu = rho_plus_p_shear_nu/(rho_nu + p_nu);
    double l3_nu = rho_l3_nu/rho_nu;
    double cs2_nu = delta_p_nu/rho_delta_nu / (-k*k);

    // /* Transformations to N-body gauge */
    // rend_interp_switch_source(&ptdat, 0, 0); // h'
    // rend_interp_switch_source(&ptdat, 1, 1); // eta'
    // double h_prime = rend_interp(k, log_tau_fin, 0);
    // double eta_prime = rend_interp(k, log_tau_fin, 1);
    // double alpha = (h_prime + 6*eta_prime)/(2*k*k);
    //
    // rend_interp_switch_source(&ptdat, 2, 0); // H_T_Nb_prime
    // rend_interp_switch_source(&ptdat, 3, 1); // t_tot
    //
    // double H_T_Nb_prime = rend_interp(k, log_tau_fin, 0);
    // double theta_shift = H_T_Nb_prime + alpha*k*k;
    // double theta_tot = rend_interp(k, log_tau_fin, 1) - theta_shift;
    //
    // /* Little h correction for theta's */
    // theta_shift *= (cosmo.h * cosmo.h);
    // theta_tot *= (cosmo.h * cosmo.h);

    /* Synchronous to N-body gauge transformation */
    // delta_nu += 3*(1+w_nu)*H_conf*theta_tot/k/k;
    // theta_nu += theta_shift;

    /* Compare with CLASS results */
    rend_interp_switch_source(&ptdat, 4, 0);
    rend_interp_switch_source(&ptdat, 5, 1);

    double class_delta_nu = rend_interp(k, log_tau_fin, 0);
    double class_theta_nu = rend_interp(k, log_tau_fin, 1);

    /* Little h correction for theta's */
    // class_theta_nu *= cosmo.h * cosmo.h;

    // class_delta_nu -= 3*(1+w_nu)*H_conf*theta_tot/k/k;
    // class_theta_nu -= theta_shift;

    /* Also determine the shear and l3*/
    rend_interp_switch_source(&ptdat, 6, 0);
    double class_shear_nu = rend_interp(k, log_tau_fin, 0);
    rend_interp_switch_source(&ptdat, 8, 1);
    double class_l3_nu = rend_interp(k, log_tau_fin, 1);

    /* Also determine the sound speed cs2 */
    rend_interp_switch_source(&ptdat, 7, 1);
    double class_cs2_nu = rend_interp(k, log_tau_fin, 1);

    printf("rel. error [delta, theta, shear, l3, cs2] = [%f, %f, %f, %f, %f]\n", delta_nu/class_delta_nu, theta_nu/class_theta_nu, shear_nu/class_shear_nu, l3_nu/class_l3_nu, cs2_nu/class_cs2_nu);
    printf("values [delta, theta, shear, l3, cs2] = [%e, %e, %e, %e, %e]\n", delta_nu, theta_nu, shear_nu, l3_nu, cs2_nu);

    printf("All done!.\n");

    rend_interp_switch_source(&ptdat, 6, 0);
    printf("\n\n %e %e\n", rend_interp(k, ptdat.log_tau[0], 1), rend_interp(k, log_tau_fin, 1));
    printf("%e %e\n", ini_shear_nu, class_shear_nu);

    // printf("\n\n");
    // for (int i=0; i<ptdat.tau_size; i++) {
    //     double shear = rend_interp(k, ptdat.log_tau[i], 1);
    //     printf("%e %e\n", exp(ptdat.log_tau[i]), shear);
    // }

    /* Release the interpolation splines */
    bg_interp_free(&bg);
    rend_interp_free(&ptdat);

    /* Clean up the remaining structures */
    cleanPerturb(&ptdat);
    cleanBackgroundTitles(&bti);
    cleanBackground(&bg);
    cleanParams(&pars);

    /* Timer */
    gettimeofday(&stop, NULL);
    long unsigned microsec = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
    printf("Time elapsed: %.3f ms\n", microsec/1000.);
}
