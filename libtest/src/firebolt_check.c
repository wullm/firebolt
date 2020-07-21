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


static inline double first_perturb_index(double k, double log_tau) {
    return perturbInterp(k, log_tau, 0);
}
static inline double second_perturb_index(double k, double log_tau) {
    return perturbInterp(k, log_tau, 1);
}

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
    struct cosmology_params cosmo;
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
    initPerturbInterp(&ptdat);
    bg_interp_init(&bg);

    /* The system to solve */
    double k = pars.kSingle;
    double tau_ini = exp(ptdat.log_tau[0]);
    double tau_fin = pars.tauFinalSingle;
    double M = cosmo.M_nu * us.ElectronVolt / (cosmo.T_nu0 * us.kBoltzmann);
    double c = us.SpeedOfLight;

    /* Size of the problem */
    int l_max = pars.MaxMultipole;
    double q_min = pars.MinMomentum;
    double q_max = pars.MaxMomentum;
    int q_steps = pars.NumberMomentumBins;
    double tol = pars.Tolerance;

    double chi = 0.;
    int steps = 1000;
    double dtau = (tau_fin - tau_ini)/steps;
    /* For each momentum bin */
    for (int j=0; j<q_steps; j++) {
        double dq = q_max/q_steps;
        double q = (j+0.5) * dq;
        double f0_eval = 1./(exp(q)+1) / log(2.);

        for (int i=0; i<steps; i++) {
            double tau = tau_ini + dtau * (i + 0.5);
            double log_tau = log(tau);
            double z = bg_z_at_log_tau(log_tau);
            double a = 1./(1+z);
            double v = q/sqrt(M*M*a*a + q*q/c/c);
            chi += v*dtau*dq*f0_eval;
        }
    }

    printf("chi = %e (M = %f)\n", chi, M);


    printf("\n");
    printf("[k] = [%f]\n", k);
    printf("[tau_ini, z_ini] = [%f, %.2e]\n", tau_ini, bg_z_at_log_tau(log(tau_ini)));
    printf("[tau_fin, z_fin] = [%f, %.2e]\n", tau_fin, bg_z_at_log_tau(log(tau_fin)));
    printf("[l_max, q_steps, q_max, tol] = [%d, %d, %.1f, %.3e]\n", l_max, q_steps, q_max, tol);
    printf("\n");

    printf("c=%f\n", c);


    // exit(0);

    /* Find indices corresponding to specific functions */
    int h_prime_index = 0, eta_prime_index = 0;
    int d_ncdm_index = 0, t_ncdm_index = 0 , shear_ncdm_index = 0,
        l3_ncdm_index = 0, cs2_ncdm_index = 0;
    int delta_shift_index = -1, theta_shift_index = -1;
    for (int i=0; i<ptdat.n_functions; i++) {
        if (strcmp(ptdat.titles[i], "h_prime") == 0) {
            h_prime_index = i;
            printf("Found 'h_prime', index = %d\n", h_prime_index);
        } else if (strcmp(ptdat.titles[i], "eta_prime") == 0) {
            eta_prime_index = i;
            printf("Found 'eta_prime', index = %d\n", eta_prime_index);
        } else if (strcmp(ptdat.titles[i], "d_ncdm[0]") == 0) {
            d_ncdm_index = i;
            printf("Found 'd_ncdm[0]', index = %d\n", d_ncdm_index);
        } else if (strcmp(ptdat.titles[i], "t_ncdm[0]") == 0) {
            t_ncdm_index = i;
            printf("Found 't_ncdm[0]', index = %d\n", t_ncdm_index);
        } else if (strcmp(ptdat.titles[i], "shear_ncdm[0]") == 0) {
            shear_ncdm_index = i;
            printf("Found 'shear_ncdm[0]', index = %d\n", shear_ncdm_index);
        } else if (strcmp(ptdat.titles[i], "l3_ncdm[0]") == 0) {
            l3_ncdm_index = i;
            printf("Found 'l3_ncdm[0]', index = %d\n", l3_ncdm_index);
        } else if (strcmp(ptdat.titles[i], "cs2_ncdm[0]") == 0) {
            cs2_ncdm_index = i;
            printf("Found 'cs2_ncdm[0]', index = %d\n", cs2_ncdm_index);
        } else if (strcmp(ptdat.titles[i], "t_cdm") == 0) {
            theta_shift_index = i;
            printf("Found '%s', index = %d\n", ptdat.titles[i], i);
        } else if (strcmp(ptdat.titles[i], "delta_shift_Nb_m") == 0) {
            delta_shift_index = i;
            printf("Found '%s', index = %d\n", ptdat.titles[i], i);
        }
    }




    // switchPerturbInterp(&ptdat, 4, 0);
    // switchPerturbInterp(&ptdat, 5, 1);
    //
    // for (int i=0; i<100; i++) {
    //     double y = exp(log(1e-5) + (log(10) - log(1e-5)) * (i+1) * 0.01);
    //     double h_prime = perturbInterp(y, log(80.), 0);
    //     double eta_prime = perturbInterp(y, log(80.), 1);
    //     // printf("%f %e %e\n", y, h_prime/6, h_prime/15 + 2*eta_prime/5);
    //     printf("%f %e %e\n", y, h_prime, eta_prime);
    // }
    //
    // printf("So far so good\t z=%f\n", bg_z_at_log_tau(log(80.)));



    double degeneracy = cosmo.Degeneracy;
    double factor_ncdm = degeneracy * 4 * M_PI * pow(cosmo.T_nu0 * us.kBoltzmann, 4) / pow(us.hPlanck/(2*M_PI), 3) / pow(us.SpeedOfLight, 7);

    /* Pre-compute the final results expected from CLASS */
    double log_tau_fin = log(tau_fin);
    // double z_fin = bg_z_at_log_tau(log_tau_fin);
    // double a_fin = 1./(1+z_fin);

    switchPerturbInterp(&ptdat, d_ncdm_index, 0);
    switchPerturbInterp(&ptdat, t_ncdm_index, 1);
    double class_delta_nu = perturbInterp(k, log_tau_fin, 0);
    double class_theta_nu = perturbInterp(k, log_tau_fin, 1);
    /* Also determine the shear and l3*/
    switchPerturbInterp(&ptdat, shear_ncdm_index, 0);
    double class_shear_nu = perturbInterp(k, log_tau_fin, 0);
    switchPerturbInterp(&ptdat, l3_ncdm_index, 1);
    double class_l3_nu = perturbInterp(k, log_tau_fin, 1);
    /* Also determine the sound speed cs2 */
    switchPerturbInterp(&ptdat, cs2_ncdm_index, 1);
    double class_cs2_nu = perturbInterp(k, log_tau_fin, 1);


    /* Find initial conditions */
    double ini_delta_nu;
    double ini_theta_nu;
    double ini_shear_nu;
    double ini_l3_nu;

    {
        double log_tau = log(tau_ini);

        /* Compare with CLASS results */
        switchPerturbInterp(&ptdat, d_ncdm_index, 0);
        switchPerturbInterp(&ptdat, t_ncdm_index, 1);

        double class_delta_nu = perturbInterp(k, log_tau, 0);
        double class_theta_nu = perturbInterp(k, log_tau, 1);

        printf("class_delta = %e, class_theta = %e\n", class_delta_nu, class_theta_nu);

        /* Also determine the shear and l3 */
        switchPerturbInterp(&ptdat, shear_ncdm_index, 0); //shear
        double class_shear_nu = perturbInterp(k, log_tau, 0);
        switchPerturbInterp(&ptdat, l3_ncdm_index, 1); //l3
        double class_l3_nu = perturbInterp(k, log_tau, 1);

        ini_delta_nu = class_delta_nu;
        ini_theta_nu = class_theta_nu;
        ini_shear_nu = class_shear_nu;
        ini_l3_nu = class_l3_nu;
    }

    /* Switch back to h' and eta' */
    switchPerturbInterp(&ptdat, h_prime_index, 0); // h'
    switchPerturbInterp(&ptdat, eta_prime_index, 1); // eta'

    printf("Initial conditions done.\n");

    printf("The initial time is %f\n", exp(ptdat.log_tau[0]));

    /* The variables that are integrated (just for this k) */
    double rho_delta_nu = 0;
    double rho_plus_p_theta_nu = 0;
    double rho_plus_p_shear_nu = 0;
    double rho_l3_nu = 0;
    double delta_p_nu = 0;

    /* The neutrino multipoles */
    double *Psi = calloc(l_max+1,sizeof(double));
    double dlogq = (log(q_max) - log(q_min)) / q_steps;

    /* For each momentum bin */
    for (int j=0; j<q_steps; j++) {
        double q = q_min * exp((j + 0.5) * dlogq);
        double dq = dlogq * q;

        /* Derivative of the distribution function (5-point stencil) */
        double y = 0.0001;
        double dlnf0_dlnq = compute_dlnf0_dlnq(q, y);
        double f0_eval = f0(q);

        /* Generate initial conditions */
        // generate_ics(&bg, &bti, q, k, exp(ptdat.log_tau[0]), &Psi, l_max);

        /* Reset the Psi vector */
        memset(Psi, 0, (l_max+1)*sizeof(double));

        /* At the initial time */
        double z0 = bg_z_at_log_tau(ptdat.log_tau[0]);
        double a0 = 1./(1+z0);
        double eps0 = hypot(q, a0*M);

        Psi[0] = - 0.25 * ini_delta_nu * dlnf0_dlnq;
        Psi[1] = - eps0 / (3*q*k) * ini_theta_nu * dlnf0_dlnq;
        Psi[2] = - 0.5 * ini_shear_nu * dlnf0_dlnq;
        Psi[3] = - 0.25 * ini_l3_nu * dlnf0_dlnq;

        double z = bg_z_at_log_tau(log_tau_fin);
        double a = 1./(1+z);
        double eps = hypot(q, a*M);
        double c_vel = us.SpeedOfLight;

        evolve_gsl(&Psi, q, k, l_max, tau_ini, tau_fin, M, c_vel, dlnf0_dlnq, bg_z_at_log_tau, first_perturb_index, second_perturb_index, tol);

        /* Do the momentum integrals */
        rho_delta_nu += q*q*eps*Psi[0]*f0_eval*dq;
        rho_plus_p_theta_nu += q*q*q*Psi[1]*f0_eval*dq;
        rho_plus_p_shear_nu += q*q*q*q/eps*Psi[2]*f0_eval*dq;
        rho_l3_nu += q*q*eps*Psi[3]*f0_eval*dq;
        delta_p_nu += q*q*q*q/eps*Psi[0]*f0_eval*dq;

        // printf("%f %e %e %e %e\n", q, Psi0, Psi1, Psi2, Psi3);
        printf("%f %e %e %e %e\n", q, Psi[0], Psi[1], Psi[2], Psi[3]);
    }

    printf("Done with integrating. Processing the moments.\n");


    /* Post-process the integrated moments */
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
    printf("w = %f\n", w_nu);
    printf("rho_nu = %e\n", rho_nu);
    printf("rho_delta_nu = %e\n", rho_delta_nu);
    printf("factor = %e\n", factor_ncdm);
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

    /* Select the transfer function needed for gauge transformation */
    switchPerturbInterp(&ptdat, delta_shift_index, 0);
    switchPerturbInterp(&ptdat, theta_shift_index, 1);

    /* N-body gauge transformation */
    delta_nu += perturbInterp(k, log_tau_fin, 0);
    theta_nu += perturbInterp(k, log_tau_fin, 1);

    /* Switch back to h' and eta' mode */
    switchPerturbInterp(&ptdat, h_prime_index, 0);
    switchPerturbInterp(&ptdat, eta_prime_index, 1);

    printf("rel. to CLASS [delta, theta, shear, l3, cs2] = [%f, %f, %f, %f, %f]\n", delta_nu/class_delta_nu, theta_nu/class_theta_nu, shear_nu/class_shear_nu, l3_nu/class_l3_nu, cs2_nu/class_cs2_nu);
    printf("values [delta, theta, shear, l3, cs2] = [%e, %e, %e, %e, %e]\n", delta_nu, theta_nu, shear_nu, l3_nu, cs2_nu);

    free(Psi);

    /* Release the interpolation splines */
    bg_interp_free(&bg);
    cleanPerturbInterp(&ptdat);

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
