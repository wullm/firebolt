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

    /* Initialize interpolation splines */
    initPerturbInterp(&ptdat);
    bg_interp_init(&bg);

    /* The system to solve */
    // double q = 1;
    // double k = 0.5;
    double tau_ini = exp(ptdat.log_tau[0]);
    // // double final_tau = bg.functions[1][bg.nrow - 1];
    // double final_tau = 2150; //approx z=40
    double M = cosmo.M_nu * us.ElectronVolt / (cosmo.T_nu0 * us.kBoltzmann);

    /* Size of the problem */
    int l_max = pars.MaxMultipole;
    double q_max = pars.MaxMomentum;
    int q_steps = pars.NumberMomentumBins;
    double tol = pars.Tolerance;

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

    printf("\n");
    printf("[tau_ini, z] = [%f, %e]\n", tau_ini, bg_z_at_log_tau(log(tau_ini)));
    printf("[l_max, q_steps, q_max, tol] = [%d, %d, %.1f, %.3e]\n", l_max, q_steps, q_max, tol);
    printf("\n");

    /* The post-processed variables go in this structure */
    struct perturb_data hires_pt;
    /* Copy the size from the input perturbation */
    hires_pt.k_size = ptdat.k_size;
    hires_pt.tau_size = ptdat.tau_size;
    hires_pt.n_functions = ptdat.n_functions;
    /* Allocate memory */
    hires_pt.k = calloc(hires_pt.k_size, sizeof(double));
    hires_pt.log_tau = calloc(hires_pt.tau_size, sizeof(double));
    hires_pt.delta = calloc(hires_pt.n_functions * hires_pt.k_size * hires_pt.tau_size, sizeof(double));
    /* Copy the titles */
    hires_pt.titles = malloc(hires_pt.n_functions * sizeof(char*));
    for (int i=0; i<hires_pt.n_functions; i++) {
        hires_pt.titles[i] = malloc(strlen(ptdat.titles[i]) + 1);
        strcpy(hires_pt.titles[i], ptdat.titles[i]);
    }

    /* Copy the wavenumbers */
    for (int index_k=0; index_k<hires_pt.k_size; index_k++) {
        hires_pt.k[index_k] = ptdat.k[index_k];
    }
    /* Copy the log conformal times */
    for (int index_tau=0; index_tau<hires_pt.tau_size; index_tau++) {
        hires_pt.log_tau[index_tau] = ptdat.log_tau[index_tau];
    }

    double degeneracy = cosmo.Degeneracy;
    double factor_ncdm = degeneracy * 4 * M_PI * pow(cosmo.T_nu0 * us.kBoltzmann, 4) / pow(us.hPlanck/(2*M_PI), 3) / pow(us.SpeedOfLight, 7);

    double *ini_delta_nu = malloc(ptdat.k_size * sizeof(double));
    double *ini_theta_nu = malloc(ptdat.k_size * sizeof(double));
    double *ini_shear_nu = malloc(ptdat.k_size * sizeof(double));
    double *ini_l3_nu = malloc(ptdat.k_size * sizeof(double));

    /* Find initial conditions */
    for (int index_k=0; index_k<ptdat.k_size; index_k++) {
        /* Wavenumber */
        double k = ptdat.k[index_k];

        double log_tau = ptdat.log_tau[0];

        /* Compare with CLASS results */
        switchPerturbInterp(&ptdat, 4, 0);
        switchPerturbInterp(&ptdat, 5, 1);

        double class_delta_nu = perturbInterp(k, log_tau, 0);
        double class_theta_nu = perturbInterp(k, log_tau, 1);

        /* Also determine the shear and l3 */
        switchPerturbInterp(&ptdat, 6, 0); //shear
        double class_shear_nu = perturbInterp(k, log_tau, 0);
        switchPerturbInterp(&ptdat, 8, 1); //l3
        double class_l3_nu = perturbInterp(k, log_tau, 1);

        ini_delta_nu[index_k] = class_delta_nu;
        ini_theta_nu[index_k] = class_theta_nu;
        ini_shear_nu[index_k] = class_shear_nu;
        ini_l3_nu[index_k] = class_l3_nu;
    }


    /* Switch back to h' and eta' */
    switchPerturbInterp(&ptdat, 0, 0); // h'
    switchPerturbInterp(&ptdat, 1, 1); // eta'

    printf("Initial conditions done.\n");

    printf("The initial time is %f\n", exp(ptdat.log_tau[0]));


    #pragma omp parallel for
    for (int index_k=0; index_k<ptdat.k_size; index_k+=5) {
        /* Wavenumber */
        double k = ptdat.k[index_k];

        if (k < 0.02) continue; //skip the easy ones

        /* The variables that are integrated (just for this k) */
        double *rho_delta_nu = calloc(ptdat.tau_size, sizeof(double));
        double *rho_plus_p_theta_nu = calloc(ptdat.tau_size, sizeof(double));
        double *rho_plus_p_shear_nu = calloc(ptdat.tau_size, sizeof(double));
        double *rho_l3_nu = calloc(ptdat.tau_size, sizeof(double));
        double *delta_p_nu = calloc(ptdat.tau_size, sizeof(double));

        /* For each momentum bin */
        for (int j=0; j<q_steps; j++) {
            double dq = q_max/q_steps;
            double q = (j+0.5) * dq;

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

            Psi[0] = - 0.25 * ini_delta_nu[index_k] * dlnf0_dlnq;
            Psi[1] = - eps0 / (3*q*k) * ini_theta_nu[index_k] * dlnf0_dlnq;
            Psi[2] = - 0.5 * ini_shear_nu[index_k] * dlnf0_dlnq;
            Psi[3] = - 0.25 * ini_l3_nu[index_k] * dlnf0_dlnq;

            /* For each large timestep */
            for (int index_tau=0; index_tau<hires_pt.tau_size; index_tau++) {
                double logtau = ptdat.log_tau[index_tau];
                double tau = exp(logtau);

                double z = bg_z_at_log_tau(logtau);
                double a = 1./(1+z);
                double eps = hypot(q, a*M);
                double c_vel = us.SpeedOfLight;

                /* Integration only need after the first time step */
                if (index_tau > 0) {
                    double logtau_prev = ptdat.log_tau[index_tau-1];
                    double tau_prev = exp(logtau_prev);
                    evolve_gsl(&Psi, &ptdat, q, k, l_max, tau_prev, tau, M, c_vel, dlnf0_dlnq, tol);
                }

                /* Do the momentum integral */
                rho_delta_nu[index_tau] += q*q*eps*Psi[0]*f0_eval*dq;
                rho_plus_p_theta_nu[index_tau] += q*q*q*Psi[1]*f0_eval*dq;
                rho_plus_p_shear_nu[index_tau] += q*q*q*q/eps*Psi[2]*f0_eval*dq;
                rho_l3_nu[index_tau] += q*q*eps*Psi[3]*f0_eval*dq;
                delta_p_nu[index_tau] += q*q*q*q/eps*Psi[0]*f0_eval*dq;

                // printf("%d %f\n", index_tau, tau1);
            }

            // printf("bin %d\n", j);
            free(Psi);
        }

        /* After summing over all momentum bins, store the unprocessed results */
        int index_delta = 4;
        int index_theta = 5;
        int index_shear = 6;
        int index_l3 = 8;
        int index_cs2 = 7;
        for (int index_tau=0; index_tau<hires_pt.tau_size; index_tau++) {
            int Nk = hires_pt.k_size;
            int Nt = hires_pt.tau_size;
            hires_pt.delta[Nt * Nk * index_delta + Nk * index_tau + index_k] = rho_delta_nu[index_tau];
            hires_pt.delta[Nt * Nk * index_theta + Nk * index_tau + index_k] = rho_plus_p_theta_nu[index_tau];
            hires_pt.delta[Nt * Nk * index_shear + Nk * index_tau + index_k] = rho_plus_p_shear_nu[index_tau];
            hires_pt.delta[Nt * Nk * index_l3 + Nk * index_tau + index_k] = rho_l3_nu[index_tau];
            hires_pt.delta[Nt * Nk * index_cs2 + Nk * index_tau + index_k] = delta_p_nu[index_tau];
        }

        free(rho_delta_nu);
        free(rho_plus_p_theta_nu);
        free(rho_plus_p_shear_nu);
        free(rho_l3_nu);
        free(delta_p_nu);

        printf("\t\t%d %f\n", index_k, k);
    }

    printf("Done with integrating. Processing the moments.\n");

    /* Post-process the integrated moments */
    for (int index_tau=0; index_tau<hires_pt.tau_size; index_tau++) {
        double log_tau = hires_pt.log_tau[index_tau];
        double z = bg_z_at_log_tau(log_tau);
        double a = 1./(1+z);
        double factor = factor_ncdm / (a*a*a*a);

        /* Obtain the background density and pressure */
        bg_interp_switch_func(&bg, bti.id_rho_ncdm[0]);
        double rho_nu = bg_func_at_log_tau(log_tau);
        bg_interp_switch_func(&bg, bti.id_p_ncdm[0]);
        double p_nu = bg_func_at_log_tau(log_tau);

        for (int index_k=0; index_k<ptdat.k_size; index_k++) {
            double k = ptdat.k[index_k];

            int index_delta = 4;
            int index_theta = 5;
            int index_shear = 6;
            int index_l3 = 8;
            int index_cs2 = 7;
            /* Retrieve the intermediate products */
            int Nk = hires_pt.k_size;
            int Nt = hires_pt.tau_size;
            double rho_delta_nu = hires_pt.delta[Nt * Nk * index_delta + Nk * index_tau + index_k];
            double rho_plus_p_theta_nu = hires_pt.delta[Nt * Nk * index_theta + Nk * index_tau + index_k];
            double rho_plus_p_shear_nu = hires_pt.delta[Nt * Nk * index_shear + Nk * index_tau + index_k];
            double rho_l3_nu = hires_pt.delta[Nt * Nk * index_l3 + Nk * index_tau + index_k];
            double delta_p_nu = hires_pt.delta[Nt * Nk * index_cs2 + Nk * index_tau + index_k];

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

            /* Compare with CLASS results */
            switchPerturbInterp(&ptdat, 4, 0);
            switchPerturbInterp(&ptdat, 5, 1);

            double class_delta_nu = perturbInterp(k, log_tau, 0);
            double class_theta_nu = perturbInterp(k, log_tau, 1);

            /* Also determine the shear and l3 */
            switchPerturbInterp(&ptdat, 6, 0); //shear
            double class_shear_nu = perturbInterp(k, log_tau, 0);
            switchPerturbInterp(&ptdat, 8, 1); //l3
            double class_l3_nu = perturbInterp(k, log_tau, 1);

            /* Also determine the sound speed cs2 */
            switchPerturbInterp(&ptdat, 7, 1);
            double class_cs2_nu = perturbInterp(k, log_tau, 1);

            printf("rel. error [delta, theta, shear, l3, cs2] = [%f, %f, %f, %f, %f]\n", delta_nu/class_delta_nu, theta_nu/class_theta_nu, shear_nu/class_shear_nu, l3_nu/class_l3_nu, cs2_nu/class_cs2_nu);

            /* Store the finished results */
            hires_pt.delta[Nt * Nk * index_delta + Nk * index_tau + index_k] = delta_nu;
            hires_pt.delta[Nt * Nk * index_theta + Nk * index_tau + index_k] = theta_nu;
            hires_pt.delta[Nt * Nk * index_shear + Nk * index_tau + index_k] = shear_nu;
            hires_pt.delta[Nt * Nk * index_l3 + Nk * index_tau + index_k] = l3_nu;
            hires_pt.delta[Nt * Nk * index_cs2 + Nk * index_tau + index_k] = cs2_nu;
        }


        // printf("%d %f\n", index_k, k);
    }

    free(ini_delta_nu);
    free(ini_theta_nu);
    free(ini_shear_nu);
    free(ini_l3_nu);

    printf("All done!.\n");


    /* Write the processed perturbation to a file */
    write_perturb(&hires_pt, &pars, &us, pars.OutputFilename);



    // /* Switch back to metric splines */
    // switchPerturbInterp(&ptdat, 0, 0); // h'
    // switchPerturbInterp(&ptdat, 1, 1); // eta'
    // /* Check the accuracy */
    // for (int index_k=0; index_k<hires_pt.k_size; index_k++) {
    //     double k = hires_pt.k[index_k];
    //     for (int index_tau=9; index_tau<100; index_tau++) {
    //
    //         int index_delta = 4;
    //         int index_theta = 5;
    //         int index_shear = 6;
    //
    //         int id_delta = hires_pt.tau_size * hires_pt.k_size * index_delta
    //                      + hires_pt.k_size * index_tau + index_k;
    //         int id_theta = hires_pt.tau_size * hires_pt.k_size * index_theta
    //                      + hires_pt.k_size * index_tau + index_k;
    //         int id_shear = hires_pt.tau_size * hires_pt.k_size * index_shear
    //                      + hires_pt.k_size * index_tau + index_k;
    //
    //         double rel_delta = hires_pt.delta[id_delta] / ptdat.delta[id_delta];
    //         double rel_theta = hires_pt.delta[id_theta] / ptdat.delta[id_theta] / (cosmo.h * cosmo.h);
    //         double rel_shear = hires_pt.delta[id_shear] / ptdat.delta[id_shear];
    //
    //         printf("%f %f %f %f\n", k, rel_delta, rel_theta, rel_shear);
    //     }
    // }
    //
    // k = 0.000006;
    //
    // final_tau = exp(ptdat.log_tau[1]);
    //
    // double z = bg_z_at_log_tau(log(final_tau));
    // double a = 1./(1+z);
    // double factor = factor_ncdm / (a*a*a*a);
    //
    //
    // switchPerturbInterp(&ptdat, 0, 0);
    // switchPerturbInterp(&ptdat, 1, 1);
    //
    // double rho_delta_nu = 0;
    // double rho_plus_p_theta_nu = 0;
    // double rho_plus_p_shear_nu = 0;
    //
    // for (int j=0; j<q_steps; j++) {
    //     double *Psi;
    //
    //     double dq = q_max/q_steps;
    //     double q = (j+0.5) * dq;
    //     double eps = hypot(q, a*M);
    //
    //     /* Trapezoid rule */
    //     if (j == q_steps - 1) {
    //         dq *= 0.5;
    //     }
    //
    //     /* Generate initial conditions */
    //     generate_ics(&bg, &bti, q, k, tau_ini, &Psi, l_max);
    //
    //     /* Print the initial conditions */
    //     // printf("\n\n");
    //     // printf("%f %e %e %e %e %e\n", tau, Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
    //
    //     /* Derivative of the distribution function (5-point stencil) */
    //     double y = 0.0001;
    //     double dlnf0_dlnq = compute_dlnf0_dlnq(q, y);
    //     double f0_eval = f0(q);
    //
    //     evolve_gsl(&Psi, &ptdat, &bg, q, k, l_max, tau_ini, final_tau, M, dlnf0_dlnq);
    //     // printf("%f %e %e %e %e %e\n", q, Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
    //     printf("%f %e %e %e %e %e\n", q, q*q*eps*Psi[0]*f0_eval*dq, q*q*q*Psi[1]*f0_eval*dq, q*q*q*q/eps*Psi[2]*f0_eval*dq, Psi[3], Psi[4]);
    //
    //     rho_delta_nu += q*q*eps*Psi[0]*f0_eval*dq;
    //     rho_plus_p_theta_nu += q*q*q*Psi[1]*f0_eval*dq;
    //     rho_plus_p_shear_nu += q*q*q*q/eps*Psi[2]*f0_eval*dq;
    //
    //     /* Free the integrated variables */
    //     free(Psi);
    // }
    //
    // rho_delta_nu *= factor;
    // rho_plus_p_theta_nu *= k*factor;
    // rho_plus_p_shear_nu *= 2*factor/3;
    //
    // bg_interp_switch_func(&bg, bti.id_rho_ncdm[1]);
    // double rho_nu = bg_func_at_log_tau(log(final_tau));
    // bg_interp_switch_func(&bg, bti.id_p_ncdm[1]);
    // double p_nu = bg_func_at_log_tau(log(final_tau));
    // double w_nu = p_nu/rho_nu;
    //
    // double delta_nu = rho_delta_nu/rho_nu;
    // double theta_nu = rho_plus_p_theta_nu/(rho_nu + p_nu);
    // double shear_nu = rho_plus_p_shear_nu/(rho_nu + p_nu);
    //
    //
    // /* Transformations to N-body gauge */
    // switchPerturbInterp(&ptdat, 0, 0); // h'
    // switchPerturbInterp(&ptdat, 1, 1); // eta'
    // double h_prime = perturbInterp(k, log(final_tau), 0);
    // double eta_prime = perturbInterp(k, log(final_tau), 1);
    // double alpha = (h_prime + 6*eta_prime)/(2*k*k);
    //
    // switchPerturbInterp(&ptdat, 2, 0); // H_T_Nb_prime
    // switchPerturbInterp(&ptdat, 3, 1); // t_tot
    //
    // double H_T_Nb_prime = perturbInterp(k, log(final_tau), 0);
    // double theta_shift = H_T_Nb_prime + alpha*k*k;
    // double theta_tot = perturbInterp(k, log(final_tau), 1)  - theta_shift;
    // /* End transformation variables */
    //
    // /* Little h correction for theta's */
    // theta_shift *= (cosmo.h * cosmo.h);
    // theta_tot *= (cosmo.h * cosmo.h);
    //
    // /* Determine a'/a = H * a */
    // bg_interp_switch_func(&bg, bti.id_H);
    // double H = bg_func_at_log_tau(log(final_tau));
    // double H_conf = H*a;
    //
    // printf("\n\n");
    // printf("alpha \t\t= %e\n", alpha);
    // printf("alpha*k^2 \t= %e\n", alpha*k*k);
    // printf("theta_tot \t= %e (%e)\n", theta_tot + theta_shift, theta_tot);
    // printf("H_T_Nb_prime \t= %e\n", H_T_Nb_prime);
    // printf("theta_shift \t= %e\n", theta_shift);
    // printf("delta_shift \t= %e\n", 3*(1+w_nu)*H_conf*theta_tot/k/k);
    // printf("w \t\t= %f\n", w_nu);
    // printf("H \t\t= %f\n\n\n", H_conf);
    //
    // /* Synchronous to N-body gauge transformation */
    // delta_nu += 3*(1+w_nu)*H_conf*theta_tot/k/k;
    // theta_nu += theta_shift;
    //
    // /* Compare with CLASS results */
    // switchPerturbInterp(&ptdat, 4, 0);
    // switchPerturbInterp(&ptdat, 5, 1);
    //
    // double class_delta_nu = perturbInterp(k, log(final_tau), 0);
    // double class_theta_nu = perturbInterp(k, log(final_tau), 1);
    //
    // /* Little h correction for theta's */
    // class_theta_nu *= cosmo.h * cosmo.h;
    //
    // /* Also determine the shear */
    // switchPerturbInterp(&ptdat, 6, 1);
    // double class_shear_nu = perturbInterp(k, log(final_tau), 1);
    //
    // printf("%f %e %e %e\n", k, delta_nu, theta_nu, shear_nu);
    // printf("%f %e %e %e\n\n", k, class_delta_nu, class_theta_nu, class_shear_nu);
    // printf("%f %f %f %f\n", k, class_delta_nu/delta_nu, class_theta_nu/theta_nu, class_shear_nu/shear_nu);


    /* Release the interpolation splines */
    bg_interp_free(&bg);
    cleanPerturbInterp(&ptdat);

    /* Clean up the remaining structures */
    cleanPerturb(&ptdat);
    cleanPerturb(&hires_pt);
    cleanBackgroundTitles(&bti);
    cleanBackground(&bg);
    cleanParams(&pars);

    /* Timer */
    gettimeofday(&stop, NULL);
    long unsigned microsec = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
    printf("Time elapsed: %.3f ms\n", microsec/1000.);
}
