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

#include <complex.h>
#include <fftw3.h>

#include "../include/firebolt.h"

const char *fname;

static inline unsigned long long binomial(unsigned long n, unsigned long k) {
    unsigned long long c = 1, i;
    if (k > n-k) // take advantage of symmetry
        k = n-k;

    for (i = 1; i <= k; i++, n--) {
        if (c/i > UINT_MAX/n) return 0;
        c = c / i * n + c % i * n / i;  // split c * n / i into (c / i * i + c % i) * n / i
    }

    return c;
}

double fbinomial(double n, double k) {
    return exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1));
}

static inline void kernel_transfer_function(struct kernel *the_kernel) {
    double k = the_kernel->k;
    double kern = (k > 0) ? rend_interp(k, log(2150), 0) : 0.f;
    the_kernel->kern = kern;
}

int main(int argc, char *argv[]) {
    printf("OK\n");


    /* Read the gaussian random field */
    int N;
    double box_len;
    double *box;
    readGRF_H5(&box, &N, &box_len, "gaussian_pure.hdf5");

    /* Determine the maximum and minimum wavenumbers */
    double dk = 2*M_PI/box_len;
    double k_max = sqrt(3)*dk*N/2;
    double k_min = dk;

    /* Ensure a safe error margin */
    k_max *= 1.2;
    k_min /= 1.2;

    printf("[kmin, kmax] = [%f, %f]\n", k_min, k_max);

    printf("Test test %f %f %f\n", fbinomial(1,2), fbinomial(2.5,1), fbinomial(10.5,5));

    /* Fourier transform */
    fftw_complex *fbox = (fftw_complex*) fftw_malloc(N*N*(N/2+1)*sizeof(fftw_complex));
    fftw_plan r2c = fftw_plan_dft_r2c_3d(N, N, N, box, fbox, FFTW_ESTIMATE);
    fft_execute(r2c);
    fft_normalize_r2c(fbox,N,box_len);

    /* Read options */
    const char *fname = "default.ini";
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
    struct multipoles m;
    struct grids grs;

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

    // /* Apply transfer function */
    // rend_interp_switch_source(&ptdat, 4, 0);
    // fft_apply_kernel(fbox, fbox, N, box_len, kernel_transfer_function);
    //
    // /* Export the real box */
    // fft_c2r_export(fbox, N, box_len, "mooi.hdf5");


    pars.GridSize = N;
    pars.BoxLen = box_len;

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

    printf("\n");
    printf("[k] = [%f]\n", k);
    printf("[tau_ini, z_ini] = [%f, %.2e]\n", tau_ini, bg_z_at_log_tau(log(tau_ini)));
    printf("[tau_fin, z_fin] = [%f, %.2e]\n", tau_fin, bg_z_at_log_tau(log(tau_fin)));
    printf("[l_max, q_steps, q_max, tol] = [%d, %d, %.1f, %.3e]\n", l_max, q_steps, q_max, tol);
    printf("\n");

    printf("The initial time is %f\n", exp(ptdat.log_tau[0]));

    double q_min = 0.01;




    /* Initialize the multipoles */
    int k_size = 30;
    int l_size = 12;
    initMultipoles(&m, k_size, q_steps, l_size, q_min, q_max, k_min, k_max);

    /* Calculate the multipoles */

    /* For each momentum bin */
    for (int i=0; i<m.q_size; i++) {
        double q = m.q[i];

        /* Derivative of the distribution function (5-point stencil) */
        double y = 0.0001;
        double dlnf0_dlnq = compute_dlnf0_dlnq(q, y);
        double f0_eval = f0(q);

        /* For each wavenumber */
        for (int j=0; j<m.k_size; j++) {
            double k = m.k[j];

            /* The neutrino multipoles */
            double *Psi = calloc(l_max+1,sizeof(double));

            /* Retrieve the stored values from the multipole array */
            for (int l=0; l<l_size; l++) {
                Psi[l] = m.Psi[l * q_steps * k_size + i * k_size + j];
            }

            evolve_gsl(&Psi, &ptdat, &bg, q, k, l_max, tau_ini, tau_fin, M, dlnf0_dlnq, tol);

            printf("%f %f %e %e %e %e %e\n", q, k, Psi[0], Psi[1], Psi[2], Psi[3], f0_eval);

            /* Store the result */
            // for (int l=0; l<l_size; l++) {
            //     m.Psi[l * q_steps * k_size + i * k_size + j] = Psi[l];
            // }

            // /* Store the result */
            // for (int l=0; l<15; l++) {
            //     /* Expand the Legendre polynomial */
            //     for (int n=0; n<=ceil(l/2); n++) {
            //         double b1 = binomial(l, n);
            //         double b2 = binomial(2*(l-n), l);
            //         double factor = pow(0.5,l) * pow(-1, n) * b1 * b2 * (2*l + 1) * pow(-1, l);
            //
            //         int this = (l - 2*n);
            //         if (this < l_size)
            //         m.Psi[this * q_steps * k_size + i * k_size + j] += factor * Psi[l] / pow(k, this);
            //     }
            // }

            /* Store the result */
            for (int n=0; n<l_size; n++) {
                /* Expand the Legendre polynomial */
                for (int l=0; l<=n; l++) {
                    double b1 = binomial(n, l);
                    double b2 = fbinomial((n+l-1)/2., n);
                    double factor = pow(2,n) * b1 * b2 * (2*n + 1) * pow(-1, n);

                    if (l==12) {
                        printf("%e\n", factor * Psi[n] / pow(k, l));
                    }

                    if (l < l_size) {
                        m.Psi[l * q_steps * k_size + i * k_size + j] += factor * Psi[n] / pow(k, l);
                    }
                }
            }

            free(Psi);
        }
    }

    // /* For each momentum bin */
    // for (int j=0; j<q_steps; j++) {
    //     double dlogq = (log(q_max) - log(q_min))/q_steps;
    //     double q = q_min * exp((j+0.5) * dlogq);
    //
    //     /* Derivative of the distribution function (5-point stencil) */
    //     double y = 0.0001;
    //     double dlnf0_dlnq = compute_dlnf0_dlnq(q, y);
    //     double f0_eval = f0(q);
    //
    //     /* The neutrino multipoles */
    //     double *Psi = calloc(l_max+1,sizeof(double));
    //
    //     /* Initial conditions are irrelevant and washed out by sources */
    //     Psi[0] = 0;
    //     Psi[1] = 0;
    //     Psi[2] = 0;
    //     Psi[3] = 0;
    //
    //     evolve_gsl(&Psi, &ptdat, &bg, q, k, l_max, tau_ini, tau_fin, M, dlnf0_dlnq, tol);
    //
    //     printf("%f %e %e %e %e %e\n", q, Psi[0], Psi[1], Psi[2], Psi[3], f0_eval);
    //     // printf("%f %e %e %e %e\n", q, Psi[0]*f0_eval, Psi[1]*f0_eval, Psi[2]*f0_eval, Psi[3]*f0_eval);
    //
    //     free(Psi);
    // }
    //
    printf("Done with integrating. Processing the moments.\n");

    /* For each wavenumber */
    for (int j=0; j<m.k_size; j++) {
        double k = m.k[j];
        printf("You what mae? %f\n", k);
    }


    /* Initialize the multipole interpolation splines */
    initMultipoleInterp(&m);

    /* Try it out */
    double kk = 0.6;
    double Psi = multipoleInterp(kk);
    printf("Interpolation gave %f %e\n", kk, Psi);

    /* Switch to a different multipole */
    switchMultipoleInterp(&m, 1, 1);
    Psi = multipoleInterp(kk);
    printf("Interpolation gave %f %e\n", kk, Psi);

    /* Switch to a different multipole */
    switchMultipoleInterp(&m, 2, 1);
    Psi = multipoleInterp(kk);
    printf("Interpolation gave %f %e\n", kk, Psi);


    initGrids(&pars, &m, &grs);
    generateGrids(&pars, &us, &cosmo, &m, fbox, &grs);

    /* Try evaluating */
    double e1,e2,f0_eval;
    for (int i=0; i<m.q_size; i++) {
        double q = m.q[i];
        f0_eval = f0(q);
        e1 = evalDensity(&grs, 1./64.*256., 14./64.*256, 60./64.*256., 1., 0., 0., i);
        e2 = evalDensity(&grs, 1./64.*256., 14./64.*256, 60./64.*256., -1., 0., 0., i);
        printf("%f\t %f\t %f\t %f\t %f\n", q, e1, e2, f0_eval*(1+e1), f0_eval*(1+e2));
    }

    /* Release the interpolation splines */
    bg_interp_free(&bg);
    rend_interp_free(&ptdat);

    /* Clean up the remaining structures */
    cleanGrids(&grs);
    cleanMultipoles(&m);
    cleanMultipoleInterp(&m);
    cleanPerturb(&ptdat);
    cleanBackgroundTitles(&bti);
    cleanBackground(&bg);
    cleanParams(&pars);

    /* Free the GRF */
    free(box);
    fftw_free(fbox);
    fftw_destroy_plan(r2c);

    /* Timer */
    gettimeofday(&stop, NULL);
    long unsigned microsec = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
    printf("Time elapsed: %.3f ms\n", microsec/1000.);
}
