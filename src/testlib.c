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

#include "../include/firebolt_min.h"
#include "../include/background.h"
#include "../include/background_interp.h"

const char *fname;

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
    struct cosmology_params cosmo;
    struct background bg;
    struct background_title_ids bti;
    struct perturb_data ptdat;
    struct multipoles m; //multipoles is standard Legendre basis
    struct multipoles mmono; //multipoles in monomial basis
    struct grids grs;

    /* Read parameters & data */
    readParams(&pars, fname);
    readUnits(&us, fname);
    readCosmology(&cosmo, fname);
    readPerturb(&pars, &us, &ptdat);
    rend_interp_init(&ptdat);

    /* Read background cosmology */
    readBackground(&pars, &us, &cosmo, &bg);
    parseBackgroundTitles(&bg, &bti);
    bg_interp_init(&bg);

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
    initMultipoles(&m, k_size, q_steps, l_max, q_min, q_max, k_min, k_max);

    /* Also initialize the multipoles in monomial basis (with much lower l_max) */
    int l_size = 12;
    initMultipoles(&mmono, k_size, q_steps, l_size, q_min, q_max, k_min, k_max);

    /* Calculate the multipoles */
    evolveMultipoles(&m, &ptdat, tau_ini, tau_fin, tol, M, bg_z_at_log_tau);

    /* Also convert to monomial basis */
    convertMultipoleBasis_L2m(&m, &mmono, l_size-1);

    printf("Done with integrating. Processing the moments.\n");

    /* Initialize the multipole interpolation splines */
    initMultipoleInterp(&mmono);

    /* Generate grids with the monomial multipoles */
    initGrids(&pars, &mmono, &grs);
    generateGrids(&pars, &us, &mmono, fbox, &grs);

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
    cleanMultipoles(&mmono);
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
