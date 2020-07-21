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

// #include "../include/firebolt_min.h"

#include <firebolt_nano.h>
#include "../include/input.h"
#include "../include/output.h"
#include "../include/perturb_data.h"
#include "../include/perturb_interp.h"

const char *fname;

static inline double sampleNorm() {
    double u = (double) rand()/RAND_MAX;
    double v = (double) rand()/RAND_MAX;

    double z0 = sqrt(-2 * log(u)) * cos(2 * M_PI * v);
    //double z1 = sqrt(-2 * log(u)) * sin(2 * M_PI * v);

    return z0;
}

static inline void tet(struct kernel *the_kernel) {
    double k = the_kernel->k;
    double logt = log(7.07);

    if (k == 0) {
        the_kernel->kern = 0.f;
    } else {
        // the_kernel->kern = perturbInterp(k, logt, 0);
        the_kernel->kern = - _Complex_I * the_kernel->kz * perturbInterp(k, logt, 0)/k/k;
    }
}

int writeGRF_H5(const double *box, int N, double boxlen, const char *fname) {
    /* Create the hdf5 file */
    hid_t h_file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create the Header group */
    hid_t h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Create dataspace for BoxSize attribute */
    const hsize_t arank = 1;
    const hsize_t adims[1] = {3}; //3D space
    hid_t h_aspace = H5Screate_simple(arank, adims, NULL);

    /* Create the BoxSize attribute and write the data */
    hid_t h_attr = H5Acreate1(h_grp, "BoxSize", H5T_NATIVE_DOUBLE, h_aspace, H5P_DEFAULT);
    double boxsize[3] = {boxlen, boxlen, boxlen};
    H5Awrite(h_attr, H5T_NATIVE_DOUBLE, boxsize);

    /* Close the attribute, corresponding dataspace, and the Header group */
    H5Aclose(h_attr);
    H5Sclose(h_aspace);
    H5Gclose(h_grp);

    /* Create the Field group */
    h_grp = H5Gcreate(h_file, "/Field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Create dataspace for the field */
    const hsize_t frank = 3;
    const hsize_t fdims[3] = {N, N, N}; //3D space
    hid_t h_fspace = H5Screate_simple(frank, fdims, NULL);

    /* Create the dataset for the field */
    hid_t h_data = H5Dcreate(h_grp, "Field", H5T_NATIVE_DOUBLE, h_fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Write the data */
    H5Dwrite(h_data, H5T_NATIVE_DOUBLE, h_fspace, h_fspace, H5P_DEFAULT, box);

    /* Close the dataset, corresponding dataspace, and the Field group */
    H5Dclose(h_data);
    H5Sclose(h_fspace);
    H5Gclose(h_grp);

    /* Close the file */
    H5Fclose(h_file);

    return 0;
}

static inline double first_perturb_index(double k, double log_tau) {
    return perturbInterp(k, log_tau, 0);
}
static inline double second_perturb_index(double k, double log_tau) {
    return perturbInterp(k, log_tau, 1);
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
    struct perturb_data ptdat;
    struct multipoles m; //multipoles is standard Legendre basis
    struct multipoles mmono; //multipoles in monomial basis
    struct grids grs;

    /* Read parameters & data */
    readParams(&pars, fname);
    readUnits(&us, fname);
    readCosmology(&cosmo, fname);
    readPerturb(&pars, &us, &ptdat);
    initPerturbInterp(&ptdat);

    /* Find indices corresponding to specific functions */
    int h_prime_index = 0, eta_prime_index = 0;
    int delta_shift_index = -1, theta_shift_index = -1;
    for (int i=0; i<ptdat.n_functions; i++) {
        if (strcmp(ptdat.titles[i], "h_prime") == 0) {
            h_prime_index = i;
            printf("Found '%s', index = %d\n", ptdat.titles[i], i);
        } else if (strcmp(ptdat.titles[i], "eta_prime") == 0) {
            eta_prime_index = i;
            printf("Found '%s', index = %d\n", ptdat.titles[i], i);
        } else if (strcmp(ptdat.titles[i], "t_cdm") == 0) {
            theta_shift_index = i;
            printf("Found '%s', index = %d\n", ptdat.titles[i], i);
        } else if (strcmp(ptdat.titles[i], "delta_shift_Nb_m") == 0) {
            delta_shift_index = i;
            printf("Found '%s', index = %d\n", ptdat.titles[i], i);
        }
    }

    /* Make sure that the interpolation functions correspond to h' and eta' */
    switchPerturbInterp(&ptdat, h_prime_index, 0);
    switchPerturbInterp(&ptdat, eta_prime_index, 1);

    printf("\n\n");

    pars.GridSize = N;
    pars.BoxLen = box_len;

    /* The system to solve */
    double k = pars.kSingle;
    double tau_ini = exp(ptdat.log_tau[0]);
    double tau_fin = pars.tauFinalSingle;
    double M = cosmo.M_nu * us.ElectronVolt / (cosmo.T_nu0 * us.kBoltzmann);
    double c_vel = us.SpeedOfLight;

    printf("Speed of light c = %f\n", c_vel);
    printf("Mass neutrino M = %f (%f)\n", M, cosmo.M_nu);

    /* Size of the problem */
    int l_max = pars.MaxMultipole;
    int l_max_convert = pars.MaxMultipoleConvert;
    double q_min = pars.MinMomentum;
    double q_max = pars.MaxMomentum;
    int q_steps = pars.NumberMomentumBins;
    double tol = pars.Tolerance;

    printf("\n");
    printf("[k] = [%f]\n", k);
    printf("[tau_ini, z_ini] = [%f, %.2e]\n", tau_ini, perturb_zAtLogTau(log(tau_ini)));
    printf("[tau_fin, z_fin] = [%f, %.2e]\n", tau_fin, perturb_zAtLogTau(log(tau_fin)));
    printf("[l_max, q_steps, q_max, tol] = [%d, %d, %.1f, %.3e]\n", l_max, q_steps, q_max, tol);
    printf("[l_max_convert] = %d\n", l_max_convert);
    printf("\n");

    printf("The initial time is %f\n", exp(ptdat.log_tau[0]));

    /* Initialize the multipoles */
    int k_size = 30;
    initMultipoles(&m, k_size, q_steps, l_max, q_min, q_max, k_min, k_max);

    /* Also initialize the multipoles in monomial basis (with much lower l_max) */
    initMultipoles(&mmono, k_size, q_steps, l_max_convert+1, q_min, q_max, k_min, k_max);

    /* Calculate the multipoles */
    evolveMultipoles(&m, tau_ini, tau_fin, tol, M, c_vel, perturb_zAtLogTau, first_perturb_index, second_perturb_index, pars.Verbose);

    /* Scale factor at final time */
    double z_fin = perturb_zAtLogTau(log(tau_fin));
    double a = 1./(1+z_fin);

    /* Select the transfer function needed for gauge transformation */
    switchPerturbInterp(&ptdat, delta_shift_index, 0);
    switchPerturbInterp(&ptdat, theta_shift_index, 1);

    /* N-body gauge transformation */
    convertMultipoleGauge_Nb(&m, log(tau_fin), a, M, c_vel, first_perturb_index, second_perturb_index);

    /* Switch back to h' and eta' mode */
    switchPerturbInterp(&ptdat, h_prime_index, 0);
    switchPerturbInterp(&ptdat, eta_prime_index, 1);

    /* Also convert to monomial basis */
    convertMultipoleBasis_L2m(&m, &mmono, l_max_convert);

    /* For each momentum bin */
    for (int i=0; i<q_steps; i++) {
        double q = mmono.q[i];
        /* For each wavenumber */
        for (int j=0; j<k_size; j++) {
            double k = mmono.k[j];

            double Psi0 = mmono.Psi[0 * q_steps * k_size + i * k_size + j];
            double Psi1 = mmono.Psi[1 * q_steps * k_size + i * k_size + j];
            double Psi2 = mmono.Psi[2 * q_steps * k_size + i * k_size + j];
            double Psi3 = mmono.Psi[3 * q_steps * k_size + i * k_size + j];

            printf("%f %f %e %e %e %e\n", q, k, Psi0, Psi1, Psi2, Psi3);
        }
    }

    printf("\n\n");

    /* For each momentum bin */
    for (int i=0; i<q_steps; i++) {
        double q = m.q[i];
        int j = 5;
        double k = m.k[j];

        double Psi0 = m.Psi[0 * q_steps * k_size + i * k_size + j];
        double Psi1 = m.Psi[1 * q_steps * k_size + i * k_size + j];
        double Psi2 = m.Psi[2 * q_steps * k_size + i * k_size + j];
        double Psi3 = m.Psi[3 * q_steps * k_size + i * k_size + j];
        double Psi4 = m.Psi[4 * q_steps * k_size + i * k_size + j];
        double Psi5 = m.Psi[5 * q_steps * k_size + i * k_size + j];
        double Psi6 = m.Psi[6 * q_steps * k_size + i * k_size + j];

        printf("%f %f %e %e %e %e %e %e %e\n", q, k, Psi0, Psi1, Psi2, Psi3, Psi4, Psi5, Psi6);
    }




    printf("Done with integrating. Processing the moments.\n");

    /* Initialize the multipole interpolation splines */
    initMultipoleInterp(&mmono);


    printf("Interpolation initialized.\n");

    /* Generate grids with the monomial multipoles */
    initGrids(pars.GridSize, pars.BoxLen, &mmono, &grs);
    generateGrids(&mmono, fbox, &grs);

    double eval_x = 16;
    double eval_y = 26;
    double eval_z = 14;
    double eval_nx1 = 0;
    double eval_ny1 = 1;
    double eval_nz1 = 0;
    double eval_nx2 = 0;
    double eval_ny2 = -1;
    double eval_nz2 = 0;

    /* Try evaluating */
    double e1,e2,f0_eval;
    for (int i=0; i<m.q_size; i++) {
        double q = m.q[i];
        f0_eval = f0(q);
        e1 = evalDensityBin(&grs, eval_x, eval_y, eval_z, eval_nx1, eval_ny1, eval_nz1, i);
        e2 = evalDensityBin(&grs, eval_x, eval_y, eval_z, eval_nx2, eval_ny2, eval_nz2, i);
        // e1 = evalDensity(&grs, 1./64.*256., 14./64.*256, 60./64.*256., 1., 0., 0., i);
        // e2 = evalDensity(&grs, 1./64.*256., 14./64.*256, 60./64.*256., -1., 0., 0., i);
        printf("%f %e %e %f %f\n", q, e1, e2, f0_eval*(1+e1), f0_eval*(1+e2));
    }


    printf("Stage two.\n");

    /* Try evaluating in-between bins */
    for (int i=0; i<10*m.q_size; i++) {
        double q = ((double) i )/10 * q_max/m.q_size;
        f0_eval = f0(q);
        e1 = evalDensity(&grs, m.q_size, log(q_min), log(q_max), eval_x, eval_y, eval_z, eval_nx1, eval_ny1, eval_nz1, q, 2);
        e2 = evalDensity(&grs, m.q_size, log(q_min), log(q_max), eval_x, eval_y, eval_z, eval_nx2, eval_ny2, eval_nz2, q, 2);
        printf("%f %e %e %f %f\n", q, e1, e2, f0_eval*(1+e1), f0_eval*(1+e2));
    }


    printf("Stage three.\n");


    /* Sample random points in phase space */
    for (int i=0; i<10000; i++) {
        double x = box_len * (double) rand() / RAND_MAX;
        double y = box_len * (double) rand() / RAND_MAX;
        double z = box_len * (double) rand() / RAND_MAX;

        /* Generate random point on the sphere */
        double nx = sampleNorm();
        double ny = sampleNorm();
        double nz = sampleNorm();
        double len = sqrt(nx*nx + ny*ny + nz*nz);

        nx /= len;
        ny /= len;
        nz /= len;

        double q = 10. * (double) rand() / RAND_MAX; //should be FD

        int bin = 1;

        double f0_eval = f0(q);
        double Psi0 = evalDensityBinZero(&grs, x, y, z, nx, ny, nz, bin);
        double Psi1 = evalDensityBinSimple(&grs, x, y, z, nx, ny, nz, bin);


        // evalDensity(&grs, m.q_size, log(q_min), log(q_max), x, y, z, nx, ny, nz, q);

        printf("%f %f %f\n", f0_eval, Psi0, Psi1);
    }

    printf("Bliep bloep.\n");

    /* Create FFT plan */
    fftw_plan c2r = fftw_plan_dft_c2r_3d(N, N, N, fbox, box, FFTW_ESTIMATE);

    /* Generate density grid */
    switchPerturbInterp(&ptdat, 7, 0);
    fft_apply_kernel(fbox, fbox, N, box_len, tet);
    printf("Stage four.\n");
    fft_execute(c2r);
    fft_normalize_c2r(box, N, box_len);


    char fff[50] = "density.hdf5";
    writeGRF_H5(box, N, box_len, fff);

    for (int i=0; i<m.q_size; i++) {
        double q = m.q[i];
        f0_eval = f0(q);
        e1 = gridCIC(box, N, box_len, eval_x, eval_y, eval_z);
        printf("%f\t %e\t %f\t %f\t %f\n", q, e1, 0., f0_eval*(1+e1), 0.);
    }

    fftw_destroy_plan(c2r);

    /* For each multipole/momentum bin pair, create the corresponding grid */
    for (int index_q=0; index_q<q_steps; index_q++) {
        for (int index_l=0; index_l<l_max_convert+1; index_l++) {
            char dbox_fname[DEFAULT_STRING_LENGTH];
            sprintf(dbox_fname, "grid_l%d_q%d.hdf5", index_l, index_q);

            double *destination = grs.grids + index_l * (N*N*N) * q_steps + index_q * (N*N*N);
            writeGRF_H5(destination, N, box_len, dbox_fname);
        }
    }

    /* Release the interpolation splines */
    cleanPerturbInterp(&ptdat);

    /* Clean up the remaining structures */
    cleanGrids(&grs);
    cleanMultipoles(&mmono);
    cleanMultipoles(&m);
    cleanMultipoleInterp();
    cleanPerturb(&ptdat);
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
