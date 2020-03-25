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
    struct perturb_data ptdat;

    /* Read parameters & data */
    readParams(&pars, fname);
    readUnits(&us, fname);
    readCosmology(&cosmo, fname);
    readPerturb(&pars, &us, &ptdat);
    readBackground(&pars, &us, &cosmo, &bg);

    /* Initialize interpolation splines */
    rend_interp_init(&ptdat);
    bg_interp_init(&bg);

    /* The system to solve */
    double q = 1;
    double k = 1;
    double tau = 0.05;
    // double final_tau = bg.functions[1][bg.nrow - 1];
    double final_tau = 2150; //approx z=40
    double M = 1190; //mass in units of neutrino temperature today (0.2 eV)
    printf("[q, k, tau, z] = [%f, %f, %f, %e]\n", q, k, tau, bg_z_at_tau(log(tau)));

    // final_tau = 100;

    /* Size of the problem */
    int l_max = 2000;
    double *Psi;

    // for (int i=0; i<100; i++) {
    //     double y = (i+1) * 0.01;
    //     printf("%f %e\n", y, rend_interp(y, log(100), 1));
    // }

    /* Generate initial conditions */
    generate_ics(&bg, q, k, tau, &Psi, l_max);

    /* Print the initial conditions */
    printf("\n\n");
    printf("%f %e %e %e %e %e\n", tau, Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

    /* Derivative of the distribution function (5-point stencil) */
    double h = 0.0001;
    double dlnf0_dlnq = compute_dlnf0_dlnq(q, h);

    evolve_gsl(&Psi, &ptdat, &bg, q, k, l_max, tau, final_tau, M, dlnf0_dlnq);
    printf("%f %e %e %e %e %e\n", final_tau, Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

    /* Release the interpolation splines */
    bg_interp_free(&bg);
    rend_interp_free(&ptdat);

    /* Free the integrated variables */
    free(Psi);

    /* Clean up the remaining structures */
    cleanPerturb(&ptdat);
    cleanBackground(&bg);
    cleanParams(&pars);
}
