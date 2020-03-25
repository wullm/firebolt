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

    struct params pars;
    struct units us;
    struct cosmology cosmo;
    struct background bg;
    struct perturb_data ptdat;

    readParams(&pars, fname);
    readUnits(&us, fname);
    readCosmology(&cosmo, fname);
    readBackground(&pars, &us, &cosmo, &bg);
    readPerturb(&pars, &us, &ptdat);

    double q = 1;
    double k = 10;
    double tau = 0.05;
    printf("[q, k, tau] = [%f, %f, %f]\n", q, k, tau);

    int l_max = 1000;
    double *Psi;

    generate_ics(&bg, q, k, tau, &Psi, l_max);


    printf("\n\n");
    printf("%f %e %e %e %e %e\n", tau, Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

    rend_interp_init(&ptdat);
    rend_interp_switch_source(&ptdat, 5);
    bg_interp_init(&bg);

    double final_tau = bg.functions[1][bg.nrow - 1];
    double M = 1190; //mass in units of neutrino temperature today (0.2 eV)

    /* Derivative of the distribution function (5-point stencil) */
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

    for (int i=0; i<20; i++) {
        tau *= pow(10, 0.25);
        final_tau = tau*pow(10, 0.25);
        evolve(&Psi, &ptdat, &bg, q, k, l_max, tau, final_tau, M, 1000000, dlnf0_dlnq);
        printf("%f %e %e %e %e %e\n", tau, Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
    }


    // final_tau = 100;
    //
    // evolve(&Psi, &ptdat, &bg, q, k, l_max, tau, final_tau, M);
    //
    // printf("We found %e %e %e %e %e\n", Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
    bg_interp_free(&bg);
    rend_interp_free(&ptdat);


    free(Psi);

    /* Clean up */
    cleanPerturb(&ptdat);
    cleanBackground(&bg);
    cleanParams(&pars);

}
