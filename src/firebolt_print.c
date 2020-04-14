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
    double tau_ini = exp(ptdat.log_tau[0]);
    double tau_fin = pars.tauFinalSingle;

    printf("\n");
    printf("[tau_ini, z_ini] = [%f, %.2e]\n", tau_ini, bg_z_at_log_tau(log(tau_ini)));
    printf("[tau_fin, z_fin] = [%f, %.2e]\n", tau_fin, bg_z_at_log_tau(log(tau_fin)));
    printf("\n");

    double log_tau_fin = log(tau_fin);

    for (int i=0; i<ptdat.k_size; i++) {
        double k = ptdat.k[i];

        /* Retrieve delta and theta results */
        rend_interp_switch_source(&ptdat, 4, 0);
        rend_interp_switch_source(&ptdat, 5, 1);

        double delta_nu = rend_interp(k, log_tau_fin, 0);
        double theta_nu = rend_interp(k, log_tau_fin, 1);

        /* Also determine the shear and l3*/
        rend_interp_switch_source(&ptdat, 6, 0);
        double shear_nu = rend_interp(k, log_tau_fin, 0);
        rend_interp_switch_source(&ptdat, 8, 1);
        double l3_nu = rend_interp(k, log_tau_fin, 1);

        /* Also determine the sound speed cs2 */
        rend_interp_switch_source(&ptdat, 7, 1);
        double cs2_nu = rend_interp(k, log_tau_fin, 1);

        printf("%f %e %e %e %e %e\n", k, delta_nu, theta_nu, shear_nu, l3_nu, cs2_nu);

    }


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
