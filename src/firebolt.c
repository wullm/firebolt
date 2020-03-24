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

    readParams(&pars, fname);
    readUnits(&us, fname);
    readCosmology(&cosmo, fname);
    readBackground(&pars, &us, &cosmo, &bg);

    double q = 1.0;
    double k = 1.0;
    double tau = exp(-0.29);
    printf("[q, k, tau] = [%f, %f, %f]\n", q, k, tau);

    int l_max = 10;
    double *Psi;

    generate_ics(&bg, q, k, tau, &Psi, l_max);

    printf("We found %f %f %f %f %f\n", Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

    free(Psi);

    /* Clean up */
    cleanBackground(&bg);
    cleanParams(&pars);

}
