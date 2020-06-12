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

#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <math.h>
#include "../include/perturb_data.h"

int readPerturb(struct params *pars, struct units *us, struct perturb_data *pt) {
    printf("Reading the perturbation from '%s'.\n", pars->PerturbFile);

    /* Open the hdf5 file (file exists error handled by HDF5) */
    hid_t h_file = H5Fopen(pars->PerturbFile, H5F_ACC_RDONLY, H5P_DEFAULT);

    /* Open the Header group */
    hid_t h_grp = H5Gopen(h_file, "Header", H5P_DEFAULT);

    /* Read the size of the perturbation */
    hid_t h_attr, h_err, h_tp;

    /* Number of wavenumbers */
    h_attr = H5Aopen(h_grp, "k_size", H5P_DEFAULT);
    h_err = H5Aread(h_attr, H5T_NATIVE_INT, &pt->k_size);
    H5Aclose(h_attr);

    /* Number of time steps */
    h_attr = H5Aopen(h_grp, "tau_size", H5P_DEFAULT);
    h_err = H5Aread(h_attr, H5T_NATIVE_INT, &pt->tau_size);
    H5Aclose(h_attr);

    /* Number of transfer functions T(k,tau) */
    h_attr = H5Aopen(h_grp, "n_functions", H5P_DEFAULT);
    h_err = H5Aread(h_attr, H5T_NATIVE_INT, &pt->n_functions);
    H5Aclose(h_attr);

    /* Check that it makes sense */
    if (pt->k_size <= 0 || pt->tau_size <= 0 || pt->n_functions <= 0) {
        printf("ERROR: reading an empty perturbation file.\n");
        return 1;
    }

    /* Allocate memory for the transfer function titles */
    pt->titles = malloc(pt->n_functions * sizeof(char*));

    /* Read the titles of the transfer functions */
    h_attr = H5Aopen(h_grp, "FunctionTitles", H5P_DEFAULT);
    h_tp = H5Aget_type(h_attr);
    h_err = H5Aread(h_attr, h_tp, pt->titles);
    H5Aclose(h_attr);
    H5Tclose(h_tp);

    /* Read the units of the perturbation data file */
    double UnitLengthCGS, UnitTimeCGS, UnitMassCGS;
    double UnitLengthMetres, UnitTimeSeconds, UnitMassKilogram;

    /* CGS Length units */
    h_attr = H5Aopen(h_grp, "Unit length in cgs (U_L)", H5P_DEFAULT);
    h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &UnitLengthCGS);
    H5Aclose(h_attr);

    /* CGS Time units */
    h_attr = H5Aopen(h_grp, "Unit time in cgs (U_t)", H5P_DEFAULT);
    h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &UnitTimeCGS);
    H5Aclose(h_attr);

    /* CGS Mass units */
    h_attr = H5Aopen(h_grp, "Unit mass in cgs (U_M)", H5P_DEFAULT);
    h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &UnitMassCGS);
    H5Aclose(h_attr);

    /* Convert to SI units */
    UnitLengthMetres = UnitLengthCGS/100;
    UnitTimeSeconds = UnitTimeCGS;
    UnitMassKilogram = UnitMassCGS/1000;

    /* Check that it makes sense */
    if (UnitLengthMetres <= 0 || UnitTimeSeconds <= 0 || UnitMassKilogram <= 0) {
        printf("ERROR: unknown units of perturbation file.\n");
        return 1;
    }

    /* Close the Header group */
    H5Gclose(h_grp);

    /* Open the data group */
    h_grp = H5Gopen(h_file, "Perturb", H5P_DEFAULT);

    /* Allocate memory */
    pt->k = calloc(pt->k_size, sizeof(double));
    pt->log_tau = calloc(pt->tau_size, sizeof(double));
    pt->redshift = calloc(pt->tau_size, sizeof(double));
    pt->delta = malloc(pt->n_functions * pt->k_size * pt->tau_size * sizeof(double));
    // pt->dydt = malloc(pt->n_functions * pt->k_size * pt->tau_size * sizeof(double));

    /* Dataspace */
    hid_t h_data;

    /* Allocation successful? */
    if (pt->k == NULL || pt->log_tau == NULL || pt->delta == NULL || pt->redshift == NULL) {
        printf("ERROR: unable to allocate memory for perturbation data.");
    }

    /* Read the wavenumbers */
    h_data = H5Dopen2(h_grp, "Wavenumbers", H5P_DEFAULT);
    h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pt->k);
    H5Dclose(h_data);

    /* Read the logs of conformal times */
    h_data = H5Dopen2(h_grp, "Log conformal times", H5P_DEFAULT);
    h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pt->log_tau);
    H5Dclose(h_data);

    /* Read the redshifts */
    h_data = H5Dopen2(h_grp, "Redshifts", H5P_DEFAULT);
    h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pt->redshift);
    H5Dclose(h_data);

    /* Read the transfer functions */
    h_data = H5Dopen2(h_grp, "Transfer functions", H5P_DEFAULT);
    h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pt->delta);
    H5Dclose(h_data);

    if (h_err < 0 || pt->k[0] == 0 || pt->log_tau[0] == 0) {
        printf("ERROR: problem with reading the perturbation data.\n");
        return 1;
    }

    /* Close the data group */
    H5Gclose(h_grp);

    /* Close the file */
    H5Fclose(h_file);

    /* Perform unit conversions for the wavenumbers */
    for (int i=0; i<pt->k_size; i++) {
        pt->k[i] *= us->UnitLengthMetres / UnitLengthMetres;
    }

    /* Perform unit conversions for the conformal times */
    for (int i=0; i<pt->tau_size; i++) {
        pt->log_tau[i] += log(UnitTimeSeconds/us->UnitTimeSeconds);
    }

    const double unit_length_factor = UnitLengthMetres / us->UnitLengthMetres;
    const double unit_time_factor = UnitTimeSeconds / us->UnitTimeSeconds;

    printf("Velocity factor = %e\n", 1./unit_time_factor);

    /* Perform unit conversions for the transfer functions */
    for (int i=0; i<pt->n_functions; i++) {
        /* Determine the unit conversion factor */
        char *title = pt->titles[i];
        double unit_factor = unitConversionFactor(title, unit_length_factor, unit_time_factor);

        printf("Unit conversion factor for '%s' is %f\n", title, unit_factor);

        /* Convert from input units to internal units */
        for (int index_k=0; index_k<pt->k_size; index_k++) {
            for (int index_tau=0; index_tau<pt->tau_size; index_tau++) {
                int index = pt->tau_size * pt->k_size * i + pt->k_size * index_tau + index_k;
                pt->delta[index] *= unit_factor;
            }
        }
    }

    return 0;
}

/* Unit conversion factor for transfer functions, depending on the title. */
double unitConversionFactor(char *title, double unit_length_factor,
                            double unit_time_factor) {

    /* Note the difference between strcmp and strncmp! */
    char *title_end = &title[strlen(title)];

    /* Most transfer functions are dimensionless, (e.g. overdensities) */
    double factor = 1.0;

    /* Energy flux transfer functions (theta = nabla.v) have dimension
    * inverse time. These have titles starting with "t_". */
    if (strncmp(title, "t_", 2) == 0) {
        factor /= unit_time_factor;
    }

    /* Functions that are time derivatives have dimension inverse time */
    if (strncmp(title_end-12, "_prime_prime", 12) == 0) {
        factor /= pow(unit_time_factor, 2);
    } else if (strncmp(title_end-6, "_prime", 6) == 0) {
        factor /= unit_time_factor;
    }

    /* Potential transfer functions have dimensions of energy per mass or
    * (Length/Time)^2. This applies to "phi", "psi", "h", "eta" and their
    * derivatives (and possibly higher order derivatives).
    */
    if (strcmp(title, "h") == 0   || strncmp(title, "h_", 2) == 0 ||
        strcmp(title, "phi") == 0 || strncmp(title, "phi_", 4) == 0 ||
        strcmp(title, "eta") == 0 || strncmp(title, "eta_", 4) == 0 ||
        strcmp(title, "psi") == 0 || strncmp(title, "psi_", 4) == 0) {
            factor *= pow(unit_length_factor / unit_time_factor, 2);
    }

    return factor;
}

int cleanPerturb(struct perturb_data *pt) {
    free(pt->k);
    free(pt->log_tau);
    free(pt->redshift);
    free(pt->delta);
    // free(pt->dydt);
    for (int i=0; i<pt->n_functions; i++) {
        free(pt->titles[i]);
    }
    free(pt->titles);

    return 0;
}
