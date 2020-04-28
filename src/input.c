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
#include <assert.h>
#include "../include/input.h"

int readParams(struct params *pars, const char *fname) {
     // pars->GridSize = ini_getl("Box", "GridSize", 64, fname);
     // pars->BoxLen = ini_getd("Box", "BoxLen", 1.0, fname);
     // pars->Homogeneous = ini_getbool("Simulation", "Homogeneous", 0, fname);

     /* Read strings */
     int len = DEFAULT_STRING_LENGTH;
     pars->OutputDirectory = malloc(len);
     pars->OutputFilename = malloc(len);
     pars->Name = malloc(len);
     pars->BackgroundFile = malloc(len);
     pars->BackgroundFormat = malloc(len);
     pars->PerturbFile = malloc(len);
     ini_gets("Output", "Directory", "./output", pars->OutputDirectory, len, fname);
     ini_gets("Simulation", "Name", "No Name", pars->Name, len, fname);
     ini_gets("Background", "File", "", pars->BackgroundFile, len, fname);
     ini_gets("Background", "Format", "Plain", pars->BackgroundFormat, len, fname);
     ini_gets("PerturbData", "File", "", pars->PerturbFile, len, fname);
     ini_gets("Output", "Filename", "perturb_out.hdf5", pars->OutputFilename, len, fname);

     pars->MaxMultipole = ini_getl("Simulation", "MaxMultipole", 50, fname);
     pars->MaxMomentum = ini_getd("Simulation", "MaxMomentum", 15, fname);
     pars->NumberMomentumBins = ini_getl("Simulation", "NumberMomentumBins", 28, fname);
     pars->InitialTime = ini_getd("Simulation", "InitialTime", 0.05, fname);
     pars->Tolerance = ini_getd("Simulation", "Tolerance", 1e-10, fname);
     pars->Verbose = ini_getl("Simulation", "Verbose", 0, fname);

     /* Single run parameters */
     pars->kSingle = ini_getd("Single", "k", 0.01, fname);
     pars->tauFinalSingle = ini_getd("Single", "tauFinal", 1000, fname);

     return 0;
}

int readUnits(struct units *us, const char *fname) {
    /* Internal units */
    us->UnitLengthMetres = ini_getd("Units", "UnitLengthMetres", 1.0, fname);
    us->UnitTimeSeconds = ini_getd("Units", "UnitTimeSeconds", 1.0, fname);
    us->UnitMassKilogram = ini_getd("Units", "UnitMassKilogram", 1.0, fname);
    us->UnitTemperatureKelvin = ini_getd("Units", "UnitTemperatureKelvin", 1.0, fname);

    /* Physical constants */
    us->SpeedOfLight = SPEED_OF_LIGHT_METRES_SECONDS * us->UnitTimeSeconds
                        / us->UnitLengthMetres;
    us->GravityG = GRAVITY_G_SI_UNITS * us->UnitTimeSeconds * us->UnitTimeSeconds
                    / us->UnitLengthMetres / us->UnitLengthMetres / us->UnitLengthMetres
                    * us->UnitMassKilogram; // m^3 / kg / s^2 to internal
    us->hPlanck = PLANCK_CONST_SI_UNITS / us->UnitMassKilogram / us->UnitLengthMetres
                    / us->UnitLengthMetres * us->UnitTimeSeconds; //J*s = kg*m^2/s
    us->kBoltzmann = BOLTZMANN_CONST_SI_UNITS / us->UnitMassKilogram / us->UnitLengthMetres
                    / us->UnitLengthMetres * us->UnitTimeSeconds * us->UnitTimeSeconds
                    * us->UnitTemperatureKelvin; //J/K = kg*m^2/s^2/K
    us->ElectronVolt = ELECTRONVOLT_SI_UNITS / us->UnitMassKilogram / us->UnitLengthMetres
                    / us->UnitLengthMetres * us->UnitTimeSeconds
                    * us->UnitTimeSeconds; // J = kg*m^2/s^2


    us->BackgroundUnitLengthMetres = ini_getd("TransferFunctions", "UnitLengthMetres", MPC_METRES, fname);

    return 0;
}

int readCosmology(struct cosmology_params *cosmo, const char *fname) {
     cosmo->h = ini_getd("Cosmology", "h", 0.70, fname);
     cosmo->T_nu0 = ini_getd("Cosmology", "T_nu0", 1.951757805, fname);
     cosmo->M_nu = ini_getd("Cosmology", "M_nu", 0.02, fname);
     cosmo->Degeneracy = ini_getd("Cosmology", "Degeneracy", 1, fname);

     return 0;
}

int cleanParams(struct params *pars) {
    free(pars->OutputDirectory);
    free(pars->OutputFilename);
    free(pars->Name);
    free(pars->BackgroundFile);
    free(pars->BackgroundFormat);
    free(pars->PerturbFile);

    return 0;
}


int readGRF_H5(double **box, int *N, double *box_len, const char *fname) {
    /* Create the hdf5 file */
    hid_t h_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

    /* Create the Header group */
    hid_t h_grp = H5Gopen(h_file, "Header", H5P_DEFAULT);

    /* Read the size of the field */
    hid_t h_attr, h_err;
    double boxsize[3];

    /* Open and read out the attribute */
    h_attr = H5Aopen(h_grp, "BoxSize", H5P_DEFAULT);
    h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &boxsize);
    assert(h_err >= 0);

    /* It should be a cube */
    assert(boxsize[0] == boxsize[1]);
    assert(boxsize[1] == boxsize[2]);
    *box_len = boxsize[0];

    /* Close the attribute, and the Header group */
    H5Aclose(h_attr);
    H5Gclose(h_grp);

    /* Open the Field group */
    h_grp = H5Gopen(h_file, "Field", H5P_DEFAULT);

    /* Open the Field dataset */
    hid_t h_data = H5Dopen2(h_grp, "Field", H5P_DEFAULT);

    /* Open the dataspace and fetch the grid dimensions */
    hid_t h_space = H5Dget_space(h_data);
    int ndims = H5Sget_simple_extent_ndims(h_space);
    hsize_t *dims = malloc(ndims * sizeof(hsize_t));
    H5Sget_simple_extent_dims(h_space, dims, NULL);

    /* We should be in 3D */
    assert(ndims == 3);
    /* It should be a cube */
    assert(dims[0] == dims[1]);
    assert(dims[1] == dims[2]);
    *N = dims[0];

    /* Allocate the array */
    *box = malloc(dims[0] * dims[1] * dims[2] * sizeof(double));

    /* Read out the data */
    h_err = H5Dread(h_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *box);

    /* Close the dataspace and dataset */
    H5Sclose(h_space);
    H5Dclose(h_data);
    free(dims);

    /* Close the Field group */
    H5Gclose(h_grp);

    /* Close the file */
    H5Fclose(h_file);

    return 0;
}
