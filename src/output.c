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

#include <string.h>
#include "../include/output.h"

int write_perturb(struct perturb_data *data, struct params *pars,
                  struct units *us, char *fname) {
    /* The memory for the transfer functions is located here */
    hid_t h_file, h_grp, h_data, h_err, h_space, h_attr;

    /* Open file */
    h_file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h_file < 0) printf("Error while opening file '%s'.\n", fname);

    printf("Writing the perturbation to '%s'.\n", fname);

    /* Open header to write simulation properties */
    h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp < 0) printf("Error while creating file header\n");

    /* We will write attributes, each consisting of a single scalar number */
    int ndims = 1;
    hsize_t dim[1] = {1};

    /* Create dataspace */
    h_space = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(h_space, ndims, dim, NULL);

    /* Write the size of the perturbation (along the k-dimension) */
    h_attr = H5Acreate1(h_grp, "k_size", H5T_NATIVE_INT, h_space, H5P_DEFAULT);
    h_err = H5Awrite(h_attr, H5T_NATIVE_INT, &data->k_size);
    H5Aclose(h_attr);

    /* Write the size of the perturbation (along the tau-direction) */
    h_attr = H5Acreate1(h_grp, "tau_size", H5T_NATIVE_INT, h_space, H5P_DEFAULT);
    h_err = H5Awrite(h_attr, H5T_NATIVE_INT, &data->tau_size);
    H5Aclose(h_attr);

    /* Write the number of transfer functions */
    h_attr = H5Acreate1(h_grp, "n_functions", H5T_NATIVE_INT, h_space, H5P_DEFAULT);
    h_err = H5Awrite(h_attr, H5T_NATIVE_INT, &data->n_functions);
    H5Aclose(h_attr);

    /* Determine the units used */
    double unit_mass_cgs = us->UnitMassKilogram * 1000;
    double unit_length_cgs = us->UnitLengthMetres * 100;
    double unit_time_cgs = us->UnitTimeSeconds;

    /* Write the internal unit system */
    h_attr = H5Acreate1(h_grp, "Unit mass in cgs (U_M)", H5T_NATIVE_DOUBLE, h_space, H5P_DEFAULT);
    h_err = H5Awrite(h_attr, H5T_NATIVE_DOUBLE, &unit_mass_cgs);
    H5Aclose(h_attr);

    h_attr = H5Acreate1(h_grp, "Unit length in cgs (U_L)", H5T_NATIVE_DOUBLE, h_space, H5P_DEFAULT);
    h_err = H5Awrite(h_attr, H5T_NATIVE_DOUBLE, &unit_length_cgs);
    H5Aclose(h_attr);

    h_attr = H5Acreate1(h_grp, "Unit time in cgs (U_t)", H5T_NATIVE_DOUBLE, h_space, H5P_DEFAULT);
    h_err = H5Awrite(h_attr, H5T_NATIVE_DOUBLE, &unit_time_cgs);
    H5Aclose(h_attr);

    /* For strings, we need to prepare a datatype */
    const hid_t h_type = H5Tcopy(H5T_C_S1);
    h_err = H5Tset_size(h_type, strlen(pars->Name)); //length of simname

    /* Create another attribute, which is the name of the simulation */
    h_attr = H5Acreate1(h_grp, "Name", h_type, h_space, H5P_DEFAULT);

    /* Write the name attribute */
    h_err = H5Awrite(h_attr, h_type, pars->Name);

    /* Done with the single entry dataspace */
    H5Sclose(h_space);


    /* Write array of column titles, corresponding to the exported functions */
    hsize_t dim_columns[1] = {data->n_functions};

    /* Create dataspace */
    h_space = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(h_space, ndims, dim_columns, NULL);

    /* Set string length to be variable */
    h_err = H5Tset_size(h_type, H5T_VARIABLE); //length of strings

    /* Create another attribute, which is the name of the simulation */
    h_attr = H5Acreate1(h_grp, "FunctionTitles", h_type, h_space, H5P_DEFAULT);

    /* Write the name attribute */
    h_err = H5Awrite(h_attr, h_type, data->titles);


    /* Done with the dataspace */
    H5Sclose(h_space);

    /* Close header */
    H5Gclose(h_grp);

    /* Open group to write the perturbation arrays */
    h_grp = H5Gcreate(h_file, "/Perturb", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp < 0) printf("Error while creating perturbation group\n");

    /* Create data space */
    h_space = H5Screate(H5S_SIMPLE);
    if (h_space < 0) printf("Error while creating data space.\n");

    /* Set the extent of the data */
    int rank = 1;
    hsize_t shape[1] = {data->k_size};
    h_err = H5Sset_extent_simple(h_space, rank, shape, shape);
    if (h_err < 0) printf("Error while changing data space shape.\n");

    /* Dataset properties */
    const hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);

    /* Create dataset */
    h_data = H5Dcreate(h_grp, "Wavenumbers", H5T_NATIVE_DOUBLE, h_space,
                       H5P_DEFAULT, h_prop, H5P_DEFAULT);
    if (h_data < 0) printf("Error while creating dataspace '%s'.\n", "Wavenumbers");

    /* Write temporary buffer to HDF5 dataspace */
    h_err = H5Dwrite(h_data, H5T_NATIVE_DOUBLE, h_space, H5S_ALL, H5P_DEFAULT, data->k);
    if (h_err < 0) printf("Error while writing data array '%s'.\n", "data->k");

    /* Close the dataset */
    H5Dclose(h_data);

    /* Set the extent of the tau data */
    rank = 1;
    hsize_t shape_tau[1] = {data->tau_size};
    h_err = H5Sset_extent_simple(h_space, rank, shape_tau, shape_tau);
    if (h_err < 0) printf("Error while changing data space shape.\n");

    /* Create dataset */
    h_data = H5Dcreate(h_grp, "Log conformal times", H5T_NATIVE_DOUBLE,
                       h_space, H5P_DEFAULT, h_prop, H5P_DEFAULT);
    if (h_data < 0)
    printf("Error while creating dataspace '%s'.\n", "Log conformal times");

    /* Write temporary buffer to HDF5 dataspace */
    h_err = H5Dwrite(h_data, H5T_NATIVE_DOUBLE, h_space, H5S_ALL, H5P_DEFAULT, data->log_tau);
    if (h_err < 0) printf("Error while writing data array '%s'.\n", "data->log_tau");

    /* Close the dataset */
    H5Dclose(h_data);

    /* Set the extent of the transfer function data */
    rank = 3;
    hsize_t shape_delta[3] = {data->n_functions, data->k_size, data->tau_size};
    h_err = H5Sset_extent_simple(h_space, rank, shape_delta, shape_delta);
    if (h_err < 0) printf("Error while changing data space shape.");

    /* Create dataset */
    h_data = H5Dcreate(h_grp, "Transfer functions", H5T_NATIVE_DOUBLE, h_space,
                       H5P_DEFAULT, h_prop, H5P_DEFAULT);
    if (h_data < 0)
    printf("Error while creating dataspace '%s'.", "Transfer functions");

    /* Write temporary buffer to HDF5 dataspace */
    h_err = H5Dwrite(h_data, H5T_NATIVE_DOUBLE, h_space, H5S_ALL, H5P_DEFAULT, data->delta);
    if (h_err < 0) printf("Error while writing data array '%s'.", "data->delta");

    /* Close the dataset */
    H5Dclose(h_data);

    /* Close the properties */
    H5Pclose(h_prop);

    /* Close the group */
    H5Gclose(h_grp);

    /* Close file */
    H5Fclose(h_file);

    return 0;
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
