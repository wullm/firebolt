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

#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "../include/grids.h"
#include "../include/fft.h"
#include "../include/multipole_interp.h"

static inline void kernel_transfer_function(struct kernel *the_kernel) {
    double k = the_kernel->k;

    if (k == 0) {
        the_kernel->kern = 0.f;
    } else {
        the_kernel->kern = multipoleInterp(k);
    }
}

/* Generate a density grid for each particle type by applying the power
 * spectrum to the random phases. The necessary transfer functions are in trs.
 */
int generateGrids(const struct multipoles *m,
                  const fftw_complex *grf,
                  struct grids *grs) {

    /* Grid dimensions */
    const int N = grs->N;
    const double boxlen = grs->boxlen;

    /* Create complex and real 3D arrays */
    fftw_complex *fbox = (fftw_complex*) fftw_malloc(N*N*(N/2+1)*sizeof(fftw_complex));
    double *box = (double*) fftw_malloc(N*N*N*sizeof(double));

    /* Create FFT plan */
    fftw_plan c2r = fftw_plan_dft_c2r_3d(N, N, N, fbox, box, FFTW_ESTIMATE);

    /* For each multipole/momentum bin pair, create the corresponding grid */
    for (int index_q=0; index_q<m->q_size; index_q++) {
        for (int index_l=0; index_l<m->l_size; index_l++) {

            /* Switch the interpolation spline to this transfer function */
            switchMultipoleInterp(m, index_l, index_q);

            /* Apply the transfer function to grf and store in fbox */
            fft_apply_kernel(fbox, grf, N, boxlen, kernel_transfer_function);

            // /* Compute the pre-factor = (2l+1) * (-1)^l */
            // double factor = (2*index_l + 1) * pow(-1, index_l);
            // /* The factor i^l is accounted for by taking derivatives later */
            //
            // /* Apply the pre-factor */
            // for (int i=0; i<N*N*(N/2+1); i++) {
            //     fbox[i] *= factor;
            // }

            /* FFT and normalize */
            fft_execute(c2r);
            fft_normalize_c2r(box,N,boxlen);

            /* Store the grid */
            double *destination = grs->grids + index_l * (N*N*N) * m->q_size + index_q * (N*N*N);
            double *source_address = box;
            memcpy(destination, source_address, (N*N*N) * sizeof(double));

            // /* Export the real box */
            // char dbox_fname[DEFAULT_STRING_LENGTH];
            // sprintf(dbox_fname, "grid_l%d_q%d.hdf5", index_l, index_q);
            // // fft_c2r_export(fbox, N, boxlen, dbox_fname);
            // writeGRF_H5(box, N, boxlen, dbox_fname);
            // printf("Density field exported to '%s'.\n", dbox_fname);
        }
    }

    /* Free all the FFT objects */
    fftw_free(fbox);
    fftw_free(box);
    fftw_destroy_plan(c2r);

    return 0;
}

int initGrids(int N, double boxlen, const struct multipoles *m, struct grids *grs) {

    /* Grid dimensions */
    grs->N = N;
    grs->boxlen = boxlen;
    int totalGridSize = (grs->N) * (grs->N) * (grs->N);

    /* Number of grids */
    grs->q_size = m->q_size;
    grs->l_size = m->l_size;
    int numberOfGrids = (grs->q_size) * (grs->l_size);

    /* Allocate memory */
    grs->grids = malloc(numberOfGrids * totalGridSize * sizeof(double));

    if (grs->grids == NULL) {
        printf("Error allocating memory for grids.");
        return 1;
    }

    return 0;
}

int cleanGrids(struct grids *grs) {
    free(grs->grids);

    return 0;
}
