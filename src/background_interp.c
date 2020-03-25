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
#include <gsl/gsl_spline.h>

#include "../include/background_interp.h"

/* GSL interpolation objects */
const gsl_interp_type *bg_interp_type;
gsl_interp_accel *bg_tau_acc;
gsl_spline *bg_spline;

int bg_interp_init(const struct background *bg) {
    /* We will use linear interpolation in log-tau space */
    bg_interp_type = gsl_interp_linear;

    /* Allocate memory for the spline */
    bg_spline = gsl_spline_alloc(bg_interp_type, bg->nrow);
    /* Note: this only copies the first transfer function from bg->delta */
    gsl_spline_init(bg_spline, bg->log_tau, bg->z, bg->nrow);


    /* Allocate memory for the accelerator objects */
    bg_tau_acc = gsl_interp_accel_alloc();

  return 0;
}

int bg_interp_free(const struct background *bg) {
    /* Done with the GSL interpolation */
    gsl_spline_free(bg_spline);
    gsl_interp_accel_free(bg_tau_acc);

    return 0;
}

double bg_z_at_tau(double log_tau) {
    return gsl_spline_eval(bg_spline, log_tau, bg_tau_acc);
}
