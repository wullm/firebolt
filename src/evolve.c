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
#include <assert.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "../include/ic.h"
#include "../include/background_interp.h"
#include "../include/perturb_interp.h"
#include "../include/evolve.h"


int func (double tau, const double Psi[], double dPsi[], void *params) {
    double *pars = (double*)params;
    double k = pars[0];
    double q = pars[1];
    double M = pars[2];
    int l_max = (int) pars[3];
    double dlnf0_dlnq = pars[4];

    double logt = log(tau);
    double z = bg_z_at_log_tau(logt);
    double a = 1./(1+z);
    double eps = hypot(q, a*M);

    double qke = q*k/eps;

    double h_prime = rend_interp(k, logt, 0);
    double eta_prime = rend_interp(k, logt, 1);

    /* Set the final Psi, using the truncation prescription of Ma & Bertschinger */
    double Psi_lmax_p1 = (2*l_max+1)/qke/tau * Psi[l_max] - Psi[l_max-1];

    /* Compute the derivatives */
    dPsi[0] = - qke * Psi[1] + h_prime/6 * dlnf0_dlnq;
    dPsi[1] =   qke/3 * (Psi[0] - 2*Psi[2]);
    dPsi[2] =   qke/5 * (2 * Psi[1] - 3 * Psi[3]) - (h_prime/15 + 2*eta_prime/5) * dlnf0_dlnq;
    for (int l=3; l<l_max; l++) {
        dPsi[l] = qke/(2*l+1.) * (l*Psi[l-1] - (l+1)*Psi[l+1]);
    }    
    dPsi[l_max] = qke/(2*l_max+1.) * (l_max*Psi[l_max-1] - (l_max+1)*Psi_lmax_p1);

    return GSL_SUCCESS;
}

int evolve_gsl(double **Psi, const struct perturb_data *ptdat, const struct background *bg,
               double q, double k, int l_max, double tau_ini, double tau_final, double M,
               double dlnf0_dlnq) {

    /* Coldsonte */
    double params[5];
    params[0] = k;
    params[1] = q;
    params[2] = M;
    params[3] = l_max;
    params[4] = dlnf0_dlnq;

    gsl_odeiv2_system sys = {func, NULL, l_max+1, params};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  1e-14, 1e-14, 1e-14);
    double t = tau_ini, t1 = tau_final;

    int status = gsl_odeiv2_driver_apply(d, &t, t1, *Psi);

    gsl_odeiv2_driver_free(d);

    return status;
}

int evolve(double **Psi, const struct perturb_data *ptdat, const struct background *bg,
           double q, double k, int l_max, double tau_ini, double tau_final, double M,
           int steps, double dlnf0_dlnq) {

    /* Create tables of h_prime & eta_prime at the needed time steps */
    double tau_factor = exp(log(tau_final/tau_ini)/(steps - 1));

    double *tau = malloc(steps * sizeof(double));
    double *h_prime = malloc(steps * sizeof(double));
    double *eta_prime = malloc(steps * sizeof(double));

    /* Create the time table */
    tau[0] = tau_ini;
    for (int i=1; i<steps; i++) {
        tau[i] = tau[i-1] * tau_factor;
    }

    /* Switch to the h_prime transfer function */
    rend_interp_switch_source(ptdat, 0, 0);
    /* Create the h_prime table */
    for (int i=0; i<steps; i++) {
        h_prime[i] = rend_interp(k, log(tau[i]), 0);
    }

    /* Switch to the eta_prime transfer function */
    rend_interp_switch_source(ptdat, 1, 1);
    /* Create the h_prime table */
    for (int i=0; i<steps; i++) {
        eta_prime[i] = rend_interp(k, log(tau[i]), 1);
    }

    /* Allocate memory for the conformal time derivatives of Psi */
    double *dPsi = malloc((l_max + 1) * sizeof(double));

    /* Integrate the equations */
    for (int i=1; i<steps; i++) {
        /* The redshift */
        double z = bg_z_at_log_tau(log(tau[i]));
        double a = 1./(1+z);
        double eps = hypot(q, a*M);

        dPsi[0] = - q*k/eps * (*Psi)[1] + h_prime[i]/6 * dlnf0_dlnq;
        dPsi[1] =   q*k/eps/3 * ((*Psi)[0] - 2*(*Psi)[2]);
        dPsi[2] =   q*k/eps/5 * (2 * (*Psi)[1] - 3 * (*Psi)[3]) - (h_prime[i]/15 + 2*eta_prime[i]/5) * dlnf0_dlnq;
        for (int l=3; l<l_max; l++) {
            dPsi[l] =   q*k/eps/(2*l+1.) * (l*(*Psi)[l-1] - (l+1)*(*Psi)[l+1]);
        }
        dPsi[l_max] = 0;

        /* Update the multipoles */
        for (int l=0; l<l_max+1; l++) {
            (*Psi)[l] += dPsi[l] * (tau[i] - tau[i-1]);
        }
    }

    free(dPsi);


    return 0;
}
