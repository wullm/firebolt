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
#include <gsl/gsl_matrix.h>
#include "../include/ic.h"
#include "../include/background_interp.h"
#include "../include/perturb_interp.h"
#include "../include/evolve.h"


int func(double tau, const double Psi[], double dPsi[], void *params) {
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

/* No longer needed... */
// int Jac(double tau, const double Psi[], double *J, double dfdt[], void *params) {
//     double *pars = (double*)params;
//     double k = pars[0];
//     double q = pars[1];
//     double M = pars[2];
//     int l_max = (int) pars[3];
//     double dlnf0_dlnq = pars[4];
//
//     double logt = log(tau);
//     double z = bg_z_at_log_tau(logt);
//     double a = 1./(1+z);
//     double eps = hypot(q, a*M);
//
//     double qke = q*k/eps;
//
//     double h_prime = rend_interp(k, logt, 0);
//     double eta_prime = rend_interp(k, logt, 1);
//
//     double h_prime_prime = rend_dydt_interp(k, logt, 0);
//     double eta_prime_prime = rend_dydt_interp(k, logt, 1);
//
//     J[0 * (l_max+1) + 1] = -qke;
//     J[1 * (l_max+1) + 0] = qke/3;
//     J[1 * (l_max+1) + 2] = -2*qke/3;
//     J[2 * (l_max+1) + 1] = 2*qke/5;
//     J[2 * (l_max+1) + 3] = -3*qke/5;
//     for (int l=3; l<l_max; l++) {
//         J[l * (l_max+1) + l-1] = l*qke/(2*l+1.);
//         J[l * (l_max+1) + l+1] = -(l+1)*qke/(2*l+1.);
//     }
//     J[l_max * (l_max+1) + l_max-1] = l_max*qke/(2*l_max+1.) + (l_max+1)*qke/(2*l_max+1.);
//     J[l_max * (l_max+1) + l_max] = -(l_max+1)/tau;
//
//     dfdt[0] = h_prime_prime/6 * dlnf0_dlnq;
//     dfdt[1] = 0.;
//     dfdt[2] = -(h_prime_prime/15 + 2*eta_prime_prime/5) * dlnf0_dlnq;
//     for (int l=3; l<l_max+1; l++) {
//         dfdt[l] = 0.;
//     }
//
//     return GSL_SUCCESS;
// }

int evolve_gsl(double **Psi, const struct perturb_data *ptdat, const struct background *bg,
               double q, double k, int l_max, double tau_ini, double tau_final, double M,
               double dlnf0_dlnq, double tolerance) {

    /* Coldsonte */
    double params[5];
    params[0] = k;
    params[1] = q;
    params[2] = M;
    params[3] = l_max;
    params[4] = dlnf0_dlnq;
    gsl_odeiv2_system sys = {func, NULL, l_max+1, params};
    double t = tau_ini, t1 = tau_final;

    double h = 1e-10;
    double err_abs = tolerance;
    double err_rel = tolerance;

    const gsl_odeiv2_step_type *type = gsl_odeiv2_step_rkck;
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(type, l_max+1);
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(err_abs, err_rel);
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(l_max+1);
    gsl_odeiv2_driver *d  = gsl_odeiv2_driver_alloc_y_new (&sys, type, h, err_abs, err_rel);
    gsl_odeiv2_step_set_driver(s,d);

    int status = 0;
    while(t<t1) {
        status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, *Psi);
        if (status != GSL_SUCCESS){
            break;
        }
    }

    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    gsl_odeiv2_driver_free(d);

    return status;
}
