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
#include "../include/evolve.h"

struct ode_pars {
    /* Array of parameters of the diferential equation */
    double *params;
    /* Reference to function that gives redshift z(log_tau) */
    double (*redshift_func)(double log_tau);
    /* Reference to function that gives h'(k, log_tau) */
    double (*h_prime_func)(double k, double log_tau);
    /* Reference to function that gives eta'(k, log_tau) */
    double (*eta_prime_func)(double k, double log_tau);
};

struct ode_pars op;

/* Right-hand side of the differential equation */
int func(double tau, const double Psi[], double dPsi[], void *ode_pars) {
    /* Unpack the parameters of the differential equation */
    struct ode_pars *ops = (struct ode_pars*) ode_pars;
    double *params = ops->params;
    double k = params[0]; //wavenumber
    double q = params[1]; //dimensionless momentum
    double M = params[2]; //mass
    int l_max = (int) params[3]; //maximum multipole
    double dlnf0_dlnq = params[4]; //logarithmic derivative of distr. func.
    double c_vel = params[5]; //speed of light

    /* Compute the necessary quantities */
    double logt = log(tau);
    // double z = perturb_zAtLogTau(logt); //redshift
    double z = ops->redshift_func(logt); //redshift
    double a = 1./(1+z); //scale factor
    double eps = hypot(q, a*M); //dimensionless energy
    double qke = (q*k/eps) * c_vel;
    double cc = c_vel * c_vel;

    /* Interpolate the source functions at (k, log tau) */
    double h_prime = ops->h_prime_func(k, logt);
    double eta_prime = ops->eta_prime_func(k, logt);

    /* Set the final Psi, using the truncation prescription of Ma & Bertschinger */
    double Psi_lmax_p1 = (2*l_max+1)/qke/tau * Psi[l_max] - Psi[l_max-1];

    /* Compute the derivatives */
    dPsi[0] = - qke * Psi[1] + h_prime/6 * dlnf0_dlnq / cc;
    dPsi[1] =   qke/3 * (Psi[0] - 2*Psi[2]);
    dPsi[2] =   qke/5 * (2 * Psi[1] - 3 * Psi[3]) - (h_prime/15 + 2*eta_prime/5) * dlnf0_dlnq / cc;
    for (int l=3; l<l_max; l++) {
        dPsi[l] = qke/(2*l+1.) * (l*Psi[l-1] - (l+1)*Psi[l+1]);
    }
    dPsi[l_max] = qke/(2*l_max+1.) * (l_max*Psi[l_max-1] - (l_max+1)*Psi_lmax_p1);

    return GSL_SUCCESS;
}

int evolve_gsl(double **Psi, double q, double k, int l_max, double tau_ini,
               double tau_final, double mass, double c_vel, double dlnf0_dlnq,
               double (*redshift_func)(double log_tau),
               double (*h_prime_func)(double k, double log_tau),
               double (*eta_prime_func)(double k, double log_tau),
               double tolerance) {

    /* Constants */
    op.params = malloc(6 * sizeof(double));
    op.params[0] = k;
    op.params[1] = q;
    op.params[2] = mass;
    op.params[3] = l_max;
    op.params[4] = dlnf0_dlnq;
    op.params[5] = c_vel;
    op.redshift_func = redshift_func;
    op.h_prime_func = h_prime_func;
    op.eta_prime_func = eta_prime_func;

    gsl_odeiv2_system sys = {func, NULL, l_max+1, &op};
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

    free(op.params);

    return status;
}
