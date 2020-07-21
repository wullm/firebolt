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

#ifndef BACKGROUND_H
#define BACKGROUND_H

#include "input.h"

struct background {
    double *z;
    double *log_tau;
    double **functions;
    char **titles;
    long int nrow;
    int ncol;
};

struct background_title_ids {
    /* Ids of various background functions f(tau) */
    int id_rho_crit, id_rho_tot, id_H;
    int id_rho_g, id_rho_b, id_rho_cdm, id_rho_ur;
    int id_rho_idr, id_rho_idm_dr, id_rho_scf, id_rho_dr;
    int id_p_scf, id_rho_dcdm;
    int id_p_tot;

    /* Whether these components are present in the table */
    char has_cdm, has_ur, has_idr, has_idm_dr, has_scf;
    char has_dr, has_dcdm;

    /* Number of ncdm species */
    int n_ncdm;

    /* Ids of the density and pressure functions for each ncdm species */
    int *id_rho_ncdm;
    int *id_p_ncdm;
};

enum background_format {
    Plain,
    CLASS
};

int readBackground(const struct params *pars, const struct units *us,
                   const struct cosmology_params *cosmo, struct background *bg);
int cleanBackground(struct background *bg);
int parseBackgroundTitles(const struct background *bg,
                          struct background_title_ids *bti);
int cleanBackgroundTitles(struct background_title_ids *bti);
#endif
