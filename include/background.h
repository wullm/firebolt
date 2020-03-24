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
    double **functions;
    char **titles;
    long int nrow;
    int ncol;
};

enum background_format {
    Plain,
    CLASS
};

int readBackground(const struct params *pars, const struct units *us,
                   const struct cosmology *cosmo, struct background *bg);
int cleanBackground(struct background *bg);

#endif
