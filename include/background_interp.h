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

#ifndef BACKGROUND_INTERP_H
#define BACKGROUND_INTERP_H

#include "background.h"

int bg_interp_init(const struct background *bg);
int bg_interp_switch_func(const struct background *bg, int index_func);
int bg_interp_free(const struct background *bg);
double bg_func_at_log_tau(double log_tau);
double bg_z_at_log_tau(double log_tau);

#endif
