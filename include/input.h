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

#ifndef INPUT_H
#define INPUT_H

#define DEFAULT_STRING_LENGTH 50

#define MPC_METRES 3.085677581282e22

/* The .ini parser library is minIni */
#include "../parser/minIni.h"

struct params {
    /* Simulation parameters */
    char *Name;
    char *TransferFunctionsFile;
    char *TransferFunctionsFormat;

    /* Output parameters */
    char *OutputDirectory;
};

struct units {
    double UnitLengthMetres;
    double UnitTimeSeconds;
    double UnitMassKilogram;

    /* Units for the transfer function input data */
    double TransferUnitLengthMetres;
    int Transfer_hExponent; //1 for h/Mpc; 0 for 1/Mpc
    int Transfer_kExponent; //0 for CLASS; -2 for CAMB/CMBFAST/Eisenstein-Hu
    int Transfer_Sign; //-1 for CLASS; +1 for CAMB/CMBFAST/Eisenstein-Hu
};

struct cosmology {
    double h;
};

int readParams(struct params *parser, const char *fname);
int readUnits(struct units *us, const char *fname);
int readCosmology(struct cosmology *cosmo, const char *fname);

#endif
