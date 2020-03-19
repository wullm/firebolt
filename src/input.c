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
#include "../include/input.h"

int readParams(struct params *pars, const char *fname) {
     // pars->GridSize = ini_getl("Box", "GridSize", 64, fname);
     // pars->BoxLen = ini_getd("Box", "BoxLen", 1.0, fname);
     // pars->Homogeneous = ini_getbool("Simulation", "Homogeneous", 0, fname);

     /* Read strings */
     int len = DEFAULT_STRING_LENGTH;
     pars->OutputDirectory = malloc(len);
     pars->Name = malloc(len);
     pars->TransferFunctionsFile = malloc(len);
     pars->TransferFunctionsFormat = malloc(len);
     ini_gets("Output", "Directory", "./output", pars->OutputDirectory, len, fname);
     ini_gets("Simulation", "Name", "No Name", pars->Name, len, fname);
     ini_gets("TransferFunctions", "File", "", pars->TransferFunctionsFile, len, fname);
     ini_gets("TransferFunctions", "Format", "Plain", pars->TransferFunctionsFormat, len, fname);

     return 0;
}

int readUnits(struct units *us, const char *fname) {
    /* Internal units */
    us->UnitLengthMetres = ini_getd("Units", "UnitLengthMetres", 1.0, fname);
    us->UnitTimeSeconds = ini_getd("Units", "UnitTimeSeconds", 1.0, fname);
    us->UnitMassKilogram = ini_getd("Units", "UnitMassKilogram", 1.0, fname);

    /* Get the transfer functions format */
    char format[DEFAULT_STRING_LENGTH];
    ini_gets("TransferFunctions", "Format", "Plain", format, DEFAULT_STRING_LENGTH, fname);

    /* Format of the transfer functions */
    int default_h_exponent;
    int default_k_exponent;
    int default_sign;
    if (strcmp(format, "CLASS") == 0) {
        default_h_exponent = 1;
        default_k_exponent = 0;
        default_sign = -1;
    } else {
        default_h_exponent = 0;
        default_k_exponent = -2;
        default_sign = +1;
    }
    us->TransferUnitLengthMetres = ini_getd("TransferFunctions", "UnitLengthMetres", MPC_METRES, fname);
    us->Transfer_hExponent = ini_getl("TransferFunctions", "hExponent", default_h_exponent, fname);
    us->Transfer_kExponent = ini_getl("TransferFunctions", "kExponent", default_k_exponent, fname);
    us->Transfer_Sign = ini_getl("TransferFunctions", "Sign", default_sign, fname);

    return 0;
}

int readCosmology(struct cosmology *cosmo, const char *fname) {
     cosmo->h = ini_getd("Cosmology", "h", 0.70, fname);

     return 0;
}
