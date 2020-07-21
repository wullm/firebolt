#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "../include/firebolt.h"

static inline void sucmsg(const char *msg) {
    printf("%s%s%s\n\n", TXT_GREEN, msg, TXT_RESET);
}

int main() {
    /* Read parameters */
    const char fname[] = "test_cosmology.ini";
    struct params pars;
    struct units us;
    struct cosmology_params cosmo;
    struct background bg;

    readParams(&pars, fname);
    readUnits(&us, fname);
    readCosmology(&cosmo, fname);

    /* Check if we can read the background data file */
    const char *background_fname = pars.BackgroundFile;
    FILE *f = fopen(background_fname, "r");
    assert(f != NULL);
    fclose(f);

    /* Read out the transfer functions */
    readBackground(&pars, &us, &cosmo, &bg);

    /* Verify that we have the expected data */
    assert(bg.nrow == 4621);
    assert(bg.ncol == 25);

    /* Verify some column titles */
    assert(strcmp(bg.titles[0], "z") == 0);
    assert(strcmp(bg.titles[1], "proper time [Gyr]") == 0);
    assert(strcmp(bg.titles[11], "(.)rho_ncdm[0]") == 0);
    assert(strcmp(bg.titles[23], "gr.fac. D") == 0);
    assert(strcmp(bg.titles[24], "gr.fac. f") == 0);

    /* Verify some data values */
    assert(fabs(bg.z[0] - 1.000000000000e+14)/bg.z[0] < 1e-3);
    assert(fabs(bg.z[4619] - 2.495252477297e-04)/bg.z[4619] < 1e-3);
    assert(fabs(bg.functions[2][0] - 7.84799257383e24)/bg.functions[2][0] < 1e-3);
    assert(fabs(bg.functions[2][4619] - 0.06905148423 )/bg.functions[2][4619] < 1e-3);
    // assert(fabs(bg.functions[2][0] - 2.561426030796e+22)/bg.functions[2][0] < 1e-3); //pre unit conversion
    // assert(fabs(bg.functions[2][4619] - 2.253700771351e-04)/bg.functions[2][4619] < 1e-3); //pre unit conversion
    assert(fabs(bg.functions[22][0] - 8.689880661743e-20)/bg.functions[22][0] < 1e-3);
    assert(fabs(bg.functions[3][0] - 1.389086530606e+04)/bg.functions[3][0] < 1e-3);
    assert(fabs(bg.functions[9][0] - 3.341e32)/bg.functions[9][0] < 1e-3);

    /* Clean up */
    cleanBackground(&bg);
    cleanParams(&pars);

    sucmsg("test_background:\t SUCCESS");
}
