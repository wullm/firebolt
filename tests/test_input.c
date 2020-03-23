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
    /* Test reading parameters */
    const char fname[] = "test_cosmology.ini";
    struct params pars;

    readParams(&pars, fname);

    assert(strcmp(pars.Name, "Test Simulation") == 0);
    assert(strcmp(pars.OutputDirectory, "../tests") == 0);
    assert(strcmp(pars.BackgroundFile, "../background/class_example_background.dat") == 0);
    assert(strcmp(pars.BackgroundFormat, "CLASS") == 0);
    assert(strcmp(pars.PerturbFile, "test_perturb.hdf5") == 0);

    /* Test reading units */
    struct units us;
    readUnits(&us, fname);

    assert(us.UnitLengthMetres == 3.086e22);
    assert(us.UnitTimeSeconds == 3.154e16);
    assert(us.UnitMassKilogram == 1.989e40);
    assert(us.BackgroundUnitLengthMetres == 3.085677581282e22);

    /* Test reading cosmology */
    struct cosmology cosmo;
    readCosmology(&cosmo, fname);

    assert(cosmo.h == 0.67556);

    /* Clean up */
    cleanParams(&pars);

    sucmsg("test_input:\t SUCCESS");
}
