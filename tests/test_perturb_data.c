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
    struct perturb_data ptdat;

    readParams(&pars, fname);
    readUnits(&us, fname);
    readCosmology(&cosmo, fname);

    /* Check if we can read the background data file */
    const char *background_fname = pars.PerturbFile;
    FILE *f = fopen(background_fname, "r");
    assert(f != NULL);
    fclose(f);

    /* Read the perturb data */
    readPerturb(&pars, &us, &ptdat);

    /* Header checks */
    assert(ptdat.k_size == 616);
    assert(ptdat.tau_size == 736);
    assert(ptdat.n_functions == 6);

    /* Check some titles */
    assert(strcmp(ptdat.titles[0], "h_prime") == 0);
    assert(strcmp(ptdat.titles[1], "eta_prime") == 0);
    assert(strcmp(ptdat.titles[2], "H_T_Nb_prime") == 0);
    assert(strcmp(ptdat.titles[3], "t_tot") == 0);
    assert(strcmp(ptdat.titles[4], "d_ncdm[0]") == 0);
    assert(strcmp(ptdat.titles[5], "t_ncdm[0]") == 0);

    /* Check some data points */
    assert((ptdat.k[0] - 4.82858e-06)/ptdat.k[0] < 1e-5);
    assert((ptdat.k[615] - 5.04006)/ptdat.k[615] < 1e-5);
    assert((ptdat.log_tau[0] - -0.295592)/ptdat.log_tau[0] < 1e-5);
    assert((ptdat.log_tau[735] - 3.82635)/ptdat.log_tau[735] < 1e-5);
    assert((ptdat.delta[0] - -352.536)/ptdat.delta[0] < 1e-5);
    assert((ptdat.delta[2720255] - -9.77855e-05)/ptdat.delta[2720255] < 1e-5);

    /* Clean up */
    cleanPerturb(&ptdat);
    cleanParams(&pars);

    sucmsg("test_perturb_data:\t SUCCESS");
}
