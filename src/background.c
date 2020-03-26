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
#include <ctype.h> //isdigit
#include <math.h>
#include "../include/background.h"

static inline char is_comment_line(char *line) {
    /* Lines starting with '#' or that don't begin with a number */
    float tmp;
    return line[0] == '#' || (sscanf(line, "%e", &tmp) <= 0);
}

void countRowsCols(FILE *f, int *nrow, int *ncol, int *leading_comment_lines) {
    char line[2000];

    /* Local counters */
    int rows = 0;
    int cols = 0;
    int comments = 0;

    /* Count the number of rows with data */
    while (fgets(line, sizeof(line), f)) {
        if (is_comment_line(line)) continue; /* skip comments */
        rows++;
    }

    /* Go back to the start */
    rewind(f);

    /* Count the number of leading comment lines */
    while (fgets(line, sizeof(line), f) && is_comment_line(line)) {
        comments++;
    }

    /* On the first non-comment line, count the number of columns */
    int read = 0, bytes;
    float tmp;
    while(sscanf(line + read, "%e%n", &tmp, &bytes) > 0) {
        read += bytes;
        cols += 1;
    }

    /* Update the counters */
    *nrow = rows;
    *ncol = cols;
    *leading_comment_lines = comments;
}


void readTitles(char *line, enum background_format format, char **titles) {
    if (format == Plain) {
        /* The Plain format is "# z name name name " */
        /* We want to read out the column names */

        /* Skip the first character, which should be '#' */
        char title[50];
        int j = 0, read = 0, bytes;

        /* Read the column titles & count the # of bytes read */
        while(sscanf(line + read, "%s%n", title, &bytes) > 0) {
            titles[j] = malloc(strlen(title) + 1);
            strcpy(titles[j], title);

            read += bytes;
            j++;
        }
    } else if (format == CLASS) {
        /* The CLASS format is "# 1:z    2:name     3:name    4:name ...." */
        /* We want to read out the column names */

        /* Skip the first character, which should be '#' */
        char string[50];
        int j = 0, read = 1, bytes;

        /* Read string until encountering a colon & count the # of bytes read */
        while(sscanf(line + read, "%[^:]%n", string, &bytes) > 0) {
            /* Parse the title from the full string */
            int string_read = 0, string_bytes; //second counter within string
            char part[50], title[50];
            strcpy(title, "");
            /* Read words (separated by spaces) until encountering a number */
            while (sscanf(string + string_read, "%s%n", part, &string_bytes) > 0
                   && !isdigit(part[0])) {
                strcat(title, part);
                strcat(title, " ");
                string_read += string_bytes;
            }

            /* If we found a non-empty title, store it */
            if (strlen(title) > 0) {
                title[strlen(title)-1] = '\0'; /* delete trailing whitespace */
                titles[j] = malloc(strlen(title) + 1);
                strcpy(titles[j], title);
                j++;
            }

            read += bytes + 1;
        }
    }
}

int readBackground(const struct params *pars, const struct units *us,
                   const struct cosmology *cosmo, struct background *bg) {
    const char *fname = pars->BackgroundFile;
    const char *formatString = pars->BackgroundFormat;
    enum background_format format;

    /* Parse the expected format of the background file */
    if (strcmp(formatString, "Plain") == 0) {
        format = Plain;
    } else if (strcmp(formatString, "CLASS") == 0) {
        format = CLASS;
    } else {
        printf("ERROR: Unknown background file format.\n");
        return 1;
    }

    printf("Reading the background data from '%s'.\n", pars->BackgroundFile);

    /* Open the data file */
    FILE *f = fopen(fname, "r");
    char line[2000];

    /* Determine the size of the table */
    int nrow;
    int ncol;
    int leading_comment_lines;
    countRowsCols(f, &nrow, &ncol, &leading_comment_lines);

    bg->nrow = nrow;
    bg->ncol = ncol;

    /* Move to the last leading comment line */
    rewind(f);
    for (int i=0; i<leading_comment_lines; i++) {
        fgets(line, sizeof(line), f);
    }

    /* Allocate memory for pointers to the column title strings */
    bg->titles = malloc(ncol * sizeof(char*));

    /* Read out the column titles */
    readTitles(line, format, bg->titles);

    /* Allocate memory for the background data */
    bg->z = malloc(nrow * sizeof(double));
    bg->log_tau = malloc(nrow * sizeof(double));
    bg->functions = malloc(ncol * sizeof(double*));
    for (int i=0; i<ncol; i++) {
        bg->functions[i] = malloc(nrow * sizeof(double));
    }

    /* Read out the data into the arrays */
    int row = 0;
    double number = 0;
    while (fgets(line, sizeof(line), f)) {
        if (is_comment_line(line)) continue; /* skip comments */
        int read = 0, bytes;
        int col = 0;
        while(sscanf(line + read, "%le%n", &number, &bytes) > 0) {
            if (col == 0) {
                bg->z[row] = number;
            } else {
                bg->functions[col-1][row] = number;
            }
            read += bytes;
            col++;
        }
        row++;
    }

    fclose(f);

    /* Unit conversions.
     */

    /* CLASS to internal units conversion factor */
    const double mpc_length_factor = MPC_METRES / us->UnitLengthMetres;
    const double mpc_time_factor = mpc_length_factor / us->SpeedOfLight;
    const double gyr_time_factor = GYR_OVER_MPC * mpc_time_factor;
    const double density_factor = 8 * M_PI * us->GravityG / 3;

    for (int i=0; i<nrow; i++) {
        /* Convert proper times to internal units [Gyr in CLASS] */
        bg->functions[0][i] *= gyr_time_factor;
        /* Convert the conformal times to internal units [Mpc in CLASS] */
        bg->functions[1][i] *= mpc_time_factor;
        /* Convert Hubble constant to internal units [1/Mpc in CLASS] */
        bg->functions[2][i] /= mpc_time_factor;
        /* Convert distances to internal units [all Mpc in CLASS] */
        bg->functions[3][i] *= mpc_length_factor;
        bg->functions[4][i] *= mpc_length_factor;
        bg->functions[5][i] *= mpc_length_factor;
        bg->functions[6][i] *= mpc_length_factor;
        /* Convert all densities and pressures to internal units [all Mpc^-2].
         * Note these quantities are multiplied by 8piG/3 and have dimensions
         * of inverse time squared.
         */
        for (int j=7; j<ncol-3; j++) {
            /* Convert to densities in internal units (MASS/LENGTH^3) */
            bg->functions[j][i] /= density_factor * mpc_time_factor * mpc_time_factor;
        }
        /* These last two columns are dimensionless growth factors */
    }

    /* Create the log-tau table */
    for (int i=0; i<nrow; i++) {
        bg->log_tau[i] = log(bg->functions[1][i]);
    }

    return 0;
}


int cleanBackground(struct background *bg) {
    free(bg->titles);
    free(bg->functions);
    free(bg->z);

    return 0;
}
char has_cdm, has_ur, has_idr, has_idm_dr, has_scf;
char has_dr, has_dcdm;


int parseBackgroundTitles(const struct background *bg,
                          struct background_title_ids *bti) {

    /* Reset the indicators */
    bti->has_cdm = 0;
    bti->has_ur = 0;
    bti->has_idr = 0;
    bti->has_idm_dr = 0;
    bti->has_scf = 0;
    bti->has_dr = 0;
    bti->has_dcdm = 0;

    /* Parse the background function titles for known ids  and
     * count the number ncdm species among the columns. */
    bti->n_ncdm = 0;
    for (int i=0; i<bg->ncol; i++) {
        char word[50];
        sscanf(bg->titles[i], "%s", word);

        /* Note, we assign id = i-1, because the first column is
         * the redshift column, which is treated separately. */

        if (strcmp(word, "H") == 0) {
            bti->id_H = i-1;
        } else if (strcmp(word, "(.)rho_crit") == 0) {
            bti->id_rho_crit = i-1;
        } else if (strcmp(word, "(.)rho_tot") == 0) {
            bti->id_rho_tot = i-1;
        } else if (strcmp(word, "(.)rho_g") == 0) {
            bti->id_rho_g = i-1;
        } else if (strcmp(word, "(.)rho_b") == 0) {
            bti->id_rho_b = i-1;
        } else if (strcmp(word, "(.)rho_cdm") == 0) {
            bti->has_cdm = 1;
            bti->id_rho_cdm = i-1;
        } else if (strcmp(word, "(.)rho_ur") == 0) {
            bti->has_ur = 1;
            bti->id_rho_ur = i-1;
        } else if (strcmp(word, "(.)rho_idr") == 0) {
            bti->has_idr = 1;
            bti->id_rho_idr = i-1;
        } else if (strcmp(word, "(.)rho_idm_dr") == 0) {
            bti->has_idm_dr = 1;
            bti->id_rho_idm_dr = i-1;
        } else if (strcmp(word, "(.)rho_scf") == 0) {
            bti->has_scf = 1;
            bti->id_rho_scf = i-1;
        } else if (strcmp(word, "(.)rho_dr") == 0) {
            bti->has_dr = 1;
            bti->id_rho_dr = i-1;
        } else if (strcmp(word, "(.)rho_dcdm") == 0) {
            bti->has_dcdm = 1;
            bti->id_rho_dcdm = i-1;
        } else if (strcmp(word, "(.)p_scf") == 0) {
            bti->has_scf = 1;
            bti->id_p_scf = i-1;
        } else if (strcmp(word, "(.)p_tot") == 0) {
            bti->id_p_tot = i-1;
        }

        /* Look for ncdm density columns (ids are assigned later) */
        int ncdm_num;
        if (sscanf(bg->titles[i], "(.)rho_ncdm[%d]", &ncdm_num) > 0) {
            bti->n_ncdm++;
        }
    }

    /* Scan for ncdm columns */
    bti->id_rho_ncdm = malloc(bti->n_ncdm * sizeof(int));
    bti->id_p_ncdm = malloc(bti->n_ncdm * sizeof(int));
    for (int i=0; i<bg->ncol; i++) {
        char word[50];
        sscanf(bg->titles[i], "%s", word);

        for (int j = 0; j<bti->n_ncdm; j++) {
            char search_rho[50];
            char search_p[50];
            sprintf(search_rho, "(.)rho_ncdm[%d]", j);
            sprintf(search_p, "(.)p_ncdm[%d]", j);

            if (strcmp(word, search_rho) == 0) {
                bti->id_rho_ncdm[j] = i-1;
            } else if (strcmp(word, search_p) == 0) {
                bti->id_p_ncdm[j] = i-1;
            }
        }
    }

    return 0;
}

int cleanBackgroundTitles(struct background_title_ids *bti) {
    free(bti->id_rho_ncdm);
    free(bti->id_p_ncdm);
    return 0;
}
