#include "cfg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// setup parameters
// ================
int LINELEN = 256;              // maximum length of line
char *SYNTAX = "%s = %s";       // string for parameter-value syntax
char COMMENT = '#';             // character for comments
int SETUPISDONE = 1;

// internal parameters
// ===================
FILE *cfgfp;
char *cfglp;
int nlines = 0;

// =============================================================================
// READY =======================================================================
// =============================================================================
//
// Get the config reader started:
// Use cfgAlloc(char *filename) to start
// Use cfgFree() when done

/* Allocate config reader with config file <filename>
 *
 * Parameters
 * ----------
 * char *filename : path to the config file to read
 *
 * Return
 * ------
 * int success : 0 if successful, 1 otherwise
 */
int cfgAlloc(char *filename)
{
    if (SETUPISDONE == 0)
    {
        printf("Cannot allocate config reader that is already allocated\n");
        cfgFree();
        exit(EXIT_FAILURE);
    }

    cfgfp = fopen(filename, "r");
    if (cfgfp == NULL)
    {
        printf("Could not open file: %s\n", filename);
        fclose(cfgfp);
        exit(EXIT_FAILURE);
    }

    cfglp = malloc(LINELEN * sizeof(char));
    if (cfglp == NULL)
    {
        printf("Failed to allocate line string");
        fclose(cfgfp);
        free(cfglp);
        exit(EXIT_FAILURE);
    }
    SETUPISDONE = 0;
    return EXIT_SUCCESS;
}

/* Free config reader
 *
 * Return
 * ------
 * int success : 0 if successful, 1 otherwise
 */
int cfgFree(void)
{
    if (SETUPISDONE == 1)
    {
        printf("Cannot free reader that is not allocated\n");
        exit(EXIT_FAILURE);
    }
    fclose(cfgfp);
    free(cfglp);
    SETUPISDONE = 1;
    return EXIT_SUCCESS;
}

// =============================================================================
// READ FILE ===================================================================
// =============================================================================

// general =====================================================================

/* Read line content with parameter key <name>
 *
 * Parameters
 * ----------
 * char *cont : pointer to save content to
 * char *name : name of the content's key
 *
 * Return
 * ------
 * int success : 0 if successful, 1 otherwise
 */
int getContent(char *cont, char *name)
{
    if (name == NULL || name == (char*)"")
    {
        printf("Cannot resolve parameter name: %s\n", name);
        return EXIT_FAILURE;
    }

    char key[LINELEN];
    while (fgets(cfglp, LINELEN, cfgfp) != NULL)
    {
        if (cfglp[0] == COMMENT || cfglp[0] == '\n')
            continue; // skip comments and empty lines

        if(sscanf(cfglp, SYNTAX, key, cont) == 2)
        {
            if (strcmp(key, name) == 0)
            {
                // line content for this key (name) is already saved in cont
                fseek(cfgfp, 0, SEEK_SET); // return to beginning of file
                return EXIT_SUCCESS;
            }
        }
    }
    // nothing found
    printf("Could not find key: %s\n", name);
    return EXIT_FAILURE;
}

// double ======================================================================

/* Return double value from config with key
 *
 * Parameters
 * ----------
 * char *name : name of the parameter key in config
 *
 * Return
 * ------
 * double result : value of parameters with key name
 */
double cfgGetd(char *name)
{
    char cont[LINELEN];
    if (getContent(cont, name) == 1)
    {
        cfgFree();
        exit(EXIT_FAILURE);
    }

    return atof(cont);
}

// int =========================================================================

/* Return int value from config with key
 *
 * Parameters
 * ----------
 * char *name : name of the parameter key in config
 *
 * Return
 * ------
 * int result : value of parameters with key name
 */
int cfgGeti(char *name)
{
    char cont[LINELEN];
    if (getContent(cont, name) == 1)
    {
        cfgFree();
        exit(EXIT_FAILURE);
    }

    return atoi(cont);
}
