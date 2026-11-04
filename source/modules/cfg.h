#ifndef CFG_H
#define CFG_H

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
int cfgAlloc(char *filename);

/* Free config reader
 *
 * Return
 * ------
 * int success : 0 if successful, 1 otherwise
 */
int cfgFree(void);

// =============================================================================
// READ FILE ===================================================================
// =============================================================================

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
double cfgGetd(char *name);

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
int cfgGeti(char *name);

#endif
