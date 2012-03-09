
#ifndef initialize
#define initialize

#include <stdio.h>

#include "/usr/include/gsl/gsl_linalg.h"

//Real Simulation
#define sample_length 10000000
//#define sample_length 1000000
//Test Mode
//#define sample_length 10000
//#define sample_length 100000
//#define sample_length 80000

#define PI 3.1415926535897932384626433832795
#define OK 1


int var_init(int *var,int k);
int var_init(int **var,int row,int col);

int var_init(long double *var,int k);

int var_init(long double **var,int row, int col);
int var_init(long double ***var,int x, int y, int z);

int var_init(char *var,int k);

#endif
