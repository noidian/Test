#ifndef GENERIC_FUNC
#define GENERIC_FUNC

#include <stdio.h>
#include <stdlib.h>
#include "init.h"


int fileopen(char *filename,int *input,int numberofvalues);

int fileopen(char *filename,long double *input,int numberofvalues);

int array_copy (long double *val1,long double *val, int numberofvalues);

int array_copy(long double *val1,long double *val, int start_val, int end_val);

int twoD_array_copy_int (int **val1,int *val, int row,int col);

int twoD_array_copy_ldouble (long double **val1,long double *val, int row,int col);

int vector_copy (int *val1,int *val, int numberofvalues);

int array_resize(long double *input,int start_val,int end_val);

int array_resize(int *input,int start_val,int end_val);

int shiftback_wd (char *cwd,int length,char *opwd,int shiftback_number);

int dot2underscore(char *str,int length);

int dot2null(char *str, int length);

long double ** matrix_calloc(int row, int col);

long double ***matrix_calloc(int x, int y,int z);

int **matrix_calloc_int(int row, int col);

long double *vector_calloc(int col);

int *vector_calloc_int(int col);

char *vector_calloc_char(int col);

int matrix_free(long double **input,int row);

int matrix_free(long double ***input,int row,int col);

int matrix_free(int **input,int row);

int filestore_ber(long double SNR,long double ber, char *filename);

int filestore_ber(long double ber, char *filename);

int filestore_prob(long double **input, char *filename, int rows);
#endif
