#ifndef MATH_FUNC
#define MATH_FUNC

# include <math.h>
#include "generic_func.h"
#include "init.h"

long double norm(long double *val, int numberofvalues);

long double norm(int *val, int numberofvalues);

long double norm(int **val, int rows, int cols);

long double norm(long double **val, int rows, int cols);

int x_corr2(int **input1,int rows1, int cols1,int **input2, int rows2, int cols2, int rshift,int cshift,int Nb_samples, long double **output);

int x_corr2(int **input1,int rows1, int cols1,long double **input2, int rows2, int cols2, int rshift,int cshift,int Nb_samples, long double **output);

int x_corr2(long double **input1,int rows1, int cols1,long double **input2, int rows2, int cols2, int rshift,int cshift,int Nb_samples, long double **output);

int  x_corr(int *x,int *y,int cshift, int numberofvals, long double *x_corr_op);

int  x_corr(int *x,long double *y,int cshift, int numberofvals, long double *x_corr_op);

int  x_corr(long double *x,long double *y,int cshift, int numberofvals, long double *x_corr_op);

long double  min(long double input1,long double input2);

int array_addition (long double *output,long double *input1,long double *input2, int numberofvalues);

int array_divide (long double *val,long double div_val, int numberofvalues);

int array_divide (long double **val,long double div_val, int rows,int numberofvalues);

int array_multiply (long double *val,long double mult_val, int numberofvalues);

int matrix_transpose(long double **input_transpose,int rows,int cols, long double **input);

int matrix_multiply (long double **output,long double **input1,int row1, int col1,long double **input2,int row2,int col2);

int matrix_multiply (long double *output,long double *input1,int col1,long double **input2,int row2,int col2);

int matrix_multiply (long double *output,long double **input1,int row1,int col1,long double *input2,int row2);

int  matrix_subtraction(long double **output,long double **input1,long double **input2,int row,int col);

int matrix_inverse(long double **output,int row,int col,long double **input);
#endif
