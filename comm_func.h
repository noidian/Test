#ifndef COMM_FUNC
#define COMM_FUNC

#include "init.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "math_func.h"
#include "generic_func.h"
#include "impulse.h"
//#include "/usr/include/gsl/gsl_matrix_long_double.h"
//#include "/usr/include/gsl/gsl_vector_long_double.h"


//Assume length2<length1 always and output = length1 - (2*(length2-1))
int convolution(long double *signal1,long double *signal2,int length1,int length2,long double *output);

int convolution(int *signal1,long double *signal2,int length1,int length2,long double *output);

int awgn(long double *output,long double *input,int Nb_samples,float SNR);

int input_ITI_adder(long double *output,int *adj_OD,long double *filter_n1,int *main,long double *filter_0, int *adj_ID, long double *filter_1,int Nb_samples,int Nb_filter);

int lagrange_err_min_1D(long double *equalizer, int no_equalizer_coeff, long double *target, int no_of_target_vals, int *a, int numberofvalues, long double *matched_main);

int lagrange_err_min_2D(long double *f_0,long double *f_n1,long double *f_1,int Nb_eq_0,int Nb_eq_n1,int Nb_eq_1,long double *G,int Nb_target,int **a_merge,int Nb_samples,long double *matched_main,long double *matched_adj1_OD,long double *matched_adj1_ID);

int lagrange_err_min_onesided(long double *f_0,long double *f_n1,int Nb_eq_0,int Nb_eq_n1,long double *G,int Nb_target,int **a_merge,int Nb_samples,long double *matched_main,long double *matched_adj1_OD);

int gene_viterbi_matrix (const int Polynom1, const int Polynom2,int **inputmat, int **nextstate, int **outmat);

int int2bin (const int number, int table[], const int length_table);

int shift_reg (int reg[], const int length_reg, const int new_value);

int bin2int (int *number, const int table[], const int length_table);

int encode_bit (const int Polynom1, const int Polynom2,const int shift_register[], const int length_register, double *out_poly1, double *out_poly2, const char flag);

int va_detection(int MemoryDepth, int Nbstates, long double *target, int target_length, long double *input, int *output,int ip_length);

int sova_detection(int MemoryDepth, int Nbstates, long double *target, int target_length,long double *input, int *output,long double *llr,long int *llr_length,int ip_length);
//int sova_detection(int MemoryDepth, int NbStates, long double *input, long double *output, int **inputmat, int **nextstate, int **outmat);

int map_detection(int MemoryDepth, int Nbstates, long double *target, int target_length, long double *input, int *output,long double *llr,long int *llr_length,int ip_length,float SNR);

long double ber_compute(int *input1,int *input2,int Nb_values);

long double oneD_equalize_oneD_SOVA(int *tx_bits,int current_sample_length, long double *rx_bits,int Nb_eq,int Nb_target,int *ber_bits,long int *llr_length, long double *llr,float SNR);

long double twoD_equalize_oneD_SOVA(int **tx_merge,int current_sample_length, long double *rx_main,long double *rx_adj1_OD,long double *rx_adj1_ID,int Nb_eq_main,int Nb_eq_adj1_OD,int Nb_eq_adj1_ID,int Nb_target,int *ber_bits,long int *llr_length, long double *llr,float SNR);

long double onesided_equalize_oneD_SOVA(int **tx_merge,int current_sample_length, long double *rx_main,long double *rx_adj1_OD,int Nb_eq_main,int Nb_eq_adj1_OD,int Nb_target,int *ber_bits,long int *llr_length, long double *llr,float SNR);



#endif
