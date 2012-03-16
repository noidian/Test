#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "LDPC_mod2.h"

//20100628pm:

//LDPC sum-product algorithm:

int LDPC_mod2::Sum_Product_Algorithm( mod2sparse *hh,
									 char *encoded,
				                     char *decoded //output
									 )
{
//input: llr_fn
//input: H

//output:



	int i,j;
	int M,N;
	float fsum;
	int signal;
	int be;
	int it;

	mod2entry *e;
	mod2entry *et;

	M=mod2sparse_rows(hh);
	N=mod2sparse_cols(hh);


	//1: initial SPA:
	//initial e->pr with llr_fn:
	for(i=0;i<N;i++)
	{
		for(e=mod2sparse_first_in_col(hh,i); !mod2sparse_at_end(e); e=mod2sparse_next_in_col(e))
		{//for all entries in a column:
			e->pr = llr_fn[i];
		}
	}//for i<N.

	//2:
	for(it=0;it<max_iter;it++)
	{
		//2.1: update all variables in a check:
		// i.e. bit-to-check process:
		// update e->lr from e-> pr:
		for(j=0;j<M;j++)
		{//for all checks:
			for(e=mod2sparse_first_in_row(hh,j);!mod2sparse_at_end(e);e=mod2sparse_next_in_row(e))
			{
				fsum=0;
				signal=1;
				//lr-j = fel( sum(fel(i))  ).
				for(et=mod2sparse_first_in_row(hh,j);!mod2sparse_at_end(et);et=mod2sparse_next_in_row(et))
				{
					if(et!=e)
						fsum +=  fel_func( absolute_value_f(et ->pr) );

					if( et->pr<0)
						signal = -1*signal;
				}
				fsum = fel_func(fsum);
				fsum = fsum*signal;
				e->lr = fsum; //store this message to another memory, in case it is wrongly updated in the same loop if in e->pr.
			}//for e in a row.
		}//for j.

		//2.2: check-to-bit process:
		be=0;
		for(i=0;i<N;i++)
		{
			fsum=0;
			for(e=mod2sparse_first_in_col(hh,i);!mod2sparse_at_end(e);e=mod2sparse_next_in_col(e) )
				fsum += e->lr;

			for(e=mod2sparse_first_in_col(hh,i);!mod2sparse_at_end(e);e=mod2sparse_next_in_col(e) )
				e->pr = fsum - e->lr;

			//make hard-decision on fsum:
			if(fsum >0)
				decoded[i] = 1;
			else
				decoded[i] = 0;

			if(decoded[i]!=encoded[i])
				be++;
		}//for i.

	}//for it.

return be;
}

double LDPC_mod2::fel_func(double x)
{
	//Need to normalize. If input x is larger than 31, the result is extremely close to 0 and can not be corrected stored.
	if(x>fel_func_max_input)
		x=fel_func_max_input;

//	return (log( (exp(x)+1)/(exp(x)-1) )/0.99999932734728);	//This is the precise form, but close to the following expression.
	return log( (exp(x)+1)/(exp(x)-1) );
}


//-------------
void LDPC_mod2::read_binary_generator_matrix_Gb(char *filename, int row_gb , int col_gb)
{
//read in generator matrix:
	int i,j;
	int it;

	FILE *fp;

	fp=fopen(filename,"r");
	if(fp==0)
	{
		printf("Wrong. %s.\n exit. \n",filename);
		exit(1);
	}

	if(flag_G_binary_allocated==0)
	{//if 0, G_binary is not allocated:
		G_binary = (char **)calloc(row_gb,sizeof(char *));
		for(i=0;i<row_gb;i++)
			G_binary[i]=(char *)calloc(col_gb,sizeof(char));
	}
	flag_G_binary_allocated=1;

	//
	for(j=0;j<row_gb;j++)
	{
		for(i=0;i<col_gb;i++)
		{
			fscanf(fp,"%1d",&it);
			G_binary[j][i]=it;
		}
	}
	fclose(fp);

}//

void LDPC_mod2::encoding_with_systemic_binary_G(char * source, char * encoded)
{
//input: G_bianry, Gb_rows, Gb_cols;

	int i,j;

	if(flag_G_binary_allocated==0)
	{
		printf("The binary generator matrix is not allocated. Exit.\n");
		exit(0);
	}//

	int it;
	for(j=0;j<Gb_rows;j++)
	{
		it=0;
		for(i=0;i<Gb_cols;i++)
			it = (it + source[i]*G_binary[j][i])%2;

		encoded[j]=it;
	}//for j.

return;
}

void LDPC_mod2::check_binary_G_for_mapping_infor_to_code()
{
//this function is only used for systemic matrix.
//input: G_binary, Gb_rows, Gb_cols;
//mapping_infor_to_code[];
	int i,j;


	//
	printf("Now to run check_binary_G_for_mapping_infor() for systemic LDPC.\n");

	if(flag_mapping_infor_to_code_allocated==0) //if 0, mapping_infor_to_code[] is not allocated.
		mapping_infor_to_code = (int *)calloc(Gb_cols,sizeof(int));
	flag_mapping_infor_to_code_allocated=1;

	//
	int it,fs;

	for(i=0;i<Gb_cols;i++)
		mapping_infor_to_code[i]=-1;

	for(j=0;j<Gb_rows;j++)
	{
		it=0;
		for(i=0;i<Gb_cols;i++)
		{
			if(G_binary[j][i]==1)
			{
				it++;
				fs=i;
			}
			if(it>1)
				break;
		}
		if(it==1)
		{
			mapping_infor_to_code[fs]=j;
		}
	}//for j.

	//check whether all infor bits are mapped in Gb_binary:
	for(i=0;i<Gb_cols;i++)
	{
		if(mapping_infor_to_code[i]==-1)
		{
			printf("Wrong. The %d-th information bit is not mapped in G_binary.\n");
			free(mapping_infor_to_code);
			exit(0);
		}
	}//for i.
return ;
}


void LDPC_mod2::read_binary_permuted_H_for_RU_encoding(char *filename_ht, int row_count, int col_count)
{
//read in permuted_H for RU encoding. The left upper part should be T
//output: permuted_H
	int i,j;
	int it;

	FILE *fp;

	fp=fopen(filename_ht,"r");
	if(fp==0)
	{
		printf("Wrong. %s does not exist.\n",filename_ht);
		exit(0);
	}

	if(flag_permuted_H_allocated==0)
	{
		permuted_H=(char**)calloc(size_T,sizeof(char*));
		for(j=0;j<row_count;j++)
			permuted_H[j]=(char *)calloc(col_count,sizeof(char));
		flag_permuted_H_allocated=1;
	}//if.

	//
	for(j=0;j<row_count;j++)
	{
		for(i=0;i<col_count;i++)
		{
			fscanf(fp,"%1d",&it);
			permuted_H[j][i]=it;
		}
	}//for j.
	fclose(fp);

}

void LDPC_mod2::read_column_swap_tracker(char *filename, int col_count)
{//read ColumnSwapTracker[]:
	int i,it;
	FILE *fp;

	fp=fopen(filename,"r");
	if(fp==0)
	{
		printf("Wrong. %s not exist. exit.\n");
		exit(0);
	}

	ColumnSwapTracker = (int *)calloc(col_count,sizeof(int));
	flag_columnswaptracke_allocated=1;
	for(i=0;i<col_count;i++)
	{
		fscanf(fp,"%d",&it);
		ColumnSwapTracker[i]=it-1;
	}
return;
}

void LDPC_mod2::encoding_in_RU_method_through_permuted_H_with_g_as_0(char * source, char*encoded)
{
//encoding through RU method, and this function only works for g=0 case.
//encoding through RU method:
//input: permuted_H;
//input: size_T, size_g,info_length;
//input: Ht_rows.
	int i,j;

	//Note: currently, this function only works for size_g=0!!

	if(size_g!=0)
	{
		printf("currently, this function only works for size_g=0!! \nExit \n");
		exit(0);
	}

	//information bits is [size_T+size_g ~ to the end]:
	char * temp_x;
	temp_x=(char*)calloc(Ht_rows,sizeof(char));

	for(j=0;j<Ht_rows;j++)
		temp_x[j] = 0;

	for(i=0;i<info_length;i++)
	{
		for(j=0;j<Ht_rows;j++)
		{
			temp_x[i] = (temp_x[i]+ permuted_H[j][i+Ht_rows]*source[i])%2;
		}//for j.
	}//for i.

	//then, to get the redundant bits through back_substitution:
	//the redundant bits are encoded[0 - Ht_rows-1]:
	for(i=Ht_rows;i<Ht_cols;i++)
		encoded[i]=source[i-Ht_rows];

	for(i=Ht_rows-1;i>=0;i--)
	{
		encoded[i]=temp_x[i];
		for(j=Ht_rows-1;j>i;j--)
			encoded[i]= (encoded[i]+permuted_H[j][i]*encoded[i])%2;
	}

	//
	free(temp_x);
	return;
}






void LDPC_mod2::trying_getting_trapping_set( mod2sparse *H, int th)
{


	int i,j;

	char *flag_chosen_error_bit;	//the flag if one bit is chosen to be incorrect.
	int count_of_chosen_bits;
	char *flag_check;				//un-satisfied checks;
	char *usc_of_bit;		//the number of un-satisfied checks for each bit.
	float bc_ratio; //the threshold ratio of error bits to error checks; if the actual ratio is larger than this, then this might be a trapping set.
	int cbc; //the count of correct bits;
	int uscc; //count of un-satisfied checks.
	int max_trapping_set; //the maximum number of bits in a trapping set.
	int flag_find_a_trapping_set;
	int it, max_usc_of_bit,cs;
	char first_bit;
	mod2entry *e;
	mod2entry *f;
	int col_h,row_h;

	col_h=H->n_cols;
	row_h=H->n_rows;

	flag_find_a_trapping_set=0;
	max_trapping_set = th; //42 for 142*255.

	bc_ratio = 1;
	//
//	bc_ratio = (max_col_weight+1)/2;

	//--
	flag_chosen_error_bit=(char *)calloc(H->n_cols,sizeof(char));
	usc_of_bit=(char *)calloc( H->n_cols,sizeof(char));
	flag_check=(char *)calloc( H->n_rows,sizeof(char));


	//
	while(flag_find_a_trapping_set==0)
	{
//		printf("Try... ");
		//initial:
		count_of_chosen_bits = 0;
		uscc=0;
		first_bit=1;
		for(i=0;i<col_h;i++)
			flag_chosen_error_bit[i]=0;
		for(i=0;i<col_h;i++)
			usc_of_bit[i]=0;
		for(i=0;i<row_h;i++)
			flag_check[i]=0;

		while(first_bit==1 || (count_of_chosen_bits<max_trapping_set && uscc*1.0 > count_of_chosen_bits*bc_ratio) )
		{
			if(first_bit==1)
			{
				cs=rand()%col_h;
				first_bit=0;
				flag_chosen_error_bit[cs]=1;
				count_of_chosen_bits++;
			}
			else
			{
				//--find the number of un-chosen bits with max_usc_of_bit un-satisfied checks:
				it=0;
				for(i=0;i<col_h;i++)
				{
					if(flag_chosen_error_bit[i]==0 && usc_of_bit[i]==max_usc_of_bit)
						it++;
				}//for i.
				cs=(rand())%it; //randomly choose one bit.
				count_of_chosen_bits++;

				//-- find the cs-th max bit:
				it=0;
				for(i=0;i<col_h;i++)
				{
					if(flag_chosen_error_bit[i]==0 && usc_of_bit[i]==max_usc_of_bit)
					{
						if(it==cs)
						{
							cs=i;
							flag_chosen_error_bit[cs]=1;
							break;
						}
						else
							it++;
					}//if.
				}//for i.

			}//else.

			//-- update un-satisfied checks:
			for (e = mod2sparse_first_in_col(H,cs); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
				flag_check[e->row] = (flag_check[e->row]+1)%2;

				//-- count un-satisfied checks:
				if(flag_check[e->row]==1)
					uscc++;
				else
					uscc--;

				//--
				for(f=mod2sparse_first_in_row(H,e->row); !mod2sparse_at_end(f); f=mod2sparse_next_in_row(f))
				{
					if(flag_check[e->row]==1)
						usc_of_bit[f->col]++;
					else
						usc_of_bit[f->col]--;

					if(usc_of_bit[f->col]<0)
						exit(0);
				}//for f.
			}

			//-- update the number of error checks that each bit partipates in:
			max_usc_of_bit=0;
			for(i=0;i<col_h;i++)
			{
				if( flag_chosen_error_bit[i]==0 && usc_of_bit[i]>max_usc_of_bit)
					max_usc_of_bit=usc_of_bit[i];
			}//for i.


		}//while 2.

		if(count_of_chosen_bits<=max_trapping_set)
		if( uscc*1.0 <= count_of_chosen_bits*bc_ratio )
		{
	//		flag_find_a_trapping_set=1;
			printf("\n\n\n\nFind one trapping set: %d bits vs %d erratic checks. \n",count_of_chosen_bits, uscc);
			printf("%d Chosen bits are: ", count_of_chosen_bits);
			for(i=0;i<col_h;i++)
			{
				if(flag_chosen_error_bit[i]==1)
					printf(" %d ",i);
			}

			printf("\n%d erratic checks are: ", uscc);
			for(i=0;i<row_h;i++)
			{
				if(flag_check[i]==1)
					printf(" %d ",i);
			}


			//-- added on 20100809:
			printf("Begin to check the trapping set... \n");
			for(j=0;j<row_h;j++)
			{
				it=0;

				for(e=mod2sparse_first_in_row(H,j); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
				{
					if(flag_chosen_error_bit[e->col]==1)
						it++;
				}
				it=it%2;
				if(flag_check[j]!=it)
				{
					printf("Wrong!\n");
					exit(0);
				}//if.
			}//for j.
			printf("Correct!\n");
		}
	}//while 1.


free(flag_chosen_error_bit);
free(usc_of_bit);
free(flag_check);
}

/*

void LDPC_mod2::QC_LDPC_encoding_initial__not_full_rank()
{
//read in all generators of G matrix.
//input:
//	generator_file:	all generators in rows of non-identity-matrix part of G matrix. There are at most (t-l)+l=t generators.
//	b: the size of submatrix;
//	t: the number of circulant in H rows.
//	l: the number of submatrix in the left of H constaining the full rank of H.
//	dimen: dimension of H;
//	sub_dimen_of_Q_circulant: the number of rows in eash semi-circulant of Q part;
//

	int i,j;
	int pq; //parts of Q.
	int info_len;
	info_len = t*b-dimen;
	int parity_len;
	parity_len = l*b;

	int column_generator_len;
	column_generator_len = t*b; //only (t-l)*b + (l*b-dimen) bits are effective.

	int ** column_generator;
	column_generator = (int **)calloc(l, sizeof(int*));
	for(i=0;i<l;i++)
		column_generator[i] = (int *)calloc(column_generator_len , sizeof(int));


	int **row_generator;
	row_generator = (int **)calloc(t+l, sizeof(int*));
	for(i=0;i<t+l;i++)
		row_generator[i] = (int *)calloc(l*b, sizeof(int));


//-----------
	//-----  1: initial column_generator[]:  ---------
	// 1.1: read in row_generator:
	pq=0;
	for(i=0;i<l;i++)
	{
		if(sub_dimen_of_Q_circulant[i]>0)
			pq++;
	}

	fp=fopen(generator_file,"r");
	for(j=0;j<t-l+pq;j++)
	{
		for(i=0;i<t*b;i++) //--------------------------------------???????????????
		{
			fscanf(fp,"%d",&ct);
			if(i<l*b)
				row_generator[j][i]=ct;
		}//for i.
	}//for j.

	for(j=0;j<l;j++)
	{
		for(i=0;i<t-l+pq;i++)
		{
			for(k=0;k<l;k++)
			{
				if(k==0)
					column_generator[j][i*b+k]=row_generator[i][j*b+k];
				else //NOTE THIS!!!
					column_generator[j][i*b+k]=row_generator[i][j*b-k+b];
			}//for k.
		}//
	}//for j.



//-----------
	free(column_generator);
	for(i=0;i<l;i++)
		free(row_generator[i]);
	free(row_generator );


}
 */

void LDPC_mod2::demo_AWGN_channel()
{
//read in already-generated machine-readable G and H matrix, AWGN, and decoding.
	int i,j;

	//initial:
	initial_matrix_file();
	read_pchk(pchk_file);
	initial_decoding_space();
	read_gen(gen_file,0,0);

	//get coding rate:
	SNR_Eb_info = 4;
	SNR_Eb_code = SNR_Eb_info-10 * log10((N)*1.0/(N-M));
	awgn_sigma = ( pow(10,-SNR_Eb_code*1.0/20) ) / sqrt(2);

	//after obtaining H parameters:
	source_data   = (char *)calloc(info_length,sizeof *source_data);
	encoded_data = (char *)calloc(N,sizeof *encoded_data);
	received_block = (double *)calloc(N,sizeof *received_block);

	//encoding:
	for(i=0;i<info_length;i++)
		source_data[i]=rand()%2;

	sparse_encode(source_data,encoded_data,0);  //

	//channel:
	for(i=0;i<N;i++)
	{
		received_block[i]
			= Gauss.AWGN_noise(awgn_sigma,encoded_data[i]*2.0-1);
	}//for i.

	//decoding:
	ldpc_decode_F();

	//count errors:

//
	free(source_data);
	free(encoded_data);
}


void LDPC_mod2::demo_platform_for_performance_simulation()
{
//  read binary G into binary G format,
//	then encoding with binary G,
//	not with complex structured machine-readable LU G matrix.


// for 4x255-by-37x255, the decoding result is very strange.

	char *filen;
	char *file_matlab_rber;
	char *file_matlab_ber;
	FILE *fp;
	int i,j;
	int bl;
	int s;
	int total_eb=10;
	int total_bl=100;
	int show_unit=10;
	int max_bl=1000000;
	int be, res;

	float snr;
	long sum_bit_error;
	long sum_block_error;
	int flag_H_type;

	//initial:
	max_iter=50;
	dec_algorithm=MS;
	flag_show_error_vs_iter==1;
//	filen="./result/sim_6x127_by_72x127_8457_9144_20100827am.txt";
//	filen="./result/sim_LDPC_2x20ex1726w4_20100827am_-3_7dB.txt";
	filen="./sim_LDPC_2x24ex1438w4_20100827pm_mit=50_two_cases_4.4dB.txt";

//	file_matlab_rber="./result/matlab_raw_ber_20100827am.txt";
//	file_matlab_ber="./result/matlab_ber_20100827am_-3_7dB.txt";


	flag_H_type=0;
	flag_G_type=flag_H_type;
	flag_count_ib_error_distr=1;
	flag_count_raw_be_for_successful_decoding=1;
	flag_count_raw_be_for_failed_decoding=1;
	max_recorded_ib=100;
	max_raw_cb_count=600;

	if(flag_H_type==0)
	{//machine-readable H:

		//r=9/10;
		//H: 2x1726 by 20*1726;
//		pchk_file="./H/LDPC_2x20ex1726w4.1";
		pchk_file="./H/LDPC_2x24ex1438w4.1";
			//initial_matrix_file();
		read_pchk(pchk_file);

	}
	if(flag_H_type==2)
	{//binary H.
		//read binary matrix first:
	//	pchk_file="./H/H_4x255_by_37x255_with_pos_in_rows_not_COC_permuted.txt";//RS_LDPC_4x255_by_37x255_8434_9435.txt";
		pchk_file="./H/H_binary_6x127_by_72x127.txt";

		H=read_binary_matrix_H_to_mod2sparse(pchk_file, 72*127, 6*127);
		//H=read_H_positions_in_rows_to_mod2sparse(pchk_file, 9435,4*255, 37);

	}

	if(flag_G_type==0)
	{
		//get machine-readable G matrix:
//		gen_file="./H/LDPC_2x20ex1726w4.g";
		gen_file="./H/LDPC_2x24ex1438w4.g";
		printf("to generate the machine-readable G-LU matrix for sparse encoding.\n");
	//	make_gen_without_reading_H(pchk_file,gen_file,"sparse",0,0,0,0);
		read_gen(gen_file,0,0);
	}//
	if(flag_G_type==2)
	{
		Gb_cols=8457;
		Gb_rows=9144;
		dimension_h=687;
//		gen_file="./H/G_binary_full_8434_9435_for_4x255_by_37x255_not_COC_permuted.dat";
		gen_file="./H/G_binary_full_8457_9144_6x127_by_72x127.dat";

		FILE *fp;
//		fp=fopen(gen_file,"r");
//		fclose(fp);
		read_binary_generator_matrix_Gb(gen_file, Gb_rows, Gb_cols);
		info_length = N-dimension_h;
	}


	fp=fopen(filen,"w");
	fprintf(fp,"H=%s\n",pchk_file);
	fprintf(fp,"G=%s\n",gen_file);
	fprintf(fp,"N=%d, M=%d,info_length = %d\n",N,M,info_length);
	fprintf(fp,"max_iter=%d\n",max_iter);
	fclose(fp);

//	fp=fopen(file_matlab_rber,"w");
//	fclose(fp);

//	fp=fopen(file_matlab_ber,"w");
//	fclose(fp);

	iter_distr=(int *)calloc(max_iter,sizeof(int));
	initial_decoding_space();

	//after obtaining H parameters:
	source_data   = (char *)calloc(info_length,sizeof *source_data);
	encoded_data = (char *)calloc(N,sizeof *encoded_data);
	received_block = (double *)calloc(N,sizeof *received_block);
	if(flag_count_ib_error_distr==1)
	{
		ib_error_length_distr=(int*)calloc(max_recorded_ib,sizeof(int));
		for(i=0;i<max_recorded_ib;i++)
			ib_error_length_distr[i]=0;
	}

	if(flag_count_raw_be_for_successful_decoding==1)
	{
		raw_cb_error_distr_for_successful_decoding
			=(int*)calloc(max_raw_cb_count,sizeof(int));

		for(i=0;i<max_raw_cb_count;i++)
			raw_cb_error_distr_for_successful_decoding[i]=0;
	}
	if(flag_count_raw_be_for_failed_decoding==1)
	{
		raw_cb_error_distr_for_failed_decoding
			=(int*)calloc(max_raw_cb_count,sizeof(int));

		for(i=0;i<max_raw_cb_count;i++)
			raw_cb_error_distr_for_failed_decoding[i]=0;
	}


	for(s=0;s<=0;s++)
	{
		snr = s*0.1+4.9;
		SNR_Eb_info = snr;
		SNR_Eb_code = SNR_Eb_info-10 * log10((N)*1.0/(info_length));
		awgn_sigma = ( pow(10,-SNR_Eb_code*1.0/20) ) / sqrt(2);

		sum_block_error=0;
		sum_bit_error=0;

		for(i=0;i<max_iter;i++)
			iter_distr[i]=0;

		//for(bl=0;bl<total_bl;bl++)
		bl=0;
		raw_error=0;
		while(sum_block_error < total_eb )
		{
			bl++;
			printf("\n%1.2f SNR,%d-th block:\n",SNR_Eb_info,bl);

			//encoding:
			for(i=0;i<info_length;i++)
				source_data[i]=rand()%2;

			if(flag_G_type==0)
				sparse_encode(source_data,encoded_data,0);  //
			if(flag_G_type==2)
				encoding_with_systemic_binary_G(source_data, encoded_data);

			i = check(H, encoded_data, parity_checks); //1 for correct; 0 for fails.
			if(i==0)
			{
				printf("Encoding Wrong!\n");
				exit(0);
			}
			else
				printf("Encoding check is correct.\n");

			//channel:
			for(i=0;i<N;i++)
			{
				received_block[i]
					= Gauss.AWGN_noise(awgn_sigma,encoded_data[i]*2.0-1);
			}//for i.

			//get raw ber:
			be=0;
			for(i=0;i<N;i++)
			{
				if( (received_block[i] >0 && encoded_data[i]==0)
					||
					(received_block[i] <=0 && encoded_data[i]==1)
				   )
				    be++;
			}
			printf("raw error is %d\n",be);
			raw_error += be;


			//decoding:
			res=ldpc_decode_F();
			sum_block_error += 1 - res;

			if(0)
			{
				int r;
				scanf("%d",&r);
				exit(0);
			}

			//record info error distr:
			if(flag_count_ib_error_distr==1 )
			{
				if( info_bit_error_in_one_block >= max_recorded_ib)
					info_bit_error_in_one_block =max_recorded_ib-1;
				ib_error_length_distr[info_bit_error_in_one_block]++;
			}

			//record raw be for successful decoding and failed decoding:
			if(res==1)
			{
				if(flag_count_raw_be_for_successful_decoding==1)
				{
					if(be>=max_raw_cb_count)
						be=max_raw_cb_count-1;
					raw_cb_error_distr_for_successful_decoding[be]++;
				}
			}
			else
			{
				if(flag_count_raw_be_for_failed_decoding==1)
				{
					if(be>=max_raw_cb_count)
						be=max_raw_cb_count-1;
					raw_cb_error_distr_for_failed_decoding[be]++;
				}
			}

			//count output error:
			sum_bit_error += code_bit_error_in_one_block;

			if(sum_block_error >= total_eb || (bl)%show_unit==0)
			{
				raw_ber=raw_error*1.0/(bl+1)/N;

				fp=fopen(filen,"a");
				fprintf(fp,"\n==\nSNR_Eb_info=%f, block=%d, raw_ber=%e,\ncbits error=%d, sum_block_error=%d, cber=%e.\n",
					SNR_Eb_info, bl+1, raw_ber,sum_bit_error, sum_block_error,
					sum_bit_error*1.0/(bl+1)/N);
				fprintf(fp,"BLER=%f\n",sum_block_error*1.0/(bl+1));
				fprintf(fp,"BER of those incorrect blocks=%e.\n",
					sum_bit_error*1.0/(sum_block_error)/N);

				fprintf(fp,"Iter distr (from 1) @ this snr is:\n");
				for(i=1;i<max_iter;i++)
//					fprintf(fp,"%d>%d. ",i,iter_distr[i]);
					fprintf(fp,"%d ",iter_distr[i]);
				fprintf(fp,"\n");

				fprintf(fp,"infor error distribution (from 1) is:\n");
				for(i=1;i<max_recorded_ib;i++)
//					fprintf(fp,"%d>%d ",i,ib_error_length_distr[i]);
					fprintf(fp,"%d ",ib_error_length_distr[i]);

			if(flag_count_raw_be_for_successful_decoding==1)
			{
				fprintf(fp,"\nraw error distribution (from 1) for successful decoding is:\n");
				for(i=1;i<max_raw_cb_count;i++)
//					fprintf(fp,"%d>%d ",i,ib_error_length_distr[i]);
					fprintf(fp,"%d ",
						raw_cb_error_distr_for_successful_decoding[i]);
			}
			if(flag_count_raw_be_for_successful_decoding==1)
			{
				fprintf(fp,"\nraw error distribution (from 1) for failed decoding is:\n");
				for(i=1;i<max_raw_cb_count;i++)
//					fprintf(fp,"%d>%d ",i,ib_error_length_distr[i]);
					fprintf(fp,"%d ",
						raw_cb_error_distr_for_failed_decoding[i]);
			}

				if(sum_block_error >= total_eb || bl>=max_bl)
				{
					fprintf(fp,"\none snr end.\n");
					fclose(fp);

//					fp=fopen(file_matlab_ber,"a");
//					fprintf(fp,"%1.2f	%e	%e\n",
//						SNR_Eb_info, raw_ber, sum_bit_error*1.0/(bl+1)/N);
//					fclose(fp);

					break; //break while.
				}//if.
				fclose(fp);
			}//if.
		}//while.

	}//for s.

	free(source_data);
	free(encoded_data);
	free(iter_distr);
	if(flag_count_ib_error_distr==1)
		free(ib_error_length_distr);
	if(flag_count_ib_error_distr==1)
		free(raw_cb_error_distr_for_successful_decoding);
	if(flag_count_raw_be_for_successful_decoding==1)
		free(raw_cb_error_distr_for_failed_decoding);

}

void LDPC_mod2::demo_change_mod2sparse_to_row_position_file()
{
	char *filenam;
	filenam="./result/col_position_for_H_2x1438_by_24x1438_LDPC_2x24ex1438w4.txt";

	pchk_file="./H/LDPC_2x24ex1438w4.1";
	read_pchk(pchk_file);
	print_mod2sparse_matrix_to_col_position_file(H,filenam);
}
