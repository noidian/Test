#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "LDPC_mod2.h"

//------------------------------------------------------------------------------------
// LDPC decoding: [11/14/2009]
//------------------------------------------------------------------------------------

//

//------------------------------------------------------------------------------------
void LDPC_mod2::demo_all_process_with_machine_readable_format()
{//example of using this class:
//This function demonstrates all large functions such as
// construct H, generate G from H, and encoding with machine-readable G,
//	and AWGN channel.


	float average_girth;

	if(0)
	{//Generate H matrix first, if H is not generate.
	   char filename[200];
	   int nf = 1;   //the number of files/ldpc codes.

	   //input:
	   int ex_factor = 512;//86; //the size of a basic square matrix.
		int index, i, j;
	   WD_vector *row_dt, *col_dt;
 	   //the element of base matrix is a small circulant matrix.

	   int sub_weight=1; //weight of circulant matrix;
	   int row_base_H=4; //row weight of base matrix;
	   int column_base_H=68; //column weight of base matrix;

	   srand( (unsigned)time(NULL) ); // Ramdomize seed
	   printf("generate %d quasi-cyclic LDPC codes: \n", nf);
	   sprintf(filename, "LDPC_4x68ex512w4_20091130pm.1");

	//	gen_quasi_ldpc_files(filename, nf, 3, 51, 1, ex_factor);
	   gen_quasi_ldpc_files(filename, nf, row_base_H, column_base_H, sub_weight, ex_factor); //H: (ex_factor*4) by (ex_factor*68);

		printf("\b begin to free memory.\b");
		//return 0;
		free(col_dt->wd);
		free(row_dt->wd);
		free(row_dt);
		free(col_dt);
	}

	if(0)
	{//if generator matrix G is not generated, need to generate it from H matrix first.
		pchk_file="./H/LDPC_4x68ex512w4_20091130pm";
		gen_file="./H/gen_file_LDPC_4x68ex512w4_20091130pm";
		make_gen_main(pchk_file,gen_file,"sparse",0,0,0,0);
	}//if 0.

	//
	initial_matrix_file(); //=> this is based on that generator matrix is already generated!!!
//	pchk_file="./H/LDPC_2x40ex863w4.1";
//	gen_file="./H/LDPC_2x40ex863w4.g";
	read_pchk(pchk_file);
	initial_decoding_space();
	read_gen(gen_file,0,0);

	find_girth(H,&Min_Girth, &average_girth,0);

	//--
	source_data   = (char *)calloc(info_length,sizeof *source_data);
	encoded_data = (char *)calloc(N,sizeof *encoded_data);


	sparse_encode(source_data,encoded_data,0);

	//--
	received_block = (double *)calloc(N,sizeof *received_block);

	ldpc_decode_F();

	free(source_data);
	free(encoded_data);

}

void LDPC_mod2::demo_reading_binary_matrix_H_to_mod2sparse_and_systemic_G()
{//20100706PM: not checked.
//This function demonstrates the process of using human-readable binary H and binary sysmtemic G to encode and decode.


//	initial_matrix_file();
	sprintf(pchk_file,"matrix_H.txt");
	sprintf(gen_file,"matrix_G.txt");

	N=2048;
	M=384;
	Gb_rows=2048;
	Gb_cols=1723;

	H=read_binary_matrix_H_to_mod2sparse(pchk_file, N, M);
	read_binary_generator_matrix_Gb(gen_file, Gb_rows, Gb_cols );

	check_binary_G_for_mapping_infor_to_code();

	initial_decoding_space();

	source_data   = (char *)calloc(Gb_cols,sizeof *source_data);
	encoded_data = (char *)calloc(N,sizeof *encoded_data);
	encoding_with_systemic_binary_G(source_data, encoded_data);

	received_block = (double *)calloc(N,sizeof *received_block);

	ldpc_decode_F();

	free(source_data);
	free(encoded_data);
	free(received_block);
return;
}

void LDPC_mod2::demo_read_binary_H_and_encode_with_RU_method()
{//20100707am: working.
//this is to demo read in binary H, and encode with permuted H through RU method.
//Note: no generator matrix is needed:

	int i,j;

	N=2048;
	M=384;

	//read H:
	pchk_file="./H/RS_LDPC_384x2048/H_10G_BaseT_nonGalois.txt";
	H=read_binary_matrix_H_to_mod2sparse(pchk_file, N, M);

	//read permuted_H:
	Ht_rows=325;
	Ht_cols=2048;
	sprintf(permuted_h_file,"./H/RS_LDPC_384x2048/6.1_row_and_column_permuted_H_with_T_and_all_dependent_rows.txt");
	read_binary_permuted_H_for_RU_encoding(permuted_h_file, Ht_rows, Ht_cols);

	//read column swap tracker:
	sprintf(columnswaptracker,"./H/RS_LDPC_384x2048/1723_2048_matrix_with_t_equals_325_ColumnSwapTracker_.txt");
	read_column_swap_tracker(columnswaptracker,2048);


	//encoding:
	source_data   = (char *)calloc(Ht_rows,sizeof *source_data);
	encoded_data = (char *)calloc(N,sizeof *encoded_data);

	size_T=325;
	size_g=0;
	encoding_in_RU_method_through_permuted_H_with_g_as_0(source_data,encoded_data);

//	check_encoded();

	free(source_data);
	free(encoded_data);
return;
}
//
unsigned LDPC_mod2::ldpc_decode_F()
{
//input: received_block
//input length:
//output:
//output length:

	int dec_over_interation;

    if (dec_algorithm == MS)
	{
		/* Min Sum algorithm */
		//fprintf(stderr, "\n Min Sum decoding begin... \n");
        dec_over_interation = ldpc_decode_MS();
    }
	else
	{
		/* BF algorithm (hard decoding, no soft info) */
		//fprintf(stderr, "\n BF hard decoding begin... \n");
        dec_over_interation = ldpc_decode_BF();
	}

	// the interation index when decoding finishes
	return dec_over_interation;
}

unsigned LDPC_mod2::ldpc_decode_MS()
{
//input: received_block;
//
	int j;
	int n,whole_valid;
	int it;
	int i;

	/* Initialize llr: probability and likelihood ratios, and find initial guess. */
	initprp_log_MS(); //=> output: H
	iteration = 0;
	whole_valid = check(H, decoded_codeword, parity_checks);


	/* Do up to abs(max_iter) iterations of probability propagation, stopping
	   early if a codeword is found, unless max_iter is negative. */
    info_bit_error_in_one_block=0;
	for (n = 0; n!=max_iter && n!=-max_iter && (max_iter<0 || !whole_valid); n++)
	{
		iteration = n+1;

		//fprintf(stderr, "\n iteration = %d", iteration);

		//input: H, not llr_fn!!!
		iterprp_log_MS(H,llr_fn,decoded_codeword,(char*)0,(char*)0);

		whole_valid = check(H, decoded_codeword, parity_checks);

/**/
		if(flag_show_error_vs_iter==1)
		{
			it=0;
			for(i=0;i<N;i++)
			{
				if(decoded_codeword[i]!=encoded_data[i])
					it++;
			}
			printf("%d iteration: %d code bits error.\n",n,it);
		}

		if(flag_count_ib_error_distr==1)
		{
			info_bit_error_in_one_block=0;
			for (j = M; j<N; j++)
			{
				if(decoded_codeword[cols[j]] != source_data[j-M])
					info_bit_error_in_one_block++;
			}
			printf("info_bit_error_in_one_block=%d\n",info_bit_error_in_one_block);
		}
			if (whole_valid)
	    {
			printf("\nIteration = %d  decoding succeed! \n", iteration);
	    }
	}

	//count bit errors:
	it=0;
	for(i=0;i<N;i++)
	{
		if(decoded_codeword[i]!=encoded_data[i])
			it++;
	}
	code_bit_error_in_one_block=it;

	//
	if(whole_valid)
		iter_distr[n]++;

	//if (!whole_valid)
	//	tot_err_blks ++;

return whole_valid;
}

/* INITIALIZE PROBABILITY PROPAGATION --- LOG.  Stores initial ratios, probabilities,
   and guess at decoding. */

void LDPC_mod2::initprp_log_MS()
{
//input:  received_block;
//input:  H: it includes n, m, llr_fn memory;

//output: H, decoded_codeword;

	mod2entry *e;
	int n;
	int j;

	n = mod2sparse_cols(H);

	for (j = 0; j<n; j++)
	{
        //Sat and Quan for LDPC decoder
//      llr_fn[j] = Sat_Quan(received_block[j],Finite_Lc_Plus_Max,Finite_Lc_Minus_Max,Finite_Lc_Tot,Finite_Lc_Frac);
		llr_fn[j] = received_block[j];///awgn_sigma;
		for (e = mod2sparse_first_in_col(H,j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->pr = llr_fn[j];
		}

		// hard decision
		decoded_codeword[j] = (received_block[j] >= 0);
	}
}

/* SEE IF A GUESS CHECKS OUT.  Returns 1 if dblk contains a valid codeword.
   Whether or not it does, the results of all the parity checks are stored
   in pchk. */

int LDPC_mod2::check ( mod2sparse *H,	/* Parity check matrix */
            char *dblk,		/* Guess for codeword */
            char *pchk		/* Place to store parity checks */
          )
{
	int M;
	int i;

	M = mod2sparse_rows(H);

	mod2sparse_mulvec (H, dblk, pchk);

	for (i = 0; i<M; i++)
	{
		if (pchk[i]!=0)
            return 0;
	}

	return 1;
}


/* DO ONE ITERATION OF PROBABILITY PROPAGATION --- LOG. */

void LDPC_mod2::iterprp_log_MS( mod2sparse *H,      /* Parity check matrix */
                     float *lratio_log, /* initial Likelihood ratios log for bits */
                     char *dblk,         /* Place to store decoding */
                     char *look_row,     /* Which rows to look at, look at all if null */
                     char *look_col      /* Which columns to look at, look at all if null */
                   )
{
//input:
//output: code_bit_error_in_one_block ;


	double min1;
	double min2;
	int min1_id;
	double sum;
	int sign;
	double abs_e_pr;
	mod2entry *e;
	int N, M;
	int i, j;
	int bitErrs;

	M = mod2sparse_rows(H);
	N = mod2sparse_cols(H);

	/* Update likelihood ratios. */
	/* CNU processing */
		for (i = 0; i<M; i++)
		{
			if (look_row!=0 && !look_row[i]) continue;

			min1 = 100000;
			min2 = 100000;
			min1_id = 0;
			sign = 1;

			for (e = mod2sparse_first_in_row(H,i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
			{
				if ( e->pr > 0 )
					abs_e_pr = e->pr;
				else
					abs_e_pr = -1 * e->pr;

				if ( abs_e_pr < min1 )
				{
					min2 = min1;
					min1 = abs_e_pr;
					min1_id = e->col;
				}
				else if ( abs_e_pr < min2 )
				{
					min2 = abs_e_pr;
				}

				sign *= ( e->pr > 0 ? 1 : -1 );
			}

			for (e = mod2sparse_first_in_row(H,i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
			{
				if ( min1_id == e->col )
					e->pr = sign*( e->pr > 0 ? 1 : -1 )* min2;
				else
					e->pr = sign*( e->pr > 0 ? 1 : -1 )* min1;
			}
		}

	/* Update probability ratios.  Also find the next guess based on the
	   individually most likely values. */
	/* VNU processing */
	bitErrs = 0;
	for (j = 0; j<N; j++)
	{ if (look_col!=0 && !look_col[j]) continue;

		sum = lratio_log[j];

		for (e = mod2sparse_first_in_col(H,j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			sum += e->pr;
			//sum = Sat_Quan(sum,Finite_Lc_Plus_Max,Finite_Lc_Minus_Max,Finite_Lc_Tot,Finite_Lc_Frac);
		 }

		dblk[j] = sum > 0; //=> 1:llr>0;  <-> 0:llr<0;

//----- count code bit errors:
		if(dblk[j]!=encoded_data[j]) {
			bitErrs++;
		//	IterBitErr[iteration][j] = 1;
		}
		else {
		//	IterBitErr[iteration][j] = 0;
		}
/**/

//----- for sova feedback info:
		 //softOut[j] = Sat_Quan((sum - rblk_f[j]),Finite_Info_Lc_Plus_Max,Finite_Info_Lc_Minus_Max,Finite_Info_Lc_Tot,Finite_Info_Lc_Frac);
	//		softOut[j] = (sum - llr_fn[j]);	//=> used when need to feedback soft info to sova.

		//IterSoft[iteration][j] = softOut[j];

//-----
		for (e = mod2sparse_first_in_col(H,j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
#ifdef sim_decoding_quantization
			e->pr = Sat_Quan((sum - e->pr),Finite_Lc_Plus_Max,Finite_Lc_Minus_Max,Finite_Lc_Tot,Finite_Lc_Frac);
#else
			e->pr = sum - e->pr;
#endif
		}
	}

#ifdef show_code_bit_error
	printf( "%d code error\n", bitErrs);
#endif

//	code_bit_error_in_one_block = bitErrs;

return;
}

//
/* saturate and quantize x to m total bits with k fractional bits */
double LDPC_mod2::Sat_Quan(double x, double Max_V, double Min_V, int m, int k)
{
  double val;
  int q,p;

  val = x;
  val = val > Max_V ? Max_V : val;
  val = val < Min_V ? Min_V : val;

  q = Quantize(val, m, k);
  Tru2intS(q,m,&p);
  return (double) p / pow(2, k);
}


int LDPC_mod2::Quantize(double x, int m, int k)
{
  int t, t1;
  t=1<<k;
  t=(int)(x*t+(x>0 ? 0.5 : -0.5));
  t1=(1<<m)-1;
  t1=t1 & t;
  return t1;
}

/*   This func to convert m-bit integer */
/*   a signed integer  */
void LDPC_mod2::Tru2intS(int x, int m, int *dx)
{
  int t1=x;
  int t2=-1;
  t1=(x>>(m-1));
  if (t1==0) (*dx)=x;
  else {
        t1=(1<<m)-1;
        t2=t2 ^ t1;  // ^ is 'xor' operation;
        *dx= t2 | x;
  }
}

unsigned LDPC_mod2::ldpc_decode_BF()
{
    char *pcfail_num; //the number of parity checks for each code bit;
	int max_pcfail;
	int  n,whole_valid;
	int  N;
	int  j;
	mod2entry *e;
	int  fail_vn;

	max_pcfail = 0;

    M = mod2sparse_rows(H);
	N = mod2sparse_cols(H);

	pcfail_num = (char *)calloc(N, sizeof(char));
	iteration = 0;
    /* Initialize initial decision infomation from likelihood. */
	for (j = 0; j<N; j++)
	{
		// hard decision
		decoded_codeword[j] = received_block[j] < 0; //: 信道是用1->-1; 0-> 1;

        // initial the parity check fail counter for each variable node
		pcfail_num[j] = 0;
	}

    whole_valid = check(H,decoded_codeword,parity_checks);

    //if (whole_valid)
	//{
	 //   fprintf(stderr, "\n iteration = 0  decoding succeed! \n");
	//}

	/* Do up to abs(max_iter) iterations of probability propagation, stopping
	   early if a codeword is found, unless max_iter is negative. */

	for (n = 0; n!=max_iter && n!=-max_iter && (max_iter<0 || !whole_valid); n++)
	{
		iteration = n+1;

        if (!whole_valid)
		{
		    // check in which row parity check failed
			// If one row failed, the parity check fail counter (pcfail_num) of every variable node in this row add 1
			for (j=0; j<M; j++)
			{
                if (parity_checks[j] != 0)
				{
                    for (e = mod2sparse_first_in_row(H,j); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
					{
						fail_vn = e->col;
						pcfail_num[fail_vn]++;
					}
				}
			}

			// choose the maximum pcfail_num
			for (j=0; j<N; j++)
			{
                if (pcfail_num[j] > max_pcfail)
					max_pcfail = pcfail_num[j];
			}

			// flip the bit with the maximum pcfail_num
			for (j=0; j<N; j++)
			{
                if (pcfail_num[j] == max_pcfail)
					decoded_codeword[j] = (decoded_codeword[j]==0);
			}
		}

		whole_valid = check(H,decoded_codeword,parity_checks);

        //if (whole_valid)
	    //{
		//	fprintf(stderr, "\n iteration = %d  decoding succeed! \n", iteration);
	    //}
	}

	//if (!whole_valid)
	//	tot_err_blks ++;
	free(pcfail_num);
	return whole_valid;

}


//----------------------------------------------------------------------------------
// LDPC encoding:
//----------------------------------------------------------------------------------
void LDPC_mod2::sparse_encode (char * source, char * encoded,int array_start)
{//use LU decomposition to encode:
//Make sure that gen_file is already generated with
//	make_gen_main() function.
// This gen_file contains special information of LU matrix used for encoding.
// If you just have binary G matrix, you can not
// directly use this encoding function--you need to generate
// the special machine-readable G matrix first (including L and U).

//Note!!: if H is not full-rank, the count of information bits is larger than N-M.
//			In this case, we only use N-M infor bits, and the left is set to be 0.

  int i, j;

  mod2entry *e;
  char *x, *y;

  FILE *fp;

  x = (char *) calloc(M, sizeof (x));
  y = (char *) calloc(M, sizeof (y));

  /* Multiply the vector of source bits by the systematic columns of the
     parity check matrix, giving x.  Also copy these bits to the coded block. */

  for (i = 0; i<M; i++) x[i] = 0;

//
#ifdef show_info1
  fp=fopen("mapping1.txt","w");
  for (j = M; j<N; j++)
  {
    fprintf(fp,"%d\n",cols[j]);
  }
  fclose(fp);
#endif
//


  for (j = M; j<N; j++)
  {
    encoded[cols[j]] = source[array_start+j-M];

    if (source[array_start+j-M]==1)
    {
		for(e = mod2sparse_first_in_col(H,cols[j]);
           !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e))
      {
		   x[mod2sparse_row(e)] ^= 1;
      }
    }
  }

  /* Solve Ly=x for y by forward substitution, then U(cblk)=y by backward
     substitution. */

  if (!mod2sparse_forward_sub(L,rows,x,y)
   || !mod2sparse_backward_sub(U,cols,y,encoded))
  {
    abort(); /* Shouldn't occur, even if the parity check matrix has
                redundant rows */
  }

  free(x);
  free(y);

}

//--------------------------------------------------------------
//  alloc.c
//--------------------------------------------------------------
/* ALLOCATE SPACE AND CHECK FOR ERROR.  Calls 'calloc' to allocate space,
   and then displays an error message and exits if the space couldn't be
   found. */

void *  LDPC_mod2::chk_alloc
( unsigned n,		/* Number of elements */
  unsigned size		/* Size of each element */
)
{
  void *p;

  p = calloc(n,size);

  if (p==0)
  { fprintf(stderr,"Ran out of memory (while trying to allocate %d bytes)\n",
      n*size);
    exit(1);
  }

  return p;
}

//
void LDPC_mod2::init_quantization()
{//initialize some values for quantization functions:
//	Finite_Info_Plus_Max = pow(2,Finite_Info_Tot-Finite_Info_Frac-1)-pow(2,-Finite_Info_Frac);
//	Finite_Info_Minus_Max = -pow(2,Finite_Info_Tot-Finite_Info_Frac-1);
	Finite_Lc_Plus_Max = pow(2,Finite_Lc_Tot-Finite_Lc_Frac-1)-pow(2,-Finite_Lc_Frac);
	Finite_Lc_Minus_Max = -pow(2,Finite_Lc_Tot-Finite_Lc_Frac-1);
//	Finite_Info_Lc_Plus_Max = pow(2,Finite_Info_Lc_Tot-Finite_Info_Lc_Frac-1)-pow(2,-Finite_Info_Lc_Frac);
//	Finite_Info_Lc_Minus_Max = -pow(2,Finite_Info_Lc_Tot-Finite_Info_Lc_Frac-1);
return;
}


void LDPC_mod2::initial_matrix_file()
{
	time_t rawtime;
	struct tm * ptm;
	time ( &rawtime );
	ptm = gmtime ( &rawtime );
	printf( "\nOS time: %d:%d:%d, %d/%d/%d\n\n", ptm->tm_hour, ptm->tm_min, ptm->tm_sec, 1+ptm->tm_mon, ptm->tm_mday, 1900+ptm->tm_year);

//	pchk_file="E:\\matrix\\matrix";
//	gen_file="E:\\matrix\\gen_file";

	max_iter=100;
    flag_show_error_vs_iter==1;
    flag_count_ib_error_distr=1;

//	pchk_file="./H/matrix";
//	gen_file="./H/gen_file";

//---------------------------------------------------------
/*
	//H: 512 by 4608;
	//G: 4608 by 4096;
	pchk_file="./H/matrix_H_512_4608";
	gen_file="./H/gen_file_H_512_4608";
*/

//---------------------------------------------------------

	//r=9/10;
	//H: 2x1726 by 20*1726;
	//G: ;
	pchk_file="./H/LDPC_2x20ex1726w4.1";
	gen_file="./H/LDPC_2x20ex1726w4.g";





//---------------------------------------------------------
/*
	//r=16/17;
	//H: 2048 by 34816;
	//G: 34816 by 32768;
	pchk_file="./H/H_4x68ex512w4_20091130pm";
	gen_file="./H/gen_file_LDPC_4x68ex512w4_20091130pm";
*/

//---------------------------------------------------------
// the LDPC used for flash-equ-ldpc paper [4/7/2010]

	//r=19/20;
	//H: 2*863=1726 by 40*863=34520;
	//G:
//	pchk_file="./H/H_LDPC_2x40ex863w4.1";
//	gen_file="./H/G_LDPC_2x40ex863w4.g";
//---------------------------------------------------------

//---------------------------------------------------------
// the LDPC used for Xie Ningde's paper [6/7/2010]
	//r=19/20;
	//H: 2*863=1726 by 40*863=34520;
	//G:
	 //!!!!!xnd: set pchk and g files:
//	pchk_file="./H/LDPC_2x40ex863w4.1";
//	gen_file="./H/LDPC_2x40ex863w4.g";
//---------------------------------------------------------


	//
	printf("The LDPC matrixes are: \n	%s | %s.\n",pchk_file,gen_file);

	dec_algorithm=MS;
//	wordlength=10;
//	max_iter=16;
//	stop_method=0;
//	Sim_Length=1000000;

#ifdef show_info1
	printf("\n wordlength=%d; max_iter=%d;\n",wordlength,max_iter);
#endif

//----------------

}

//

/* READ PARITY CHECK MATRIX.  Sets the H, M, and N global variables.  If an
   error is encountered, a message is displayed on standard error, and the
   program is terminated. */

void LDPC_mod2::read_pchk
( char *pchk_file
)
{
  FILE *f;
  FILE *f2;

  printf("\n Read machine readable matrix information.\n");

  f = fopen(pchk_file,"rb");
//  f = fopen("./H/matrix","rb");
  if (f==NULL)
  { fprintf(stderr,"Can't open parity check file: %s\n",pchk_file);
    exit(1);
  }

  if (intio_read(f)!=('P'<<8)+0x80)
  { fprintf(stderr,"File %s doesn't contain a parity check matrix\n",pchk_file);
    exit(1);
  }

  H = mod2sparse_read(f);

#ifdef show_info1
  //output file:
  f2=fopen("H.txt","w");
  mod2sparse_print(f2,H);
  fclose(f2);
//  exit(0);
  //
#endif

  if (H==0)
  { fprintf(stderr,"Error reading parity check matrix from %s\n",pchk_file);
    exit(1);
  }

  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  //==
  	info_rate = (float) (N-M)/N;
	/* Initialize the COM_PK defined in ldpc.h */
	info_length = N - M;

#ifdef show_info1
  printf("rows num: %d;  columns num: %d \n",M,N);
#endif

  fclose(f);

  if (N<=M)
    {
      fprintf(stderr, "Error: N = %d but M = %d\n",N,M);
      exit(1);
    }
}

/* READ GENERATOR MATRIX.  The parity check matrix must have already been
   read, unless the last argument is set to 1.  The generator matrix must be
   compatible with the parity check matrix, if it has been read.  If the
   second argument is 1, only the column ordering (the last N-M of which are
   the indexes of the message bits) is read, into the 'cols' global variable.
   Otherwise, everything is read, into the global variables appropriate
   to the representation.  The 'type' global variable is set to a letter
   indicating which represention is used.

   If an error is encountered, a message is displayed on standard error,
   and the program is terminated. */

void LDPC_mod2::read_gen
( char *gen_file,	/* Name of generator matrix file */
  int cols_only,	/* Read only column ordering? */
  int no_pchk_file	/* No parity check file used? */
)
{
  int M2, N2;
  FILE *f;
  int i;

  f = fopen(gen_file,"rb");
  if (f==NULL)
  { fprintf(stderr,"Can't open generator matrix file: %s\n",gen_file);
    exit(1);
  }

  if (intio_read(f)!=('G'<<8)+0x80)
  { fprintf(stderr,"File %s doesn't contain a generator matrix\n",gen_file);
    exit(1);
  }

  if (fread (&type, 1, 1, f) != 1) goto error;

  M2 = intio_read(f);
  N2 = intio_read(f);

  if (feof(f) || ferror(f)) goto error;

  if (no_pchk_file)
  { M = M2;
    N = N2;
  }
  else
  { if (M2!=M || N2!=N)
    { fprintf(stderr,
              "Generator matrix and parity-check matrix are incompatible\n");
      exit(1);
    }
  }

  cols = (int *)calloc (N, sizeof *cols);
  mapping_code_bit_to_infor = cols;

  for (i = 0; i<N; i++)
  { cols[i] = intio_read(f);
    if (feof(f) || ferror(f)) goto error;
  }

  if (!cols_only)
  {
    switch (type)
    {
      case 's':
      {
        rows = (int *) calloc (M, sizeof *rows);

        for (i = 0; i<M; i++)
        { rows[i] = intio_read(f);
          if (feof(f) || ferror(f)) goto error;
        }

        if ((L = mod2sparse_read(f)) == 0) goto error;
        if ((U = mod2sparse_read(f)) == 0) goto error;

        if (mod2sparse_rows(L)!=M || mod2sparse_cols(L)!=M) goto garbled;
        if (mod2sparse_rows(U)!=M || mod2sparse_cols(U)<M) goto garbled;

        break;
      }

      case 'd':
      {
        if ((G = mod2dense_read(f)) == 0) goto error;

        if (mod2dense_rows(G)!=M || mod2dense_cols(G)!=N-M) goto garbled;

        break;
      }

      case 'm':
      {
        if ((G = mod2dense_read(f)) == 0) goto error;

        if (mod2dense_rows(G)!=M || mod2dense_cols(G)!=M) goto garbled;

        break;
      }

      default:
      { fprintf(stderr,
         "Unknown type of generator matrix in file %s\n",gen_file);
        exit(1);
      }
    }
  }

  fclose(f);

//  going on--
//	  此处出错,因为此时没有生成G,而是用的另外一个case,得到的是L,U
 // f=fopen("G.txt","w");
 // mod2dense_print(f,G);


  return;

error:
  fprintf(stderr,"Error reading generator matrix from file %s\n",gen_file);
  exit(1);

garbled:
  fprintf(stderr,"Garbled generator matrix in file %s\n",gen_file);
  exit(1);
}

//--------------------------------------------------------------------
//--graph.c
//--------------------------------------------------------------------

void LDPC_mod2::init_node_vector( node_vector *v,int n)
{
  v->n = n;
  v->entry = (node_entry *) chk_alloc(n, sizeof(node_entry));
}

/* check the number of connected component in the current graph */
int LDPC_mod2::Comp_Num
(
 mod2sparse *H,            /* Parity check matrix */
 node_vector *subcode,     /* The subcode node array */
 node_vector *digit       /* The digit node array */
)
{
  int M,N,i,j,k;
  node_entry *queue_head, *queue_tail, *current_node;
  mod2entry *e;
  int tmp_index,parent_index,current_level,count,components;
  int root_index;


  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  /* check the dimension */
  if ((subcode->n != M) || (digit->n != N))
    {
      fprintf(stderr, "ERROR: Diemnsion mismatch!\n");
      abort();
    }

  /* Initialize each node */
  for (i = 0; i < M; i++)
    {
      subcode->entry[i].color = white;
      /*  subcode->entry[i].parent = NULL; */
      subcode->entry[i].depth = -1;
      /*  subcode->entry[i].next = NULL; */
      subcode->entry[i].index = i;
    }
  for (i = 0; i < N; i++)
    {
      digit->entry[i].color = white;
      /*  digit->entry[i].parent = NULL; */
      digit->entry[i].depth = -1;
      /*  digit->entry[i].next = NULL; */
      digit->entry[i].index = i;
    }

  /* Begin from the root */
  root_index = 0;
  subcode->entry[root_index].color = gray;
  subcode->entry[root_index].depth = 0;
  queue_head = &subcode->entry[root_index];
  queue_tail = &subcode->entry[root_index];

  k = 0; /* k = 0: subcode node level; k = 1: digit node level */

  count = 1;

  components = 1;

  /* perform the Bread-first search untill all nodes are not white
     in color */
  while (count < (M+N))
    {
      while (queue_head != NULL)
	{
	  switch (k)
	    {
	    case 0:
	      /* if current node is in subcode node level */
	      e = mod2sparse_first_in_row(H, queue_head->index);
	      break;
	    case 1:
	      /* if current node is in digit node level */
	      e = mod2sparse_first_in_col(H, queue_head->index);
	      break;
	    default:
	      abort();
	    }

	  while (!mod2sparse_at_end(e))
	    {
	      switch (k)
		{
		case 0:
		  /* if current node is in subcode node level */
		  /* get the digit node index of e */
		  tmp_index = mod2sparse_col(e);
		  current_node = &digit->entry[tmp_index];
		  break;
		case 1:
		  /* if current node is in digit node level */
		  /* get the subcode node index of e */
		  tmp_index = mod2sparse_row(e);
		  current_node = &subcode->entry[tmp_index];
		  break;
		default:
		  abort();
		}

	      if (current_node->color == white)
		{
		  /* if current node color is white, then:
		     1. change color to gray;
		     2. set depth & parent;
		     3. insert digit node to queue;
		     4. increase count by 1. */

		  current_node->color = gray;
		  current_node->parent = queue_head;
		  current_node->depth = queue_head->depth + 1;
		  if (queue_head->depth == 0)
		    current_node->group_ID = current_node->index;
		  else
		    current_node->group_ID = queue_head->group_ID;
		  queue_tail->next = current_node;
		  queue_tail = current_node;
		  count++;
		}
	      switch (k)
		{
		case 0:
		  /* if current node is in subcode node level */
		  /* get the next digit node */
		  e = mod2sparse_next_in_row(e);
		  break;
		case 1:
		  /* if current node is in digit node level */
		  /* get the next subcode node */
		  e = mod2sparse_next_in_col(e);
		  break;
		}
	    }
	  queue_head->color = black;
	  queue_head = queue_head->next;

	  /* check whether next node is in subcode node level or digit node level */
	  if (queue_head != NULL)
	    k = (queue_head->depth % 2);
	}

      /* if there are more than 1 component */
      if (count < (M+N))
	{
	  components ++;

	  /* find the next root_index */
	  for (i = 0; i < M; i++)
	    if (subcode->entry[i].color == white)
	      break;

	  if (i < M) /* successfully find the next root_index */
	    {
	      root_index = i;
	      subcode->entry[root_index].color = gray;
	      subcode->entry[root_index].depth = 0;
	      queue_head = &subcode->entry[root_index];
	      queue_tail = &subcode->entry[root_index];
	      count++;
	      k = 0;
	    }
	  else
	    {
	      fprintf(stderr,"Error with the Comp_Num function!\n");
	      abort();
	    }
	}
    }
  return components;
}


/* Perform the Bread-first search operation based on the LDPC
   parity check matrix to find the first cycle, the root is
   always a subcode node.
   Return value: value of girth passing the node with coordinate of index
*/

int LDPC_mod2::BFS
(
 mod2sparse *H,            /* Parity check matrix */
 node_vector *subcode,     /* The subcode node array */
 node_vector *digit,       /* The digit node array */
 int root_index            /* the index of root */
 )
{
  int i,j,M,N,find_cycle,k;
  node_entry *queue_head, *queue_tail, *current_node;
  mod2entry *e;
  int tmp_index,parent_index,current_level;
  int trace[max_girth];
  int record_girth = 0;
  int found=0, girth=0;
  int L=64,p=6,tmp;
  FILE *fp;

  if (record_girth == 1)
    {
      fp = fopen("girth_record","a");
      fprintf(fp,"\n");
    }

  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  /* check the dimension */
  if ((subcode->n != M) || (digit->n != N))
    {
      fprintf(stderr, "ERROR: Diemnsion mismatch!\n");
      abort();
    }

  /* Initialize each node */
  for (i = 0; i < M; i++)
    {
      subcode->entry[i].color = white;
      subcode->entry[i].parent = NULL;
      subcode->entry[i].depth = -1;
      subcode->entry[i].next = NULL;
      subcode->entry[i].index = i;
    }
  for (i = 0; i < N; i++)
    {
      digit->entry[i].color = white;
      digit->entry[i].parent = NULL;
      digit->entry[i].depth = -1;
      digit->entry[i].next = NULL;
      digit->entry[i].index = i;
    }

  /* Begin from the root */
  subcode->entry[root_index].color = gray;
  subcode->entry[root_index].depth = 0;
  queue_head = &subcode->entry[root_index];
  queue_tail = &subcode->entry[root_index];

  k = 0; /* k = 0: subcode node level; k = 1: digit node level */

  /* perform the Bread-first search untill a gray node will color
     another gray node or all the graph has been searched. */
  while (queue_head != NULL)
    {

      switch (k)
	{
	case 0:
	  /* if current node is in subcode node level */
	  e = mod2sparse_first_in_row(H, queue_head->index);
	  break;
	case 1:
	  /* if current node is in digit node level */
	  e = mod2sparse_first_in_col(H, queue_head->index);
	  break;
	default:
	  abort();
	}

      while (!mod2sparse_at_end(e))
	{
	  switch (k)
	    {
	    case 0:
	      /* if current node is in subcode node level */
	      /* get the digit node index of e */
	      tmp_index = mod2sparse_col(e);
	      current_node = &digit->entry[tmp_index];
	      break;
	    case 1:
	      /* if current node is in digit node level */
	      /* get the subcode node index of e */
	      tmp_index = mod2sparse_row(e);
	      current_node = &subcode->entry[tmp_index];
	      break;
	    default:
	      abort();
	    }

	  if (current_node->color == white)
	    {
	      /* if current node color is white, then:
		 1. change color to gray;
		 2. set depth & parent
		 3. insert digit node to queue */

	      current_node->color = gray;
	      current_node->parent = queue_head;
	      current_node->depth = queue_head->depth + 1;
	      if (queue_head->depth == 0)
		current_node->group_ID = current_node->index;
	      else
		current_node->group_ID = queue_head->group_ID;

	      queue_tail->next = current_node;
	      queue_tail = current_node;
	    }
	  else if (current_node->color == gray)
	    /* if current node is already gray, then a circle has been found */
	    /* if current node groupd ID is not equal to head node ID, then
	       girth pass through root */
	    if (current_node->group_ID != queue_head->group_ID)
	      {
		/* trace back find the cycle */
		if (record_girth == 1)
		  {
		    if ((current_node->depth % 2) == 0)
		      {
			fprintf(fp,"subcode ");
			current_level = 0;
		      }
		    else
		      {
			fprintf(fp,"digit ");
			current_level = 1;
		      }

		    /*  fprintf(fp,"%d ", current_node->index); */
		    if (current_level == 0)
		      {
			tmp = current_node->index/L;
			tmp = (tmp >= p)? tmp-p : tmp;
			fprintf(fp,"(%d %d)", tmp, current_node->index%L);
		      }

		    while ((current_node->index != root_index) || (current_level != 0))
		      {
			current_node = current_node->parent;
			/* fprintf(fp,"%d", current_node->index); */
			if (current_level == 1)
			  {
			    tmp = current_node->index/L;
			    tmp = (tmp >= p)? tmp-p : tmp;
			    fprintf(fp,"(%d %d)", tmp, current_node->index%L);
			  }
			current_level = (current_level + 1) % 2;
		      }

		    current_node = queue_head;

		    for (i = 0; i < queue_head->depth; i++)
		      {
			trace[i] = current_node->index;
			current_node = current_node->parent;
		      }

		    current_level = 1;

		    for (i = queue_head->depth-1; i > -1; i--)
		      {
			if (current_level == 0)
			  {
			    tmp = trace[i]/L;
			    tmp = tmp >= p? tmp-p : tmp;
			    fprintf(fp,"(%d %d)", tmp, trace[i]%L);
			  }
			current_level = (current_level + 1) % 2;
			/*  fprintf(fp,"%d ",trace[i]); */
		      }
		    fprintf(fp,"\n");
		    /* fclose(fp); */
		  }

		/* recrod all the cycles instead of only first one */
		if (found == 0)
		  {
		    found = 1;
		    girth = 2*(queue_head->depth+1);
		  }
		/* return 2*(queue_head->depth+1); */
	      }

	  switch (k)
	    {
	    case 0:
	      /* if current node is in subcode node level */
	      /* get the next digit node */
	      e = mod2sparse_next_in_row(e);
	      break;
	    case 1:
	      /* if current node is in digit node level */
	      /* get the next subcode node */
	      e = mod2sparse_next_in_col(e);
	      break;
	    }
	}

      queue_head->color = black;
      queue_head = queue_head->next;
      /* recrod all the cycles instead of only first one */
      if (found == 1)
	if (girth < 2*(queue_head->depth+1))
	  {
	    if (record_girth == 1)
	      {
		fprintf(fp,"\n");
		fclose(fp);
	      }
	    return girth;
	  }

      /* check whether next node is in subcode node level or digit node level */
      if (queue_head != NULL)
	k = (queue_head->depth % 2);
    }
  return -1;
}


/* find the girth distribution of a bipartite graph */
void LDPC_mod2::find_girth
(
 mod2sparse *H,            /* Parity check matrix */
 int *min_girth,
 float *avg_girth,
 FILE *fp
 )
{
  int i;//,M2,N2;
  node_vector *subcode;     /* The subcode node array */
  node_vector *digit;       /* The digit node array */
  int girth[max_girth],tmp,tot_girth,tot_value,comp;
  int record;

  if (fp != 0)
    record = 1;
  else
    record = 0;

  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  /* allocate the subcode node vector and digit node vector */
  /*    init_node_vector(&subcode,M);  */
  /*    init_node_vector(&digit,N); */
  subcode = (node_vector *)calloc(1, sizeof *subcode);
  subcode->n = M;
  subcode->entry = (node_entry *)calloc(M, sizeof(node_entry));

  digit = (node_vector *)calloc(1, sizeof *digit);
  digit->n = N;
  digit->entry = (node_entry*)calloc(N, sizeof(node_entry));

  /* check the number of connected components */
  comp = Comp_Num(H, subcode, digit);

  fprintf(stderr,"LDPC bipartite graph contains %d connected components!\n",comp);

  for (i = 0; i < max_girth; i++)
    girth[i] = 0;

  /* do the Breadth-first search from each subcode */
  for (i = 0; i < M; i++)
    {
      tmp = BFS(H, subcode, digit, i);

      if (tmp > 0)
	{
	  tmp = tmp/2 - 2;
	  if (tmp == -1)
	    {
	      fprintf(stderr,"Error: girth is equal to 2!\n");
	      abort();
	    }
	  girth[tmp]++;
	}
    }

  /* find the min_girth */
  for (i = 0; i < max_girth; i++)
    if (girth[i] > 0)
      {
	min_girth[0] = (i+2)*2;
	break;
      }

  /* calculate the average girth */
  tot_girth = 0;
  tot_value = 0;
  for (i = 0; i < max_girth; i++)
    {
      if (girth[i] != 0)
	{
	  fprintf(stderr,"%4d points with girth of %d\n", girth[i],(i+2)*2);
	  if (record == 1)
	    fprintf(fp,"%4d points with girth of %d\n", girth[i],(i+2)*2);
	}

      tot_girth +=girth[i];
      tot_value += girth[i]*(i+2)*2;
    }
  avg_girth[0] = (float)  tot_value/tot_girth;
  if (record == 1)
    {
      fprintf(fp,"LDPC matrix min_girth: %d; Average girth:* %f.\n",min_girth[0],avg_girth[0]);
    }

  free(subcode->entry);
  free(digit->entry);

}

//-------------------------
LDPC_mod2::LDPC_mod2()
{
//initial value:
	M=0;
	N=0;
	max_iter=-1;
	flag_G_binary_allocated=0;
	flag_mapping_infor_to_code_allocated=0;
	flag_permuted_H_allocated=0;
	flag_columnswaptracke_allocated=0;

#ifdef sim_decoding_quantization
	init_quantization();
#endif

	Gauss.Gauss_initial();
}


//------------------------------
void LDPC_mod2::initial_decoding_space()
{
	if(M==0 || N==0)
	{
		printf("Error: M=N=0! Please run read_pchk() first to get N and M.\n");
		exit(1);
	}
	decoded_codeword = (char *)calloc(N,sizeof *decoded_codeword);
	llr_fn = (float *)calloc(N, sizeof *llr_fn);
	parity_checks = (char *)calloc(M,sizeof*parity_checks);
}


//------------------------------
unsigned LDPC_mod2::ldpc_decode_MS_without_llr_init()
{//this function works as MinSum decoding, without llr initial for AWGN channel.
//Suppose the input data to llr_fn is already given as llr.

	//input: llr_fn; encoded_data;
	//output: decoded_codeword;

	int n,whole_valid;

//	/* Initialize llr: probability and likelihood ratios, and find initial guess. */
//	initprp_log_MS();

//	show_array_float("llr_2.txt",llr_fn,N);

	//init decoded_codeword:
	for(n=0;n<N;n++)
	{
		if(llr_fn[n]>0) decoded_codeword[n]=1;
		else 0;
	}
//	show_array("hard-dec_2.txt",decoded_codeword,N);

	iteration = 0;
	whole_valid = check(H, decoded_codeword, parity_checks);

	//
//	show_array("hard-dec_3.txt",decoded_codeword,N);
	initprp_log_MS_without_llr_init();


	/* Do up to abs(max_iter) iterations of probability propagation, stopping
	   early if a codeword is found, unless max_iter is negative. */

	if(max_iter<0)
	{
		printf("Wrong! max_iter is not setted!\b");
		exit(1);
	}

	for (n = 0; n!=max_iter && n!=-max_iter && (max_iter<0 || !whole_valid); n++)
	{
		iteration = n+1;

		//fprintf(stderr, "\n iteration = %d", iteration);

		//input: encoded_data;
		iterprp_log_MS(H,llr_fn,decoded_codeword,(char*)0,(char*)0);
		//printf("%d code bit error.\n",code_bit_error_in_one_block);

		whole_valid = check(H, decoded_codeword, parity_checks);

        //if (whole_valid)
	    //{
		//	fprintf(stderr, "\nIteration = %d  decoding succeed! \n", iteration);
	    //}
	}

	//if (!whole_valid)
	//	tot_err_blks ++;

return whole_valid;
}

unsigned LDPC_mod2::count_info_error(char * decoded, char * source)
{//return the number of info error.
//input: decoded, info_length, N=block_len;
//output: info error count.

	int i;
	int code_len;
	code_len = N;


	int temp;
	temp=0;

	//mapping rule: the [i+M]-th code bit is the cols[i+M] info bit.
	for(i=M;i<N;i++)
	{
		if(decoded[cols[i]]!= source[i-M])
			temp++;
	}

	return temp;
}

unsigned LDPC_mod2::count_code_bit_error(char * decoded, char * encoded)
{
	int i,j;

	j=0;
	for(i=0;i<N;i++)
	{
		if(decoded[i]!=encoded[i])
			j++;
	}
return j;
}


// check ldpc_decode_ms() [11/18/2009]
void LDPC_mod2::check_ldpc_decoding()
{
	int i,j;

	initial_matrix_file();
	read_pchk(pchk_file);
	initial_decoding_space();
	read_gen(gen_file,0,0);

//	find_girth(H,&Min_Girth, &average_girth,0);

	//--
	source_data   = (char *)calloc(info_length,sizeof *source_data);
	encoded_data = (char *)calloc(N,sizeof *encoded_data);

	//random data:
	for(i=0;i<info_length;i++)
	{
		source_data[i]=(rand()%2);
	}
	sparse_encode(source_data,encoded_data,0);

	//--
	received_block = (double *)calloc(N,sizeof *received_block);

	//channel:
	for(i=0;i<N;i++)
		llr_fn[i] = encoded_data[i]=1? 10:-10;
	//input: llr_fn:

	ldpc_decode_MS_without_llr_init();
	return;
	ldpc_decode_F();

}

void LDPC_mod2::initprp_log_MS_without_llr_init()
{
//input:  received_block;
//output: llr, decoded_codeword;

	mod2entry *e;
	int n;
	int j;

	n = mod2sparse_cols(H);

	for (j = 0; j<n; j++)
	{
#ifdef sim_decoding_quantization
        //Sat and Quan for LDPC decoder
        llr_fn[j] = Sat_Quan(llr_fn[j],Finite_Lc_Plus_Max,Finite_Lc_Minus_Max,Finite_Lc_Tot,Finite_Lc_Frac);
#else
		//llr_fn[j] = llr_fn[j];
#endif
		for (e = mod2sparse_first_in_col(H,j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->pr = llr_fn[j];
		}

		// hard decision
		decoded_codeword[j] = (llr_fn[j] >= 0);
	}
}


void LDPC_mod2::show_array(char *fname,char *input, int len)
{
	int i;

	FILE *fp;
	fp=fopen(fname,"w");
	for(i=0;i<len;i++)
	{
		fprintf(fp,"%d\n",input[i]);
	}
	fclose(fp);

}

void LDPC_mod2::show_flash_array(char *fname,char *flash_array, int len)
{
	int i;

	FILE *fp;
	fp=fopen(fname,"w");
	for(i=0;i<len;i++)
	{
		if(flash_array[2*i]==0)
			fprintf(fp,"0\n0\n");
		else if(flash_array[2*i]==1)
			fprintf(fp,"0\n1\n");
		else if(flash_array[2*i]==2)
			fprintf(fp,"1\n1\n");
		else //if(flash_array[2*i]==0)
			fprintf(fp,"1\n0\n");
	}
	fclose(fp);
return;
}


void LDPC_mod2::show_array_float(char *fname,double*input, int len)
{
	int i;

	FILE *fp;
	fp=fopen(fname,"w");
	for(i=0;i<len;i++)
	{
		fprintf(fp,"%f\n",(float)input[i]);
	}
	fclose(fp);

}


unsigned LDPC_mod2::ldpc_decode_MS_with_defect_detection( int iter_defect_detect) // dgq on [1/19/2010]
{//this function works as MinSum decoding, without llr initial for AWGN channel.
//Suppose the input data to llr_fn is already given as llr.
//input is llr.
//The llr is detected to see whether input llr is completed opposite to
// llr from check equations.
// If yes, then to set this initial llr to be 0.

	//input: llr_fn; encoded_data;
	//output: decoded_codeword;

	int n,whole_valid;


//------------ defect detection first -----------------------
// added on [1/19/2010]
// detect defect: to see whether there is one initial LLR is inverse to all those
//		check equations' infor.
	check_initial_LLR_for_defect_detection(
		H,llr_fn,decoded_codeword,(char*)0,(char*)0, iter_defect_detect); //=> modifying llr_fn;



//------------start to LDPC decode --------------------------

	//init decoded_codeword:
	for(n=0;n<N;n++)
	{
		if(llr_fn[n]>0) decoded_codeword[n]=1;
		else 0;
	}
//	show_array("hard-dec_2.txt",decoded_codeword,N);

	iteration = 0;
	whole_valid = check(H, decoded_codeword, parity_checks);

	//
//	show_array("hard-dec_3.txt",decoded_codeword,N);
	initprp_log_MS_without_llr_init();


	/* Do up to abs(max_iter) iterations of probability propagation, stopping
	   early if a codeword is found, unless max_iter is negative. */

	if(max_iter<0)
	{
		printf("Wrong! max_iter is not setted!\b");
		exit(1);
	}

	for (n = 0; n!=max_iter && n!=-max_iter && (max_iter<0 || !whole_valid); n++)
	{
		iteration = n+1;

		//fprintf(stderr, "\n iteration = %d", iteration);

		//input: encoded_data;
		iterprp_log_MS(H,llr_fn,decoded_codeword,(char*)0,(char*)0);
		//printf("%d code bit error.\n",code_bit_error_in_one_block);

		whole_valid = check(H, decoded_codeword, parity_checks);

        //if (whole_valid)
	    //{
		//	fprintf(stderr, "\nIteration = %d  decoding succeed! \n", iteration);
	    //}
	}

	//if (!whole_valid)
	//	tot_err_blks ++;

return whole_valid;
}


void LDPC_mod2::check_initial_LLR_for_defect_detection( mod2sparse *H,      /* Parity check matrix */
                     float *lratio_log, /* initial Likelihood ratios log for bits */
                     char *dblk,         /* Place to store decoding */
                     char *look_row,     /* Which rows to look at, look at all if null */
                     char *look_col,     /* Which columns to look at, look at all if null */
					 int iter_for_dd
                   )
{
// added by dgq on [1/19/2010]: for defect detection initializing LLR:
//input: lratio_log:
//output: lratio_log.

	double min1;
	double min2;
	int min1_id;
	double sum;
	int sign;
	double abs_e_pr;
	mod2entry *e;
	int N, M;
	int i, j;
	int bitErrs;
	int flag1;
	int iter;

	M = mod2sparse_rows(H);
	N = mod2sparse_cols(H);

	//1, update a few iterations:
	for(iter=0;iter<iter_for_dd;iter++)
	{
		for (i = 0; i<M; i++)// for row.
		{
			if (look_row!=0 && !look_row[i]) continue;

			min1 = 100000;
			min2 = 100000;
			min1_id = 0;
			sign = 1;

			for (e = mod2sparse_first_in_row(H,i);
				 !mod2sparse_at_end(e);
				 e = mod2sparse_next_in_row(e)
				 )
				{
					if ( e->pr > 0 )
						abs_e_pr = e->pr;
					else
						abs_e_pr = -1 * e->pr;

					if ( abs_e_pr < min1 )
					{
						min2 = min1;
						min1 = abs_e_pr;
						min1_id = e->col;
					}
					else if ( abs_e_pr < min2 )
					{
						min2 = abs_e_pr;
					}

					sign *= ( e->pr > 0 ? 1 : -1 );
				}

				for (e = mod2sparse_first_in_row(H,i);
					!mod2sparse_at_end(e);
					e = mod2sparse_next_in_row(e)
					 )
				{//update for each entry:
					if ( min1_id == e->col )
						e->pr = sign*( e->pr > 0 ? 1 : -1 )* min2;
					else
						e->pr = sign*( e->pr > 0 ? 1 : -1 )* min1;
				}
		}//for i.

		//CVN UPDATE:
		for (j = 0; j<N; j++)
		{ if (look_col!=0 && !look_col[j]) continue;
			sum = lratio_log[j];
			for (e = mod2sparse_first_in_col(H,j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
				sum += e->pr;
				//sum = Sat_Quan(sum,Finite_Lc_Plus_Max,Finite_Lc_Minus_Max,Finite_Lc_Tot,Finite_Lc_Frac);
			 }
			//-----
			for (e = mod2sparse_first_in_col(H,j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
	#ifdef sim_decoding_quantization
				e->pr = Sat_Quan((sum - e->pr),Finite_Lc_Plus_Max,Finite_Lc_Minus_Max,Finite_Lc_Tot,Finite_Lc_Frac);
	#else
				e->pr = sum - e->pr;
	#endif
			}
		}
	}//


	//------------------------------------
	//2, to check whether initial llr is inverse to all check equation's output:
	for (j = 0; j<N; j++)//for each variable note:
	{ if (look_col!=0 && !look_col[j]) continue;

		sum = lratio_log[j];

		if(lratio_log[j]>0) flag1=1;
		else if(lratio_log[j]<0) flag1=-1;
		else flag1=0;

		for (e = mod2sparse_first_in_col(H,j);
			!mod2sparse_at_end(e);
			e = mod2sparse_next_in_col(e))
		{
			if(e->pr > 0 && flag1==-1) flag1=0; //inconsistency.
			if(e->pr < 0 && flag1==+1) flag1=0; //inconsistency.
		}
		if(flag1==0) //inconsistency means a defect, then to reset llr to be 0.
			lratio_log[j]=0;
	}

return;
}

float LDPC_mod2::absolute_value_f(float a)
{
	if(a>=0)
		return a;
	else
		return (a*(-1));

}

LDPC_mod2::~LDPC_mod2()
{
	int i;

	if(flag_G_binary_allocated==1)
	{
		for(i=0;i<Gb_rows;i++)
			free(G_binary[i]);
		free(G_binary);
	}

	if(flag_mapping_infor_to_code_allocated==1)
	{
		free(mapping_infor_to_code);
	}

	if(flag_permuted_H_allocated==1)
	{
		for(i=0;i<Ht_rows;i++)
			free(permuted_H[i]);
		free(permuted_H);
	}//if.

	if(flag_columnswaptracke_allocated==1)
		free(ColumnSwapTracker);
}

