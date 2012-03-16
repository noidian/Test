#include "LDPC_mod2_define.h"
#include "class_AWGN.h"

typedef unsigned long mod2word;	/* Data type that holds packed bits */

//make-gen.c
typedef enum { Sparse, Dense, Mixed } make_method;      /* Ways of making it */



//---- mod2sparse:
#define mod2sparse_first_in_row(m,i) ((m)->rows[i].right) /* Find the first   */ //????
#define mod2sparse_first_in_col(m,j) ((m)->cols[j].down)  /* or last entry in */
#define mod2sparse_last_in_row(m,i) ((m)->rows[i].left)   /* a row or column  */
#define mod2sparse_last_in_col(m,j) ((m)->cols[j].up)

#define mod2sparse_next_in_row(e) ((e)->right)  /* Move from one entry to     */
#define mod2sparse_next_in_col(e) ((e)->down)   /* another in any of the four */
#define mod2sparse_prev_in_row(e) ((e)->left)   /* possible directions        */
#define mod2sparse_prev_in_col(e) ((e)->up)

#define mod2sparse_at_end(e) ((e)->row<0) //mod2sparse_at_end(e)==1 means this entry does not exist./* See if we've reached the end. For the last element in row or column, the value of row and col will be -1     */

#define mod2sparse_row(e) ((e)->row)      /* Find out the row or column index */
#define mod2sparse_col(e) ((e)->col)      /* of an entry (indexes start at 0) */

#define mod2sparse_rows(m) ((m)->n_rows)  /* Get the number of rows or columns*/
#define mod2sparse_cols(m) ((m)->n_cols)  /* in a matrix                      */

//-- mod2dense:
#define mod2dense_rows(m) ((m)->n_rows)  /* Get the number of rows or columns */
#define mod2dense_cols(m) ((m)->n_cols)  /* in a matrix                       */

/* Extract the i'th bit of a mod2word. */

#define mod2_getbit(w,i) (((w)>>(i))&1)

/* Make a word like w, but with the i'th bit set to 1 (if it wasn't already). */

#define mod2_setbit1(w,i) ((w)|(1<<(i)))

/* Make a word like w, but with the i'th bit set to 0 (if it wasn't already). */

#define mod2_setbit0(w,i) ((w)&(~(1<<(i))))

enum decode_mode {BP, MS};

typedef enum{
	white, gray, black
}colar_type;

typedef struct node_entry /* Structure representing a node in the bipartite
			     graph */
{
  /* used for Breadth-first search */
  colar_type color;
  struct node_entry *parent;
  int depth;
  struct node_entry *next; /* used to keep track the queue in Breadth-first search  */
  int index;
  int group_ID;  /* used to find the first girth passing the root node
		    the ID of second level nodes */
} node_entry;

typedef struct node_vector
{
  int n;        /* length of vector */
  node_entry *entry;
} node_vector;


typedef enum
{ Mod2sparse_first,
  Mod2sparse_mincol,
  Mod2sparse_minprod
} mod2sparse_strategy;


typedef struct mod2entry /* Structure representing a non-zero entry, or
			      the header for a row or column               */
{
  int row, col;		  /* Row and column indexes of this entry, starting
                             at 0, and with -1 for a row or column header  */

  struct mod2entry *left, *right,  /* Pointers to entries adjacent in row  */
                   *up, *down;     /*   and column, or to headers.  Free   */
                                   /*   entries are linked by 'left'.      */

//#ifdef with_LDPC_DECODING
  double pr;	  /* Probability and likelihood ratios - not used  */
			           /*   by the mod2sparse module itself             */
  double lr;

  int pr_F;         /* finite representation of log likehood ratio */
  int NEG;
//#endif

#ifdef high_order_LDPC
  char ev; //entry value. Used only for GF(q) LDPC. added on 20091210
#endif

} mod2entry;  //Information of each entry in H.

typedef struct mod2block /* Block of entries allocated all at once */
{
  struct mod2block *next;  /* Next block that has been allocated */

  mod2entry entry[Mod2sparse_block]; /* Entries in this block */

} mod2block;

typedef struct		/* Representation of a sparse matrix */
{
  int n_rows;		  /* Number of rows in the matrix */
  int n_cols;		  /* Number of columns in the matrix */

  mod2entry *rows;	  /* Pointer to array of row headers */ //NOTE: this is actually not a entry, just a indicator, and its row and col value are always -1.

  mod2entry *cols;	  /* Pointer to array of column headers */ //NOTE: this is actually not a entry, just a indicator, and its row and col value are always -1.


  mod2block *blocks;	  /* Blocks that have been allocated */
  mod2entry *next_free;	  /* Next free entry */

} mod2sparse; //for mod2entry, row-connection is separate to column-connection.

//

#define Mod2sparse_block 10  /* Number of entries to block together for


/* STRUCTURE REPRESENTING A DENSE MATRIX.  These structures are dynamically
   allocated using mod2dense_allocate (or by other procedures that call
   mod2dense_allocate).  They should be freed with mod2dense_free when no
   longer required.

   Direct access to this structure should be avoided except in low-level
   routines.  Use the macros and procedures defined below instead. */

typedef struct
{
  int n_rows;		/* Number of rows in the matrix */
  int n_cols;		/* Number of columns in the matrix */

  int n_words;		/* Number of words used to store a column of bits */

  mod2word **col;	/* Pointer to array of pointers to columns */

  mod2word *bits;	/* Pointer to storage block for bits in this matrix
                           (pieces of this block are pointed to from col) */
} mod2dense;
                             //    memory allocation */







//----------------------------------------------------------------------------
// for making H matrix: [12/1/2009]
#define MinQueueSize ( 5 )
#define Error( Str )        FatalError( Str )
#define FatalError( Str )   fprintf( stderr, "%s\n", Str ), exit( 1 )

typedef int ElementType;
struct QueueRecord
{
   int Capacity;
   int Front;
   int Rear;
   int Size;
   ElementType *Array;
};
struct QueueRecord;
typedef struct QueueRecord *Queue;

//

#define max_depth 100 // maximum girth = 2*max_depth
//#define max_girth 200

#define MAX_CYCLE_LENGTH 8

// Structure representing a node in the bipartite graph

typedef struct tanner_graph  // bipartite graph
{
   int n_nodes;
   node_entry *check_node;
   node_entry *bit_node;
   mod2sparse *matrix;
} tanner_graph;


// structure to record nodes forming 4-cycle and shift values in base matrix
typedef struct cyc4_entry
{
	// for each 4-cycle ther are two check node: r1, r2, two bits node: c1, c2,
	// and 4 shift values: s1, s2, s3, s4
	int r1, r2, c1, c2;
	int *s1; // correspond to row r1 and col c1 position
	int *s2; // correspond to row r1 and col c2 position
	int *s3; // correspond to row r2 and col c1 position
	int *s4; // correspond to row r2 and col c2 position
} cyc4_entry;

// structure to record nodes forming 6-cycle and shift values in base matrix
typedef struct cyc6_entry
{
	// for each 6-cycle ther are 3 check node: r1, r2, r3; 3 bits node: c1, c2, c3
	// and 6 shift values: s1, s2, s3, s4, s5, s6
	int r1, r2, r3, c1, c2, c3;
	int *s1, *s2, *s3, *s4, *s5, *s6;
   int *ss1, *ss2, *ss3, *ss4, *ss5, *ss6;
} cyc6_entry;

// structure to record the sets of 4 dependant colomns
typedef struct weight4_entry
{
	int c[4];

} weight4_entry;


/* weight distribution pair */
typedef struct WD_pair
{
   int weight;
   int num;   //I guess this is the number of columns with this column weight.
} WD_pair;

typedef struct WD_vector
{
   int n;
   WD_pair *wd;
} WD_vector;
//end for making H matrix.
//----------------------------------------------------------------------------

class LDPC_mod2
{
public:
//
double Finite_Lc_Plus_Max; //pow(2,Finite_Lc_Tot-Finite_Lc_Frac-1)-pow(2,-Finite_Lc_Frac)
double Finite_Lc_Minus_Max; // -pow(2,Finite_Lc_Tot-Finite_Lc_Frac-1)



//---------------------------------------------------------------------------------
	//parameter:
//---------------------------------------------------------------------------------
	enum decode_mode dec_algorithm;
	char *pchk_file;   /* Parity check matrix file name */
	char *gen_file;    /* Generator matrix file name */

	char type;		/* Type of generator matrix representation (s/d/m) */

	int M;			/* Number of rows in parity check matrix */
	int N;			/* Number of columns in parity check matrix */
	//M and N are obtained in read_pchk() function.

	int flag_G_type; //0: sparse machine-readable G; 1: dense G; 2: binary G (mine).
	mod2sparse *H;		/* Parity check matrix */
	mod2dense *G;		/* Dense or mixed representation of generator matrix,
				   if type=='d' or type=='m' */

	mod2sparse *L, *U;	/* Sparse LU decomposition, if type=='s' */
	int *rows;		/* Ordering of rows in generator matrix (type 's') */

	int *cols;	//mapping: /* Ordering of columns in generator matrix */
				//for (j = M; j<N; j++)  encoded[cols[j]] = source[j-M];

	//cols[N]:
	int * mapping_code_bit_to_infor; //=cols. mapping_code_bit[i] is the position of coded[i] in infor data.

	//-added on 20100706am by Dong Guiqiang:
	char ** G_binary; //
	char flag_G_binary_allocated; //Flag showing G_binary space is allocated, thus need to free if over.
	int Gb_rows; //the row count of G_binary;
	int Gb_cols; //the column count of G_binary;
	int *mapping_infor_to_code; //coded[ mapping_infor_to_code[i] ] = infor[i];
	char flag_mapping_infor_to_code_allocated; //to show whether mapping_infor_to_cod[] is allocated.

	void print_mod2sparse_matrix_to_col_position_file(mod2sparse * h,char *filen);
	void read_binary_generator_matrix_Gb(char *filename, int row_gb , int col_gb);
	void check_binary_G_for_mapping_infor_to_code();
	void encoding_with_systemic_binary_G(char * source, char * encoded);

	//added on 20100706pm:
	char permuted_h_file[200];
	char columnswaptracker[200];
	char ** permuted_H; //the new H got through column and row swapping from original H.
	int * ColumnSwapTracker; //the ColumnSwapTracker[i] -th column in permuted_H is the i-th column of original H.
								//I.e. the ColumnSwapTracker[i]-th encoded bit is the i-th bit for original H.
	int size_T;
	int size_g;
	int Ht_rows; //size_g+size_T;
	int Ht_cols;

	char flag_permuted_H_allocated; //
	char flag_columnswaptracke_allocated;//
	void read_binary_permuted_H_for_RU_encoding(char *filename_ht, int size_t, int col_count);
//	void encoding_in_RU_method_through_permuted_H();
	void encoding_in_RU_method_through_permuted_H_with_g_as_0(char * source, char*encoded);
	void read_column_swap_tracker(char *filename, int col_count);




/* cols:
	//in current program, info_bits is mapped to the back end of coded sequence.
	//mapping : info[i] -> coded[ mapping[i+(block_len - info_len)].

	//cols is actually the mapping from info-bits to encoded bits:
	//the infor bits are at the back of encoded bits:
	//encoded block = [ parity bits | infor bits].

		//	for(i=0;i<info_len;i++)
		//		if(output[offset+mapping[i+(N-info_len)] ]!= info[i] ) error++;


		fp=fopen("3.Mapping_from_info_to_coded.txt","w");
		for (j = 0; j<N; j++)
		{
			if(j<M)
				fprintf(fp,"-\n");
			else
				fprintf(fp,"%d\n",cols[j]);
		}
*/


	int Min_Girth; //minimum girth;

	float info_rate;
	int info_length;
	int dimension_h; //the rank of H.

	int iteration; //recording current iteration number.
	int  max_iter;     /* maximum iteration number;
                      if max_iter<0, then doesn't use earlier stopping in the decoding */
	int flag_show_error_vs_iter; //1: show number of bit errors during each iteration.
	float * llr_fn;   //finite precision llr. log(p1/p0)  ::  1: llr>0; 0: llr<0;
	int *iter_distr; //the distribution of iteration.

	int code_bit_error_in_one_block; //
	int info_bit_error_in_one_block;
	long raw_error;
	float raw_ber;
	int flag_count_ib_error_distr; //1 for valid.
	int * cb_error_length_distr; //the distribution of count of coded bits errors.
	int * ib_error_length_distr; //the distribution of count of erratic infor bits.
	int max_recorded_ib;

	int flag_count_raw_be_for_successful_decoding; //1: to count the distribution of raw code bit errors for successful decoding;
	int flag_count_raw_be_for_failed_decoding; //1: to count the distribution of raw code bit errors for failed decoding;
	int * raw_cb_error_distr_for_successful_decoding;  //1: the distribution of raw code bit errors for successful decoding;
	int * raw_cb_error_distr_for_failed_decoding; //1: the distribution of raw code bit errors for failed decoding;
	int max_raw_cb_count;

	char **IterBitErr;  //[iter][j]: for getting the error pattern in each iteration;
	double	*softOut;   //for sova feedback info.
	double **IterSoft;  //[iter][j]: the j-th soft output in iter iteration;

	char *source_data;  //source data;
	char *encoded_data;		//encoded block;
	double *received_block;  //soft info. from SOVA detector.
	char *decoded_codeword; //the decode output;   // 1->1; 0 -> -1;		decoded_codeword[j] = (received_block[j] >= 0);

	char *parity_checks; //point for parity check result in LDPC.


//---------------------------------------------------------------------------------
	//function:
//---------------------------------------------------------------------------------
	LDPC_mod2();
	~LDPC_mod2();
	//basic function:
	void *  chk_alloc
		( unsigned n,		/* Number of elements */
		  unsigned size		/* Size of each element */
		);

//-------------------------------------------------------------------

	//mod2convert.c:
	void mod2sparse_to_dense (mod2sparse *, mod2dense *);
	void mod2dense_to_sparse (mod2dense *, mod2sparse *);

	//mod2sparse.c
	//static mod2entry *alloc_entry ( mod2sparse *m );
	int  intio_read  (FILE *);	/* Read an integer */
	void intio_write (FILE *, int);	/* Write an integer */
	mod2entry *alloc_entry ( mod2sparse *m );
	void initial_matrix_file();

	void init_quantization();
	mod2sparse *mod2sparse_allocate (int, int);
	void mod2sparse_free            (mod2sparse *);

	void mod2sparse_clear    (mod2sparse *);
	void mod2sparse_copy     (mod2sparse *, mod2sparse *);
	void mod2sparse_copycols (mod2sparse *, mod2sparse *, int *);

	void mod2sparse_print       (FILE *, mod2sparse *);
	int  mod2sparse_write       (FILE *, mod2sparse *);
	mod2sparse *mod2sparse_read (FILE *);

	mod2entry *mod2sparse_find   (mod2sparse *, int, int);
	mod2entry *mod2sparse_insert (mod2sparse *, int, int);
	void *mod2sparse_remove (mod2sparse *, int, int);
	void mod2sparse_delete       (mod2sparse *, mod2entry *);

	void mod2sparse_transpose (mod2sparse *, mod2sparse *);
	void mod2sparse_add       (mod2sparse *, mod2sparse *, mod2sparse *);
	void mod2sparse_multiply  (mod2sparse *, mod2sparse *, mod2sparse *);
	void mod2sparse_mulvec    (mod2sparse *, char *, char *);

	int mod2sparse_equal (mod2sparse *, mod2sparse *);

	int mod2sparse_count_row (mod2sparse *, int);
	int mod2sparse_count_col (mod2sparse *, int);

	void mod2sparse_add_row (mod2sparse *, int, mod2sparse *, int);
	void mod2sparse_add_col (mod2sparse *, int, mod2sparse *, int);

	int mod2sparse_decomp (mod2sparse *, int, mod2sparse *, mod2sparse *,
						   int *, int *, mod2sparse_strategy, int, int);

	int mod2sparse_forward_sub  (mod2sparse *, int *, char *, char *);
	int mod2sparse_backward_sub (mod2sparse *, int *, char *, char *);


//----------------------------------------------------------------------
	//mod2dense.c

	mod2dense *mod2dense_allocate (int, int);
	void mod2dense_free           (mod2dense *);

	void mod2dense_clear    (mod2dense *);
	void mod2dense_copy     (mod2dense *, mod2dense *);
	void mod2dense_copycols (mod2dense*, mod2dense *, int *);

	void mod2dense_print      (FILE *, mod2dense *);
	int  mod2dense_write      (FILE *, mod2dense *);
	mod2dense *mod2dense_read (FILE *);

	int  mod2dense_get (mod2dense *, int, int);
	void mod2dense_set (mod2dense *, int, int, int);
	int  mod2dense_flip(mod2dense *, int, int);

	void mod2dense_transpose (mod2dense *, mod2dense *);
	void mod2dense_add       (mod2dense *, mod2dense *, mod2dense *);
	void mod2dense_multiply  (mod2dense *, mod2dense *, mod2dense *);

	int mod2dense_equal (mod2dense *, mod2dense *);

	int mod2dense_invert          (mod2dense *, mod2dense *);
	int mod2dense_forcibly_invert (mod2dense *, mod2dense *, int *, int *);
	int mod2dense_invert_selected (mod2dense *, mod2dense *, int *);


//----------------------------------------------------------------------
	//-- graph.c
//----------------------------------------------------------------------
	/* check the number of connected component in the current graph */
	int Comp_Num (mod2sparse *, node_vector *, node_vector *);

	/* the Breadth-first search to find first cycle,
	   always start from subcode node */
	int BFS( mod2sparse *, node_vector *, node_vector *, int);

	/* find the girth distribution of a bipartite graph */
	void find_girth( mod2sparse *,int *,float *,FILE *);

	void init_node_vector( node_vector *v,int n);


//----------------------------------------------------------------------
	//LDPC encoding function:
//----------------------------------------------------------------------
    void sparse_encode (char * source, char * encoded, int array_start);
	void read_pchk( char *pchk_file);
	void read_gen
		( char *gen_file,	/* Name of generator matrix file */
		  int cols_only,	/* Read only column ordering? */
		  int no_pchk_file	/* No parity check file used? */
		);



//------------------------------------------------------------------------
	//LDPC decoding function:
//------------------------------------------------------------------------
	unsigned count_info_error(char * decoded, char * source);
	unsigned count_code_bit_error(char * decoded, char * encoded);

	void initprp_log_MS_without_llr_init(); //input: llr_fn; output: H;
	unsigned ldpc_decode_MS_without_llr_init(); //MinSum, without llr init func.

	void initial_decoding_space();
	unsigned ldpc_decode_MS();
	unsigned ldpc_decode_F();
//	void mod2sparse_mulvec ( mod2sparse *m,	/* The sparse matrix, with M rows and N columns */
//							  char *u,		/* The input vector, N long */
//							  char *v		/* Place to store the result, M long */
//							);


	/* DO ONE ITERATION OF PROBABILITY PROPAGATION --- LOG. */
	void iterprp_log_MS( mod2sparse *H,      /* Parity check matrix */
                     float *lratio_log, /* Likelihood ratios log for bits */
                     char *dblk,         /* Place to store decoding */
                     char *look_row,     /* Which rows to look at, look at all if null */
                     char *look_col      /* Which columns to look at, look at all if null */
                   );

	void initprp_log_MS(); /* INITIALIZE PROBABILITY PROPAGATION --- LOG.  Stores initial ratios, probabilities,
   and guess at decoding. */

	double Sat_Quan(double x, double Max_V, double Min_V, int m, int k);
	int Quantize(double x, int m, int k);
	void Tru2intS(int x, int m, int *dx);

	int check ( mod2sparse *H,	/* Parity check matrix */
            char *dblk,		/* Guess for codeword */
            char *pchk		/* Place to store parity checks */
          );

	unsigned ldpc_decode_BF();
	void check_ldpc_decoding();
	//--demo:
	void demo_all_process_with_machine_readable_format(); //This function demonstrates the process of using machine-readable H and machine-readable G to encode and decode.
	void demo_reading_binary_matrix_H_to_mod2sparse_and_systemic_G(); //with read in binary H and G.
	void demo_read_binary_H_and_encode_with_RU_method(); //this is to demo read in binary H, and encode with permuted H through RU method.

	void show_array(char *fname,char *input, int len);
	void show_array_float(char *fname,double *input, int len);
	void show_flash_array(char *fname,char *flash_array, int len);




//----- make-generator matrix:  -------
	// added on [11/30/2009]
	void make_dense_mixed (FILE *f, make_method, char *);     /* Procs to make it */
	void make_sparse (FILE *, mod2sparse_strategy, int, int);


	int make_gen_main
	(
	  char *argv1,
	  char *argv2,
	  char *argv3,
	  char *argv4,
	  char *argv5,
	  char *argv6,
	  char *argv7
	);
	void demo_make_gen();
//------------------------------------



//------------------------------------
// make generator matrix: [12/1/2009]
//	int init_addr[DMMB_NUM];
	int init_addr[144];

	int main_construct_H_matrix();
	void gen_quasi_ldpc_files(char *codefile,  /* File name for the LDPC codes to be constructed */
								  int filenum,     /* Number of LDPC codes to be constructed */
                          int b_rows,      /* Number of rows in base matrix */ //=>number of basic matrix in rows of H.
						  // base matrix is the mother matrix, with element as a small circulant matrix.
								  int b_cols,      /* Number of columns in base matrix */  //=>number of basic matrix in columns of H.
								  int sub_w,       /* Weight for the circulant */
								  int ex_factor);
	int quasi_cyclic ( mod2sparse *exPCH,      /* place to store the expanded ldpc matrix */
                    int rownum,             /* Number of rows in base matrix */
                    int colnum,             /* Number of columns in base matrix */
						  int sub_w,              /* weight of each circulant submatrix */
                    int ex_factor           /* expansion factor */
						 );
	void insert_sub_matrix( mod2sparse *mainMatrix,
                        int startRow,
                        int startCol,
                        int subDim,
                        int shift );
	void remove_sub_matrix( mod2sparse *mainMatrix,
                        int startRow,
                        int startCol,
                        int subDim,
                        int shift );
	int check_cycle(tanner_graph *g, int root_index, int length);
	int IsEmpty( Queue Q );
	int IsFull( Queue Q );

	ElementType FrontAndDequeue( Queue Q );
	void Dequeue( Queue Q );
	ElementType Front( Queue Q );

	void Enqueue( ElementType X, Queue Q );
	 int Succ( int Value, Queue Q );
	void DisposeQueue( Queue Q );
	void MakeEmpty( Queue Q );
	Queue CreateQueue( int MaxElements );

	FILE * open_file_std ( char *fname,	/* Name of file to open, or "-" for stdin/stdout */
                      char *mode	 );/* Mode for opening: eg, "r" or "w" */

	node_entry trace_back(node_entry u, int n);
	int demo_make_H_matrix();
//-------------------------------------------------------------------
	int decode_high_order_LDPC();

	void check_initial_LLR_for_defect_detection( mod2sparse *H,      /* Parity check matrix */
                     float *lratio_log, /* initial Likelihood ratios log for bits */
                     char *dblk,         /* Place to store decoding */
                     char *look_row,     /* Which rows to look at, look at all if null */
                     char *look_col,      /* Which columns to look at, look at all if null */
					 int iter_defect_detect
                   );
	unsigned ldpc_decode_MS_with_defect_detection(int iter_defect_detect);

	//20100630:
	float absolute_value_f(float a);
	double fel_func(double x);


	int Sum_Product_Algorithm( mod2sparse *hh,
									 char *encoded,
				                     char *decoded //output
									 );

	mod2sparse * read_binary_matrix_H_to_mod2sparse(char *filename, int n_col, int m_row);

	mod2sparse * read_H_positions_in_rows_to_mod2sparse(char *filename, int n_col, int m_row, int row_weight);

	int make_gen_without_reading_H
	(
	  char *argv1,  //parity check file => not used!
	  char *argv2,  //generator file.
	  char *argv3,  //method flag.
	  char *argv4,  //0.
	  char *argv5,  //0.
	  char *argv6,  //0.
	  char *argv7   //0.
	); //make generator matrix into L and U.

	//--
	void demo_dense_to_sparse();
	mod2dense mod2_read_binary_matrix_file_to_dense(char *filename);

	mod2sparse * read_RS_LDPC_base_matrix_to_sparse_matrix(
				char * input_filename,
				int b,
				char * output_filename,
				int row,
				int col
				);

	void trying_getting_trapping_set(mod2sparse *H, int th);

	//20100825:
	AWGN Gauss;
	float SNR_Eb_code;
	float SNR_Eb_info; // = SNR_Eb_code + 10 * log10(row_g*1.0/col_g);
	float awgn_sigma;

	void demo_AWGN_channel( );

	void demo_platform_for_performance_simulation();
	void demo_change_mod2sparse_to_row_position_file();

	//------------for generating H matrix:
	//for girth check:
	int * touched_bit_d1;
	int * touched_bit_d2;
	int * touched_bit_d3;
	int * touched_bit_d4;
	int * touched_bit_d5;

	void construct_H_matrix_try_1();
	int cycle_4_check(mod2sparse *hh);
	int cycle_6_check(mod2sparse *hh);
	int cycle_8_check(mod2sparse *hh);
	int cycle_10_check(mod2sparse *hh);
	int cycle_12_check(mod2sparse *hh);
	int detect_bits_with_distance1(int b, int *touched_bit, mod2sparse *hh);
	int row_containing_both_bits(mod2sparse *hh, int b1, int b2);


}; //end of LDPC_mod2.

/*
 	COM_PK.rblk => received_block;
	COM_PK.dblk,COM_PK.chks => decoded_codeword, parity_checks
	rblk_f => llr_fn
 */
