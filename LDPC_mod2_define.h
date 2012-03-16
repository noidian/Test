//#define show_code_bit_error	

//#define sim_decoding_quantization

/* Decoding Algorithms LLR_BP or Min-Sum */
#define max_girth 200 /* maximum girth = 2*max_girth */
//#define with_LDPC_DECODING

//
#define Mod2sparse_block 10  /* Number of entries to block together for
                                memory allocation */

#define Finite_Lc_Tot 16
#define Finite_Lc_Frac 3

//-----------------------------------------------------------------------------
//mod2dense.c
#define mod2_wordsize 32	/* Number of bits that fit in a mod2word. Can't
				   be increased without changing intio module */

#define mod2_wordsize_shift 5	/* Amount to shift by to divide by wordsize */
#define mod2_wordsize_mask 0x1f /* What to and with to produce mod wordsize */
//-----------------------------------------------------------------------------

	//20100630pm:
#define fel_func_max_input 64 //the maximum input value of fel_func(). If input is larger than this, the output will be almost the same and too small.