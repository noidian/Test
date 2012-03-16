#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "LDPC_mod2.h"

//  [11/30/2009]

//this file is mainly from make-gen.c file;
/* MAKE-GEN.C - Make generator matrix from parity-check matrix. */
/* Copyright (c) 2000 by Radford M. Neal */



void usage(void);


/* MAIN PROGRAM. */

/*
int LDPC_mod2::make_gen_main
( int argc,
  char **argv
)
*/

void LDPC_mod2::demo_make_gen()
{
//	pchk_file="./H/matrix";
	pchk_file="./H/LDPC_2x24ex382w4_20100830pm_3.1";

	gen_file="./H/LDPC_2x24ex382w4_20100830pm_3.g";
	
	make_gen_main(pchk_file,gen_file,"sparse",0,0,0,0);
}

int LDPC_mod2::make_gen_main
( 
  char *argv1,  //parity check file.
  char *argv2,  //generator file.
  char *argv3,  //method flag.
  char *argv4,  //0.
  char *argv5,  //0.
  char *argv6,  //0.
  char *argv7   //0.
)
{
//NOTE that the input includes sparse_mode H matrix file, which is machine readable!!

  char *pchk_file, *gen_file, *other_gen_file;
  mod2sparse_strategy strategy;
  int abandon_when, abandon_number;
  make_method method;
  char *meth;
  char junk;
  FILE *f;

  /* Look at arguments. */

  if (!(pchk_file = argv1)
   || !(gen_file = argv2)
   || !(meth = argv3))
  { usage();
  }
  
  if (strcmp(meth,"sparse")==0)     
  { method = Sparse;
    strategy = Mod2sparse_minprod;
    abandon_number = 0;
    if (argv4)
    { if (strcmp(argv4,"first")==0)        strategy = Mod2sparse_first;
      else if (strcmp(argv4,"mincol")==0)  strategy = Mod2sparse_mincol;
      else if (strcmp(argv4,"minprod")==0) strategy = Mod2sparse_minprod;
      else 
      { usage();
      }
      if (argv5)
      { if (sscanf(argv5,"%d%c",&abandon_number,&junk)!=1 || abandon_number<=0
         || !argv6 
         || sscanf(argv6,"%d%c",&abandon_when,&junk)!=1 || abandon_when<=0
         || argv7)
        { usage();
        }
      }
    }
  }
  else if (strcmp(meth,"dense")==0) 
  { method = Dense;
    other_gen_file = argv4;
    if (other_gen_file && argv5)
    { usage();
    }
  }
  else if (strcmp(meth,"mixed")==0) 
  { method = Mixed;
    other_gen_file = argv4;
    if (other_gen_file && argv5)
    { usage();
    }
  }
  else 
  { usage();
  }

  /* Read parity check matrix. */

  read_pchk(pchk_file);

  if (N<=M)
  { fprintf(stderr,
     "Can't encode if number of bits (%d) isn't greater than number of checks (%d)\n",N,M);
    exit(1);
  }

  /* Create generator matrix file. */

  f = fopen(gen_file,"wb");
  if (f==NULL)
  { fprintf(stderr,"Can't create generator matrix file: %s\n",gen_file);
    exit(1);
  }

  /* Allocate space for column permutation. */

  cols = (int *)calloc (N, sizeof *cols);

  /* Create generator matrix with specified method. */

  switch (method)
  { case Sparse: 
    { make_sparse(f,strategy,abandon_number,abandon_when); 
      break;
    }
    case Dense: case Mixed:
    { make_dense_mixed(f,method,other_gen_file);
      break;
    }
    default: abort();
  }

  /* Check for error writing file. */

  if (ferror(f) || fclose(f)!=0)
  { fprintf(stderr,"Error writing to generator matrix file\n");
    exit(1);
  }

  return 0;
}


int LDPC_mod2::make_gen_without_reading_H
( 
  char *argv1,  //parity check file => not used!
  char *argv2,  //generator file.
  char *argv3,  //method flag.
  char *argv4,  //0.
  char *argv5,  //0.
  char *argv6,  //0.
  char *argv7   //0.
)
{
//NOTE that the input includes sparse_mode H matrix file, which is machine readable!!
//

  char *pchk_file, *gen_file, *other_gen_file;
  mod2sparse_strategy strategy;
  int abandon_when, abandon_number;
  make_method method;
  char *meth;
  char junk;
  FILE *f;

  /* Look at arguments. */

  if (!(pchk_file = argv1)
   || !(gen_file = argv2)
   || !(meth = argv3))
  { usage();
  }
  
  if (strcmp(meth,"sparse")==0)     
  { method = Sparse;
    strategy = Mod2sparse_minprod;
    abandon_number = 0;
    if (argv4)
    { if (strcmp(argv4,"first")==0)        strategy = Mod2sparse_first;
      else if (strcmp(argv4,"mincol")==0)  strategy = Mod2sparse_mincol;
      else if (strcmp(argv4,"minprod")==0) strategy = Mod2sparse_minprod;
      else 
      { usage();
      }
      if (argv5)
      { if (sscanf(argv5,"%d%c",&abandon_number,&junk)!=1 || abandon_number<=0
         || !argv6 
         || sscanf(argv6,"%d%c",&abandon_when,&junk)!=1 || abandon_when<=0
         || argv7)
        { usage();
        }
      }
    }
  }
  else if (strcmp(meth,"dense")==0) 
  { method = Dense;
    other_gen_file = argv4;
    if (other_gen_file && argv5)
    { usage();
    }
  }
  else if (strcmp(meth,"mixed")==0) 
  { method = Mixed;
    other_gen_file = argv4;
    if (other_gen_file && argv5)
    { usage();
    }
  }
  else 
  { usage();
  }

  /* No need to Read parity check matrix. It is already read in outside of this function.*/
  //  read_pchk(pchk_file);

  if (N<=M)
  { fprintf(stderr,
     "Can't encode if number of bits (%d) isn't greater than number of checks (%d)\n",N,M);
    exit(1);
  }

  /* Create generator matrix file. */

  f = fopen(gen_file,"wb");
  if (f==NULL)
  { fprintf(stderr,"Can't create generator matrix file: %s\n",gen_file);
    exit(1);
  }

  /* Allocate space for column permutation. */

  cols = (int *)calloc (N, sizeof *cols);

  /* Create generator matrix with specified method. */

  switch (method)
  { case Sparse: 
    { make_sparse(f,strategy,abandon_number,abandon_when); 
      break;
    }
    case Dense: case Mixed:
    { make_dense_mixed(f,method,other_gen_file);
      break;
    }
    default: abort();
  }

  /* Check for error writing file. */

  if (ferror(f) || fclose(f)!=0)
  { fprintf(stderr,"Error writing to generator matrix file\n");
    exit(1);
  }

  return 0;
}



/* MAKE DENSE OR MIXED REPRESENTATION OF GENERATOR MATRIX. */

void LDPC_mod2::make_dense_mixed
( FILE *f,
  make_method method,
  char *other_gen_file
)
{ 
  mod2dense *DH, *A, *AI, *B;
  int i, j, c, c2, n;

  DH = mod2dense_allocate(M,N);
  AI = mod2dense_allocate(M,M);
  B  = mod2dense_allocate(M,N-M);
  G  = mod2dense_allocate(M,N-M);

  mod2sparse_to_dense(H,DH);

  /* If another generator matrix was specified, invert using the set of
     columns it specifies. */

  if (other_gen_file)
  { 
    read_gen(other_gen_file,1,0);

    A  = mod2dense_allocate(M,M);
    mod2dense_copycols(DH,A,cols);

    if (!mod2dense_invert(A,AI))
    { fprintf(stderr,
       "Couldn't invert sub-matrix with column order given in other file\n");
      exit(1);
    }
  }

  /* If no other generator matrix was specified, invert using whatever 
     selection of columns is needed to get a non-singular sub-matrix. */

  if (!other_gen_file)
  {
    A  = mod2dense_allocate(M,N);

    n = mod2dense_invert_selected(DH,A,cols);
    mod2sparse_to_dense(H,DH);  /* DH was destroyed by invert_selected */

    if (n>0)
    { fprintf(stderr,"Note: Parity check matrix has %d redundant checks\n",n);
    }

    mod2dense_copycols(A,AI,cols);
  }

  mod2dense_copycols(DH,B,cols+M);

  /* Form final generator matrix. */

  if (method==Dense) 
  { mod2dense_multiply(AI,B,G);
  }
  else if (method==Mixed)
  { G = AI;
  }
  else
  { abort();
  }

  /* Compute and print number of 1s. */

  if (method==Dense)  
  { c = 0;
    for (i = 0; i<M; i++)
    { for (j = 0; j<N-M; j++)
      { c += mod2dense_get(G,i,j);
      }
    }
    fprintf(stderr,
      "Number of 1s per check in Inv(A) X B is %.1f\n", (double)c/M);
  }

  if (method==Mixed)
  { c = 0;
    for (i = 0; i<M; i++)
    { for (j = 0; j<M; j++)
      { c += mod2dense_get(G,i,j);
      }
    }
    c2 = 0;
    for (i = M; i<N; i++) 
    { c2 += mod2sparse_count_col(H,cols[i]);
    }
    fprintf(stderr,
     "Number of 1s per check in Inv(A) is %.1f, in B is %.1f, total is %.1f\n",
     (double)c/M, (double)c2/M, (double)(c+c2)/M);
  }

  /* Write the represention of the generator matrix to the file. */

  intio_write(f,('G'<<8)+0x80);

  if (method==Dense)      
  { fwrite ("d", 1, 1, f);
  }
  if (method==Mixed) 
  { fwrite ("m", 1, 1, f);
  }

  intio_write(f,M);
  intio_write(f,N);

  for (i = 0; i<N; i++) 
  { intio_write(f,cols[i]);
  }

  mod2dense_write (f, G);
}


/* MAKE SPARSE REPRESENTATION OF GENERATOR MATRIX. */

void LDPC_mod2::make_sparse
( FILE *f,
  mod2sparse_strategy strategy,
  int abandon_number,
  int abandon_when
)
{
  int n, cL, cU, cB;
  int i, j;

  /* Find LU decomposition. */

  rows = (int *)calloc(M, sizeof *rows);

  L = mod2sparse_allocate(M,M);
  U = mod2sparse_allocate(M,N);

  n = mod2sparse_decomp(H,M,L,U,rows,cols,strategy,abandon_number,abandon_when);

  if (n!=0 && abandon_number==0)
  { fprintf(stderr,"Note: Parity check matrix has %d redundant checks\n",n);
  }
  if (n!=0 && abandon_number>0)
  { fprintf(stderr,
  "Note: Have %d dependent columns, but this could be due to abandonment.\n",n);
    fprintf(stderr,
  "      Try again with lower abandonment number.\n");
    exit(1);
  }

  /* Compute and print number of 1s. */

  cL = cU = cB = 0;

  for (i = 0; i<M; i++) cL += mod2sparse_count_row(L,i);
  for (i = 0; i<M; i++) cU += mod2sparse_count_row(U,i);
  for (i = M; i<N; i++) cB += mod2sparse_count_col(H,cols[i]);

  fprintf(stderr,
   "Number of 1s per check in L is %.1f, U is %.1f, B is %.1f, total is %.1f\n",
    (double)cU/M, (double)cL/M, (double)cB/M, (double)(cL+cU+cB)/M);

  /* Write it all to the generator matrix file. */

  intio_write(f,('G'<<8)+0x80);

  fwrite ("s", 1, 1, f);

  intio_write(f,M);
  intio_write(f,N);

  for (i = 0; i<N; i++) 
  { intio_write(f,cols[i]);
  }

  for (i = 0; i<M; i++) 
  { intio_write(f,rows[i]);
  }

  mod2sparse_write (f, L);
  mod2sparse_write (f, U);
}


/* PRINT USAGE MESSAGE AND EXIT. */

void usage(void)
{ fprintf (stderr, 
   "Usage:  make-gen pchk-file gen-file method\n");
  fprintf (stderr, 
   "Method: sparse [ \"first\" | \"mincol\" | \"minprod\" ] [ abandon_num abandon_when ]\n");
  fprintf (stderr, 
   "    or: dense [ other-gen-file ]\n");
  fprintf (stderr, 
   "    or: mixed [ other-gen-file ]\n");
  exit(1);
}

void LDPC_mod2::print_mod2sparse_matrix_to_col_position_file(mod2sparse * h,char *filen)
{
	mod2entry *e;
	int i,j;

	FILE *fp;
	fp=fopen(filen,"w");
	for(j=0;j<M;j++)
	{
		for(e=mod2sparse_first_in_row(h,j);!mod2sparse_at_end(e);e=mod2sparse_next_in_row(e))
		{
			fprintf(fp,"%d\t",e->col);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

}