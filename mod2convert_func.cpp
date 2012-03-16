/* MOD2CONVERT.C - Routines converting between sparse and dense mod2 matrices.*/

/* Copyright (c) 1996 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "LDPC_mod2.h"

/* CONVERT A MOD2 MATRIX FROM SPARSE TO DENSE FORM.  The dense matrix 
   must already have been allocated, and must have at least as many
   rows and colums as the sparse matrix to be converted. */

void LDPC_mod2::mod2sparse_to_dense 
( mod2sparse *m, 	/* Sparse matrix to convert */
  mod2dense *r		/* Place to store result */
)
{
  mod2entry *e;
  int i;

  if (mod2sparse_rows(m)>mod2dense_rows(r) 
   || mod2sparse_cols(m)>mod2dense_cols(r))
  { fprintf(stderr,
      "mod2sparse_to_dense: Dimension of result matrix is less than source\n");
    exit(1);
  }

  mod2dense_clear(r);

  for (i = 0; i<mod2sparse_rows(m); i++)
  { e = mod2sparse_first_in_row(m,i);
    while (!mod2sparse_at_end(e))
    { mod2dense_set(r,i,mod2sparse_col(e),1);
      e = mod2sparse_next_in_row(e);
    }
  }
}


/* CONVERT A MOD2 MATRIX FROM DENSE TO SPARSE FORM.  The sparse matrix 
   must already have been allocated, and must have at least as many
   rows and colums as the dense matrix to be converted. */

void LDPC_mod2::mod2dense_to_sparse 
( mod2dense *m, 	/* Dense matrix to convert */
  mod2sparse *r		/* Place to store result */
)
{
  int i, j;

  if (mod2dense_rows(m)>mod2sparse_rows(r) 
   || mod2dense_cols(m)>mod2sparse_cols(r))
  { fprintf(stderr,
      "mod2dense_to_sparse: Dimension of result matrix is less than source\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2dense_rows(m); i++)
  { for (j = 0; j<mod2dense_cols(m); j++)
    { if (mod2dense_get(m,i,j))
      { mod2sparse_insert(r,i,j);
      }
    }
  }

return;
}
