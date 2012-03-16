/* MOD2DENSE.C - Procedures for handling dense mod2 matrices. */

/* Copyright (c) 1996, 2000 by Radford M. Neal 
 *	modified by dgq
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

/* ALLOCATE SPACE FOR A DENSE MOD2 MATRIX.  Prints a message and exits if there 
   is not enough memory available for the matrix. */

mod2dense * LDPC_mod2::mod2dense_allocate 
( int n_rows, 		/* Number of rows in matrix */
  int n_cols		/* Number of columns in matrix */
)
{
  mod2dense *m;
  int j;

  if (n_rows<=0 || n_cols<=0)
  { fprintf(stderr,"mod2dense_allocate: Invalid number of rows or columns\n");
    exit(1);
  }

  m = (mod2dense *)calloc (1, sizeof *m);

  m->n_rows = n_rows;
  m->n_cols = n_cols;

  m->n_words = (n_rows+mod2_wordsize-1) >> mod2_wordsize_shift;

  m->col = (mod2word **)calloc (m->n_cols, sizeof *m->col);

  m->bits = (mod2word *)calloc(m->n_words*m->n_cols, sizeof *m->bits);

  for (j = 0; j<m->n_cols; j++)
  { m->col[j] = m->bits + j*m->n_words;
  }

  return m;
}


/* FREE SPACE OCCUPIED BY A DENSE MOD2 MATRIX.  The pointer passed should no
   longer be used after mod2dense_free returns. */

void LDPC_mod2::mod2dense_free
( mod2dense *m		/* Matrix to free */
)
{ free(m->bits);
  free(m->col);
  free(m);
}


/* CLEAR A DENSE MOD2 MATRIX.  Sets all the elements of the matrix to 0. */

void LDPC_mod2::mod2dense_clear
( mod2dense *r
)
{
  int k, j;

  for (j = 0; j<mod2dense_cols(r); j++)
  { for (k = 0; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }

}


/* COPY A DENSE MOD2 MATRIX.  The first matrix is copied to the second, which
   must already have been allocated, and must have dimensions at least as
   big as the first.  If the second matrix is larger, elements outside the
   boundaries of the first matrix are set to zero. */

void LDPC_mod2::mod2dense_copy
( mod2dense *m,		/* Matrix to copy */
  mod2dense *r		/* Place to store copy of matrix */
)
{ 
  int k, j;

  if (mod2dense_rows(m)>mod2dense_rows(r) 
   || mod2dense_cols(m)>mod2dense_cols(r))
  { fprintf(stderr,"mod2dense_copy: Destination matrix is too small\n");
    exit(1);
  }

  for (j = 0; j<mod2dense_cols(m); j++)
  { for (k = 0; k<m->n_words; k++)
    { r->col[j][k] = m->col[j][k];
    }
    for ( ; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }

  for ( ; j<mod2dense_cols(r); j++)
  { for (k = 0; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }
}


/* COPY COLUMNS OF A DENSE MOD2 MATRIX.  Copies selected columns of the 
   first matrix to the second matrix (which must already exist, and have 
   at least as many rows as the first matrix).  The columns to copy
   are given in order as an array of length the same as the number of 
   columns in the second matrix; duplicates are allowed. */

void LDPC_mod2::mod2dense_copycols
( mod2dense *m,		/* Matrix to copy */
  mod2dense *r,		/* Place to store copy of matrix */
  int *cols		/* Index of columns to copy, from 0 */
)
{ 
  int k, j;

  if (mod2dense_rows(m)>mod2dense_rows(r))
  { fprintf(stderr,
      "mod2dense_copycols: Destination matrix has fewer rows than source\n");
    exit(1);
  }

  for (j = 0; j<mod2dense_cols(r); j++)
  { if (cols[j]<0 || cols[j]>=mod2dense_cols(m))
    { fprintf(stderr,"mod2dense_copycols: Column index out of range\n");
      exit(1);
    }
    for (k = 0; k<m->n_words; k++)
    { r->col[j][k] = m->col[cols[j]][k];
    }
    for ( ; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }
}


/* PRINT A DENSE MOD2 MATRIX IN HUMAN-READABLE FORM.  The matrix is printed
   as "0" and "1" characters, with one line of "0"s and "1"s for each row 
   of the matrix. */

void LDPC_mod2::mod2dense_print     
( FILE *f,
  mod2dense *m
)
{ 
  int i, j;

  for (i = 0; i<mod2dense_rows(m); i++)
  { for (j = 0; j<mod2dense_cols(m); j++)
    { fprintf(f," %d",mod2dense_get(m,i,j));
    }
    fprintf(f,"\n");
  }
}


/* WRITE A DENSE MOD2 MATRIX TO A FILE IN MACHINE-READABLE FORM.  Returns 
   one if successful, zero if an error of some sort occurred. 

   The data written to the file consists of the number of rows and the
   number of columns (in machine-readable, not human-readable, form),
   followed by the bits in each column, packed into words. 

   Data is written using intio_write, so that it will be readable on a machine
   with a different byte-ordering.  At present, this assumes that the words 
   used to pack bits into are no longer than 32 bits. */

int LDPC_mod2::mod2dense_write     
( FILE *f, 
  mod2dense *m
)
{ 
  int j, k;

  intio_write(f,m->n_rows);
  if (ferror(f)) return 0;

  intio_write(f,m->n_cols);
  if (ferror(f)) return 0;

  for (j = 0; j<mod2dense_cols(m); j++)
  {
    for (k = 0; k<m->n_words; k++)
    { intio_write(f,m->col[j][k]);
      if (ferror(f)) return 0;
    }
  }

  return 1;
}


/* READ A DENSE MOD2 MATRIX STORED IN MACHINE-READABLE FORM FROM A FILE.  
   The data is expected to be in the format written by mod2dense_write.  The 
   value returned is zero if some sort of error occurred (either an error 
   reading the file, or data not in the right format).

   The file is left positioned after the last of the data read.  Other data
   (such as another matrix) may follow. */

mod2dense *LDPC_mod2::mod2dense_read  
( FILE *f
)
{ 
  int n_rows, n_cols;
  mod2dense *m;
  int j, k;
  
  n_rows = intio_read(f);
  if (feof(f) || ferror(f) || n_rows<=0) return 0;

  n_cols = intio_read(f);
  if (feof(f) || ferror(f) || n_cols<=0) return 0;

  m = mod2dense_allocate(n_rows,n_cols);

  for (j = 0; j<mod2dense_cols(m); j++)
  {
    for (k = 0; k<m->n_words; k++)
    { m->col[j][k] = intio_read(f);
      if (feof(f) || ferror(f)) 
      { mod2dense_free(m);
        return 0;
      }
    }
  }

  return m;
}


/* GET AN ELEMENT FROM A DENSE MOD2 MATRIX. */

int LDPC_mod2::mod2dense_get  
( mod2dense *m, 	/* Matrix to get element from */
  int row,		/* Row of element (starting with zero) */
  int col		/* Column of element (starting with zero) */
)
{
  if (row<0 || row>=mod2dense_rows(m) || col<0 || col>=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_get: row or column index out of bounds\n");
    exit(1);
  }

  return mod2_getbit (m->col[col][row>>mod2_wordsize_shift], 
                      row&mod2_wordsize_mask);
}


/* SET AN ELEMENT IN A DENSE MOD2 MATRIX. */

void LDPC_mod2::mod2dense_set 
( mod2dense *m, 	/* Matrix to modify element of */
  int row,		/* Row of element (starting with zero) */
  int col,		/* Column of element (starting with zero) */
  int value		/* New value of element (0 or 1) */
)
{ 
  mod2word *w;

  if (row<0 || row>=mod2dense_rows(m) || col<0 || col>=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_set: row or column index out of bounds\n");
    exit(1);
  }

  w = &m->col[col][row>>mod2_wordsize_shift];

  *w = value ? mod2_setbit1(*w,row&mod2_wordsize_mask) 
             : mod2_setbit0(*w,row&mod2_wordsize_mask);
}


/* FLIP AN ELEMENT OF A DENSE MOD2 MATRIX.  Returns the new value of the bit. */

int LDPC_mod2::mod2dense_flip  
( mod2dense *m, 	/* Matrix to flip element in */
  int row,		/* Row of element (starting with zero) */
  int col		/* Column of element (starting with zero) */
)
{
  mod2word *w;
  int b;

  if (row<0 || row>=mod2dense_rows(m) || col<0 || col>=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_flip: row or column index out of bounds\n");
    exit(1);
  }

  b = 1 ^ mod2_getbit (m->col[col][row>>mod2_wordsize_shift], 
                       row&mod2_wordsize_mask);

  w = &m->col[col][row>>mod2_wordsize_shift];

  *w = b ? mod2_setbit1(*w,row&mod2_wordsize_mask) 
         : mod2_setbit0(*w,row&mod2_wordsize_mask);

  return b;
}


/* COMPUTE THE TRANSPOSE OF A DENSE MOD2 MATRIX.  Stores the transpose of 
   the first matrix in the second matrix (which must already exist, and
   have as many rows as the first matrix has columns, and as many columns 
   as the first matrix has rows).

   The result matrix must not be the same as the operand. */

void LDPC_mod2::mod2dense_transpose
( mod2dense *m,		/* Matrix to compute transpose of (left unchanged) */
  mod2dense *r		/* Result of transpose operation */
)
{
  mod2word w, v, *p;
  int k1, j1, i2, j2;

  if (mod2dense_rows(m)!=mod2dense_cols(r) 
   || mod2dense_cols(m)!=mod2dense_rows(r))
  { fprintf(stderr,
     "mod2dense_transpose: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (r==m)
  { fprintf(stderr, 
     "mod2dense_transpose: Result matrix is the same as the operand\n");
    exit(1);
  }

  mod2dense_clear(r);

  for (j1 = 0; j1<mod2dense_cols(m); j1++)
  { 
    i2 = j1 >> mod2_wordsize_shift;
    v = 1 << (j1 & mod2_wordsize_mask);

    p = m->col[j1];
    k1 = 0;

    for (j2 = 0; j2<mod2dense_cols(r); j2++)
    { if (k1==0)
      { w = *p++;
        k1 = mod2_wordsize;
      }
      if (w&1)
      { r->col[j2][i2] |= v;
      }
      w >>= 1;     
      k1 -= 1;
    }
  }
}


/* ADD TWO DENSE MOD2 MATRICES.  Adds the first matrix and the second matrix, 
   storing the result in the third matrix.  Neither of the first two matrices
   is changed by this operation.  The three matrices must have the same numbers 
   of rows and columns.  

   It is permissible for the result matrix to be the same as one of the 
   operands. */

void LDPC_mod2::mod2dense_add
( mod2dense *m1,	/* Left operand of add */
  mod2dense *m2,	/* Right operand of add */
  mod2dense *r		/* Place to store result of add */
)
{
  int j, k;

  if (mod2dense_rows(m1)!=mod2dense_rows(r) 
   || mod2dense_cols(m1)!=mod2dense_cols(r) 
   || mod2dense_rows(m2)!=mod2dense_rows(r)
   || mod2dense_cols(m2)!=mod2dense_cols(r))
  { fprintf(stderr,"mod2dense_add: Matrices have different dimensions\n");
    exit(1);
  }

  for (j = 0; j<mod2dense_cols(r); j++)
  { for (k = 0; k<r->n_words; k++)
    { r->col[j][k] = m1->col[j][k] ^ m2->col[j][k];
    }
  }
}


/* MULTIPLY TWO DENSE MOD2 MATRICES.  Multiplies the first matrix by the second 
   matrix, storing the result in the third matrix.  Neither of the first 
   two matrices is changed by this operation.  The three matrices must have 
   compatible numbers of rows and columns. 

   The algorithm used runs faster if the second matrix (right operand of the
   multiply) is sparse, but it is also appropriate for dense matrices.  This
   procedure could be speeded up a bit by replacing the call of mod2dense_get
   with in-line code that avoids division, but this doesn't seem worthwhile
   at the moment. 

   The result matrix must not be the same as either operand. */

void LDPC_mod2::mod2dense_multiply 
( mod2dense *m1, 	/* Left operand of multiply */
  mod2dense *m2,	/* Right operand of multiply */
  mod2dense *r		/* Place to store result of multiply */
)
{
  int i, j, k;

  if (mod2dense_cols(m1)!=mod2dense_rows(m2) 
   || mod2dense_rows(m1)!=mod2dense_rows(r) 
   || mod2dense_cols(m2)!=mod2dense_cols(r))
  { fprintf(stderr,
     "mod2dense_multiply: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (r==m1 || r==m2)
  { fprintf(stderr,
      "mod2dense_multiply: Result matrix is the same as one of the operands\n");
    exit(1);
  }

  mod2dense_clear(r);

  for (j = 0; j<mod2dense_cols(r); j++)
  { for (i = 0; i<mod2dense_rows(m2); i++)
    { if (mod2dense_get(m2,i,j))
      { for (k = 0; k<r->n_words; k++)
        { r->col[j][k] ^= m1->col[i][k];
        }
      }
    }
  }
}


/* SEE WHETHER TWO DENSE MOD2 MATRICES ARE EQUAL.  Returns one if every 
   element of the first matrix is equal to the corresponding element of the
   second matrix.  The two matrices must have the same number of rows and the
   same number of columns. */

int LDPC_mod2::mod2dense_equal
( mod2dense *m1,
  mod2dense *m2
)
{
  int k, j, w;
  mod2word m;

  if (mod2dense_rows(m1)!=mod2dense_rows(m2) 
   || mod2dense_cols(m1)!=mod2dense_cols(m2))
  { fprintf(stderr,"mod2dense_equal: Matrices have different dimensions\n");
    exit(1);
  }

  w = m1->n_words;

  /* Form a mask that has 1s in the lower bit positions corresponding to
     bits that contain information in the last word of a matrix column. */

  m = (1 << (mod2_wordsize - (w*mod2_wordsize-m1->n_rows))) - 1;
  
  for (j = 0; j<mod2dense_cols(m1); j++)
  {
    for (k = 0; k<w-1; k++)
    { if (m1->col[j][k] != m2->col[j][k]) return 0;
    }

    if ((m1->col[j][k]&m) != (m2->col[j][k]&m)) return 0;
  }

  return 1;
}


/* INVERT A DENSE MOD2 MATRIX.  Inverts the first matrix passed, destroying 
   its contents in the process (though it remains a valid matrix for storing
   into later).  The matrix to be inverted must have the same number of
   rows as columns.  The inverse of this matrix is stored in the matrix 
   passed as the second argument (which must already exist, and have the 
   same number of rows and columns as the first).

   The value returned by mod2dense_invert is one if the inversion was 
   successful, zero if the matrix turned out to be singular (in which 
   case the contents of both the original matrix and the result matrix 
   will be garbage). 

   The result matrix must not be the same as the operand matrix. 

   The algorithm used is based on inverting A by transforming the equation
   AI = A to the equation AB = I using column operations, at which point B 
   is the inverse of A.  The representation of matrices used allows easy
   swapping of columns as needed by fiddling pointers.
*/

int LDPC_mod2::mod2dense_invert 
( mod2dense *m,		/* The matrix to find the inverse of (destroyed) */
  mod2dense *r		/* Place to store the inverse */
)
{
  mod2word *s, *t;
  int i, j, k, n, w, k0, b0;

  if (mod2dense_rows(m)!=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_invert: Matrix to invert is not square\n");
    exit(1);
  }

  if (r==m)
  { fprintf(stderr, 
      "mod2dense_invert: Result matrix is the same as the operand\n");
    exit(1);
  }

  n = mod2dense_rows(m);
  w = m->n_words;

  if (mod2dense_rows(r)!=n || mod2dense_cols(r)!=n)
  { fprintf(stderr,
     "mod2dense_invert: Matrix to receive inverse has wrong dimensions\n");
    exit(1);
  }

  mod2dense_clear(r);
  for (i = 0; i<n; i++) 
  { mod2dense_set(r,i,i,1);
  }

  for (i = 0; i<n; i++)
  { 
    k0 = i >> mod2_wordsize_shift;
    b0 = i & mod2_wordsize_mask;

    for (j = i; j<n; j++) 
    { if (mod2_getbit(m->col[j][k0],b0)) break;
    }

    if (j==n) return 0;

    if (j!=i)
    {
      t = m->col[i];
      m->col[i] = m->col[j];
      m->col[j] = t;

      t = r->col[i];
      r->col[i] = r->col[j];
      r->col[j] = t;
    }

    for (j = 0; j<n; j++)
    { if (j!=i && mod2_getbit(m->col[j][k0],b0))
      { s = m->col[j];
        t = m->col[i];
        for (k = k0; k<w; k++) s[k] ^= t[k];
        s = r->col[j];
        t = r->col[i];
        for (k = 0; k<w; k++) s[k] ^= t[k];
      }
    }
  }

  return 1;
}


/* INVERT A DENSE MOD2 MATRIX WITH COLUMNS SELECTED FROM A BIGGER MATRIX. 
   Inverts the matrix obtained by selecting certain columns from the first
   matrix passed (which must have at least as many columns as rows), storing
   the inverse in the second matrix passed.  The second matrix must already 
   exist, and must have the same number of rows and columns as the first 
   matrix.  The result is stored in the selected columns of the result matrix,  
   with the other columns being set to garbage. (Normally, one would extract 
   just the relevant columns afterwards using mod2dense_copycols.)  The 
   contents of the first matrix are destroyed (though it remains a valid 
   matrix for storing into later).

   Indexes of the columns selected are stored, in order, in cols, followed 
   by the columns not selected (in arbitrary order).  The value returned is 
   zero if an invertible matrix was found, and otherwise the number of
   columns left to find when inversion failed.  If inversion fails, the
   partial results are stored, with the selection of columns after that
   point being arbitrary (but valid).

   Note that when the first matrix is square, and non-singular, the result 
   is NOT in general the same as that obtained by calling mod2dense_invert, 
   which orders the columns of the inverse so that it applies to the original 
   ordering of the columns of the first matrix.

   The result matrix must not be the same as the operand matrix. 

   See mod2dense_invert for remarks on the algorithm.
*/

int LDPC_mod2::mod2dense_invert_selected 
( mod2dense *m,		/* Matrix from which to pick a submatrix to invert */
  mod2dense *r,		/* Place to store the inverse */
  int *cols		/* Set to indexes of columns used and not used */
)
{
  mod2word *s, *t;
  int i, j, k, n, n2, w, k0, b0, c, R;

  if (mod2dense_rows(m)>mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_invert_selected: Matrix has too few columns\n");
    exit(1);
  }

  if (r==m)
  { fprintf(stderr, 
      "mod2dense_invert_selected: Result matrix is the same as the operand\n");
    exit(1);
  }

  n = mod2dense_rows(m);
  w = m->n_words;

  n2 = mod2dense_cols(m);

  if (mod2dense_rows(r)!=n || mod2dense_cols(r)!=n2)
  { fprintf(stderr,
 "mod2dense_invert_selected: Matrix to receive inverse has wrong dimensions\n");
    exit(1);
  }

  mod2dense_clear(r);

  for (j = 0; j<n2; j++)
  { cols[j] = j;
  }

  R = 0;

  for (i = 0; i<n; i++)
  { 
    k0 = i >> mod2_wordsize_shift;
    b0 = i & mod2_wordsize_mask;

    for (j = i; j<n2; j++) 
    { if (mod2_getbit(m->col[cols[j]][k0],b0)) break;
    }

    if (j==n2) 
    { R += 1;
      continue;
    }

    c = cols[j];
    cols[j] = cols[i];
    cols[i] = c;

    mod2dense_set(r,i,c,1);

    for (j = 0; j<n2; j++)
    { if (j!=c && mod2_getbit(m->col[j][k0],b0))
      { s = m->col[j];
        t = m->col[c];
        for (k = k0; k<w; k++) s[k] ^= t[k];
        s = r->col[j];
        t = r->col[c];
        for (k = 0; k<w; k++) s[k] ^= t[k];
      }
    }
  }

  return R;
}


/* FORCIBLY INVERT A DENSE MOD2 MATRIX.  Inverts the first matrix passed, 
   destroying its contents in the process (though it remains a valid matrix 
   for storing into later).  The matrix to be inverted must have the same 
   number of rows as columns.  The inverse of this matrix is stored in the 
   matrix passed as the second argument (which must already exist, and have 
   the same number of rows and columns as the first).

   If the matrix to be inverted is singular, the inversion proceeds anyway,
   with bits of the matrix being changed as needed to create a non-singular
   matrix.  The value returned by mod2dense_invert is the number of elements 
   that had to be changed to make inversion possible (zero, if the original 
   matrix was non-singular). 

   The row and column indexes of the elements of the original matrix 
   that were changed are stored in the arrays passed as the last two
   elements.  These arrays must have as many elements as the dimension
   of the matrix.  (This is so even if it is known that fewer elements
   than this will be changed, as these arrays are also used as temporary
   storage by this routine.) 

   The result matrix must not be the same as the operand matrix. */

int LDPC_mod2::mod2dense_forcibly_invert 
( mod2dense *m, 	/* The matrix to find the inverse of (destroyed) */
  mod2dense *r,		/* Place to store the inverse */
  int *a_row,		/* Place to store row indexes of altered elements */
  int *a_col		/* Place to store column indexes of altered elements */
)
{
  mod2word *s, *t;
  int i, j, k, n, w, k0, b0;
  int u, c;

  if (mod2dense_rows(m)!=mod2dense_cols(m))
  { fprintf(stderr,
      "mod2dense_forcibly_invert: Matrix to invert is not square\n");
    exit(1);
  }

  if (r==m)
  { fprintf(stderr, 
      "mod2dense_forcibly_invert: Result matrix is the same as the operand\n");
    exit(1);
  }

  n = mod2dense_rows(m);
  w = m->n_words;

  if (mod2dense_rows(r)!=n || mod2dense_cols(r)!=n)
  { fprintf(stderr,
 "mod2dense_forcibly_invert: Matrix to receive inverse has wrong dimensions\n");
    exit(1);
  }

  mod2dense_clear(r);
  for (i = 0; i<n; i++) 
  { mod2dense_set(r,i,i,1);
  }

  for (i = 0; i<n; i++)
  { a_row[i] = -1;
    a_col[i] = i;
  }

  for (i = 0; i<n; i++)
  { 
    k0 = i >> mod2_wordsize_shift;
    b0 = i & mod2_wordsize_mask;

    for (j = i; j<n; j++) 
    { if (mod2_getbit(m->col[j][k0],b0)) break;
    }

    if (j==n)
    { j = i;
      mod2dense_set(m,i,j,1);
      a_row[i] = i;
    }

    if (j!=i)
    { 
      t = m->col[i];
      m->col[i] = m->col[j];
      m->col[j] = t;

      t = r->col[i];
      r->col[i] = r->col[j];
      r->col[j] = t;

      u = a_col[i];
      a_col[i] = a_col[j];
      a_col[j] = u;
    }

    for (j = 0; j<n; j++)
    { if (j!=i && mod2_getbit(m->col[j][k0],b0))
      { s = m->col[j];
        t = m->col[i];
        for (k = k0; k<w; k++) s[k] ^= t[k];
        s = r->col[j];
        t = r->col[i];
        for (k = 0; k<w; k++) s[k] ^= t[k];
      }
    }
  }

  c = 0;
  for (i = 0; i<n; i++)
  { if (a_row[i]!=-1)
    { a_row[c] = a_row[i];
      a_col[c] = a_col[i];
      c += 1;
    }
  }

  return c;
}
