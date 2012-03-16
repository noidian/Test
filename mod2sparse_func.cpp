/* MOD2SPARSE_func.Cpp - Procedures for handling sparse mod2 matrices. */

/* Copyright (c) 2000 by Radford M. Neal
 * modified a little by dgq on 20091113.

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


/* ALLOCATE AN ENTRY WITHIN A MATRIX. */

//static mod2entry * LDPC_mod2::alloc_entry
mod2entry * LDPC_mod2::alloc_entry
( mod2sparse *m
)
{
  mod2block *b;
  mod2entry *e;
  int k;

  if (m->next_free==0)
  {
    b = (mod2block *) calloc(1, sizeof (*b));

    b->next = m->blocks;
    m->blocks = b;

    for (k = 0; k<Mod2sparse_block; k++)
    { b->entry[k].left = m->next_free;
      m->next_free = &b->entry[k];
    }
  }

  e = m->next_free;
  m->next_free = e->left;

  e->pr = 0;
//  e->lr = 0;

  return e;
}


/* ALLOCATE SPACE FOR A SPARSE MOD2 MATRIX.  Prints a message and exits if
   there is not enough memory available for the matrix.  The matrix will
   be set to all zeros (ie, have no entries). */

mod2sparse *LDPC_mod2::mod2sparse_allocate
( int n_rows, 		/* Number of rows in matrix */
  int n_cols		/* Number of columns in matrix */
)
{
  mod2sparse *m;
  mod2entry *e;
  int i, j;

  if (n_rows<=0 || n_cols<=0)
  { fprintf(stderr,"mod2sparse_allocate: Invalid number of rows or columns\n");
    exit(1);
  }

  m = (mod2sparse *)calloc (1, sizeof *m);

  m->n_rows = n_rows;
  m->n_cols = n_cols;

  m->rows = (mod2entry *)calloc (n_rows, sizeof *m->rows);
  m->cols = (mod2entry *)calloc (n_cols, sizeof *m->cols);

  m->blocks = 0;
  m->next_free = 0;

  for (i = 0; i<n_rows; i++)
  { e = &m->rows[i];
    e->left = e->right = e->up = e->down = e;
    e->row = e->col = -1;
  }

  for (j = 0; j<n_cols; j++)
  { e = &m->cols[j];
    e->left = e->right = e->up = e->down = e;
    e->row = e->col = -1;
  }

  return m;
}


/* FREE SPACE OCCUPIED BY A SPARSE MOD2 MATRIX.  The pointer passed and all
   pointers to entries within this matrix should no longer be used after
   mod2sparse_free returns. */

void LDPC_mod2::mod2sparse_free
( mod2sparse *m		/* Matrix to free */
)
{
  mod2block *b;

  //before free m->rows and ->cols, those entries should be freed first...to do...20100705pm.
  free(m->rows);
  free(m->cols);

  while (m->blocks!=0)
  { b = m->blocks;
    m->blocks = b->next;
    free(b);
  }
}


/* CLEAR A SPARSE MATRIX TO ALL ZEROS.  The space occupied by entries in
   the matrix is freed for general use. */

void LDPC_mod2::mod2sparse_clear
( mod2sparse *r
)
{
  mod2block *b;
  mod2entry *e;
  int i, j;

  for (i = 0; i<mod2sparse_rows(r); i++)
  { e = &r->rows[i];
    e->left = e->right = e->up = e->down = e;
  }

  for (j = 0; j<mod2sparse_cols(r); j++)
  { e = &r->cols[j];
    e->left = e->right = e->up = e->down = e;
  }

  while (r->blocks!=0)
  { b = r->blocks;
    r->blocks = b->next;
    free(b);
  }
}


/* COPY A SPARSE MATRIX.  The first matrix is copied to the second, which
   must already have been allocated, and must have dimensions at least as
   big as the first.  The space occupied by entries in the second matrix is
   freed for general use (which may include being reused immediately for
   the copies of the entries in the first matrix). */

void LDPC_mod2::mod2sparse_copy
( mod2sparse *m,	/* Matrix to copy */
  mod2sparse *r		/* Place to store copy of matrix */
)
{
  mod2entry *e;
  int i;

  if (mod2sparse_rows(m)>mod2sparse_rows(r)
   || mod2sparse_cols(m)>mod2sparse_cols(r))
  { fprintf(stderr,"mod2sparse_copy: Destination matrix is too small\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2sparse_rows(m); i++)
  {
    e = mod2sparse_first_in_row(m,i);

    while (!mod2sparse_at_end(e))
    { mod2sparse_insert(r,e->row,e->col);
      e = mod2sparse_next_in_row(e);
    }
  }
}


/* COPY COLUMNS OF A SPARSE MOD2 MATRIX.  Copies selected columns of the
   first matrix to the second matrix (which must already exist, and have
   at least as many rows as the first matrix).  The columns to copy
   are given in order as an array of length the same as the number of
   columns in the second matrix; duplicates are allowed. */

void LDPC_mod2::mod2sparse_copycols
( mod2sparse *m,	/* Matrix to copy */
  mod2sparse *r,	/* Place to store copy of matrix */
  int *cols		/* Index of columns to copy, from 0 */
)
{
  mod2entry *e;
  int j;

  if (mod2sparse_rows(m)>mod2sparse_rows(r))
  { fprintf(stderr,
      "mod2sparse_copycols: Destination matrix has fewer rows than source\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (j = 0; j<mod2sparse_cols(r); j++)
  { if (cols[j]<0 || cols[j]>=mod2sparse_cols(m))
    { fprintf(stderr,"mod2sparse_copycols: Column index out of range\n");
      exit(1);
    }
    e = mod2sparse_first_in_col(m,cols[j]);
    while (!mod2sparse_at_end(e))
    { mod2sparse_insert(r,e->row,e->col);
      e = mod2sparse_next_in_col(e);
    }
  }
}


/* PRINT A SPARSE MOD2 MATRIX IN HUMAN-READABLE FORM.  The output consists of
   one line per row, of the form

      row: col col col ...

   where 'row' is the index of the row, and the 'col' entries are the indexes
   of columns that are non-zero in that row.  Indexes start at zero.  The
   number of columns is not indicated in the output. */

void LDPC_mod2::mod2sparse_print
( FILE *f,
  mod2sparse *m
)
{
  int rdigits, cdigits;
  mod2entry *e;
  int i;

  rdigits = mod2sparse_rows(m)<=10 ? 1
          : mod2sparse_rows(m)<=100 ? 2
          : mod2sparse_rows(m)<=1000 ? 3
          : mod2sparse_rows(m)<=10000 ? 4
          : mod2sparse_rows(m)<=100000 ? 5
          : 6;

  cdigits = mod2sparse_cols(m)<=10 ? 1
          : mod2sparse_cols(m)<=100 ? 2
          : mod2sparse_cols(m)<=1000 ? 3
          : mod2sparse_cols(m)<=10000 ? 4
          : mod2sparse_cols(m)<=100000 ? 5
          : 6;

  for (i = 0; i<mod2sparse_rows(m); i++)
  {
    fprintf(f,"%*d:",rdigits,i);

    e = mod2sparse_first_in_row(m,i);
    while (!mod2sparse_at_end(e))
    { fprintf(f," %*d",cdigits,mod2sparse_col(e));
      e = mod2sparse_next_in_row(e);
    }

    fprintf(f,"\n");
  }
}


/* WRITE A SPARSE MOD2 MATRIX TO A FILE IN MACHINE-READABLE FORM.  Returns
   one if successful, zero if an error of some sort occurred.

   The data written to the file is machine but not human readable.  It
   consists of negative integers giving row indexes (starting at 1), which
   apply until the next row index, and positive integers giving column
   indexes (starting at 1) for a non-zero entry in the matrix.  A zero index
   is written at the end.

   Data is written one byte at a time using intio_write so that it will
   be readable on a machine with a different byte ordering. */

int LDPC_mod2::mod2sparse_write
( FILE *f,
  mod2sparse *m
)
{
  mod2entry *e;
  int i;

  intio_write(f,m->n_rows);
  if (ferror(f)) return 0;

  intio_write(f,m->n_cols);
  if (ferror(f)) return 0;

  for (i = 0; i<mod2sparse_rows(m); i++)
  {
    e = mod2sparse_first_in_row(m,i);

    if (!mod2sparse_at_end(e))
    {
      intio_write (f, -(i+1));
      if (ferror(f)) return 0;

      while (!mod2sparse_at_end(e))
      {
        intio_write (f, mod2sparse_col(e)+1);
        if (ferror(f)) return 0;

        e = mod2sparse_next_in_row(e);
      }
    }
  }

  intio_write(f,0);
  if (ferror(f)) return 0;

  return 1;
}


/* READ A SPARSE MOD2 MATRIX STORED IN MACHINE-READABLE FORM FROM A FILE.
   The data is expected to be in the format written by mod2sparse_write.
   The value returned is zero if some sort of error occurred (either an
   error reading the file, or data not in the right format).

   The file is left positioned after the last of the data read.  Other data
   (such as another matrix) may follow. */

mod2sparse *LDPC_mod2::mod2sparse_read
( FILE *f
)
{
  int n_rows, n_cols;
  mod2sparse *m;
  int v, row, col;

  n_rows = intio_read(f);
  if (feof(f) || ferror(f) || n_rows<=0)
	  return 0;

  n_cols = intio_read(f);
  if (feof(f) || ferror(f) || n_cols<=0)
	  return 0;

  m = mod2sparse_allocate(n_rows,n_cols);

  row = -1;

  for (;;)
  {
    v = intio_read(f);
    if (feof(f) || ferror(f))
		break;

    if (v==0)
    { return m;
    }
    else if (v<0)
    { row = -v-1;
      if (row>=n_rows)
		  break;
    }
    else
    { col = v-1;
      if (col>=n_cols)
		  break;
      if (row==-1)
		  break;
      mod2sparse_insert(m,row,col);
    }
  }

  /* Error if we get here. */

  mod2sparse_free(m);
  return 0;
}


/* LOOK FOR AN ENTRY WITH GIVEN ROW AND COLUMN.  If the given element is
   non-zero, a pointer to its entry is returned.  If the element is zero
   (and hence has no entry), a zero pointer is returned.

   This procedure first looks to see if the entry is the last in its row or
   in its column (or would be if it existed).  If it isn't, it searches in
   parallel from the beginning of the row and the beginning of the column,
   so that the operation will be fast if either the row or the column is
   sparse. */

mod2entry *LDPC_mod2::mod2sparse_find
( mod2sparse *m,
  int row,
  int col
)
{
  mod2entry *re, *ce;

  if (row<0 || row>=mod2sparse_rows(m) || col<0 || col>=mod2sparse_cols(m))
  { fprintf(stderr,"mod2sparse_find: row or column index out of bounds\n");
    exit(1);
  }

  /* Check last entries in row and column. */

  re = mod2sparse_last_in_row(m,row);
  if (mod2sparse_at_end(re) || mod2sparse_col(re)<col)
  { return 0;
  }
  if (mod2sparse_col(re)==col)
  { return re;
  }

  ce = mod2sparse_last_in_col(m,col);
  if (mod2sparse_at_end(ce) || mod2sparse_row(ce)<row)
  { return 0;
  }
  if (mod2sparse_row(ce)==row)
  { return ce;
  }

  /* Search row and column in parallel, from the front. */

  re = mod2sparse_first_in_row(m,row);
  ce = mod2sparse_first_in_col(m,col);

  for (;;)
  {
    if (mod2sparse_at_end(re) || mod2sparse_col(re)>col)
    { return 0;
    }
    if (mod2sparse_col(re)==col)
    { return re;
    }

    if (mod2sparse_at_end(ce) || mod2sparse_row(ce)>row)
    { return 0;
    }
    if (mod2sparse_row(ce)==row)
    { return ce;
    }

    re = mod2sparse_next_in_row(re);
    ce = mod2sparse_next_in_col(ce);
  }
}


/* INSERT AN ENTRY WITH GIVEN ROW AND COLUMN.  Adds a new entry at the given
   row and column.  If such an entry already exists, nothing is done (this
   is not considered to be an error).

   This procedure first looks to see if the entry belongs at the end of the
   the row and/or column.  If it doesn't, it searches for its place in the
   row or column starting from the beginning.

   Inserting an entry may require allocating space from the general memory
   pool, which may fail, in which case a message is printed and the program
   terminated. */

mod2entry *LDPC_mod2::mod2sparse_insert
( mod2sparse *m,
  int row,
  int col
)
{
  mod2entry *re, *ce, *ne;

  if (row<0 || row>=mod2sparse_rows(m) || col<0 || col>=mod2sparse_cols(m))
  { fprintf(stderr,"mod2sparse_insert: row or column index out of bounds\n");
    exit(1);
  }

  /* Find old entry and return it, or allocate new entry and insert into row. */

  re = mod2sparse_last_in_row(m,row); //the last element in row-th row.

  if (!mod2sparse_at_end(re) && mod2sparse_col(re)==col) //if ((e)->row>0) => this is valid entry.
  {//if this entry already exist, then return;
	  return re; //
  }

  if (mod2sparse_at_end(re) || mod2sparse_col(re)<col) //if ((e)->row<0), this entry is not valid entry, and it should be the H->rows[i], which is not entry, but a indicator/pointer.
  {//if there is no element in current chain, then insert it.
	  re = re->right;
  }
  else
  {//if there are already entries in the chain, then find the position to insert the entry:
    re = mod2sparse_first_in_row(m,row);

    for (;;)
    {
      if (!mod2sparse_at_end(re) && mod2sparse_col(re)==col)
      { return re;
      }

      if (mod2sparse_at_end(re) || mod2sparse_col(re)>col) //find the first element with col larger than the input col. The new element should be added before this found element.
      { break;
      }

	  //if this element's col value smaller than expected, then go to check the next element to see whether its col larger than expected.
	  //this shows that ->right ->right is the order of col increasing, as common sense.
      re = mod2sparse_next_in_row(re);
    }
  }

  //if reach here, re is the entry on the right of inserted entry:
  //then insert a new entry to the left of re:
  ne = alloc_entry(m);

  ne->row = row;
  ne->col = col;

  ne->left = re->left;
  ne->right = re;
  ne->left->right = ne;
  ne->right->left = ne;

  /* Insert new entry into column.  If we find an existing entry here,
     the matrix must be garbled, since we didn't find it in the row. */

  ce = mod2sparse_last_in_col(m,col); //((m)->cols[j].up)

  if (!mod2sparse_at_end(ce) && mod2sparse_row(ce)==row)
  {//if ce is valid entry, i.e. not pointer, and its row already is expected, which means the entry already exists. If so, this function should return already.
	fprintf(stderr,"mod2sparse_insert: Garbled matrix\n");
	exit(1);
  }

  if (mod2sparse_at_end(ce) || mod2sparse_row(ce)<row)
  { //if ce is the pointer/indicator, then just insert this new entry:
	  ce = ce->down;
  }
  else
  {	//find the position to insert this new entry:
    ce = mod2sparse_first_in_col(m,col);

    for (;;)
    {
      if (!mod2sparse_at_end(ce) && mod2sparse_row(ce)==row)
      { //if this entry already exist, it should not reach here.
		fprintf(stderr,"mod2sparse_insert: Garbled matrix\n");
        exit(1);
      }

	  //to insert this new entry before the first entry with the row value larger than input row
      if (mod2sparse_at_end(ce) || mod2sparse_row(ce)>row)
      { break;
      }

      ce = mod2sparse_next_in_col(ce);
    }
  }

  ne->up = ce->up;
  ne->down = ce;
  ne->up->down = ne;
  ne->down->up = ne;

  /* Return the new entry. */

  return ne;
}


/* DELETE AN ENTRY FROM A SPARSE MATRIX.  Deletes a non-zero entry - ie,
   effectively sets the element of the matrix to zero.  The entry is freed
   for future use in the same matrix, but not (immediately, at least) for
   use in other matrices, or generally. */

void LDPC_mod2::mod2sparse_delete
( mod2sparse *m,
  mod2entry *e
)
{
  if (e==0)
  { fprintf(stderr,"mod2sparse_delete: Trying to delete a null entry\n");
    exit(1);
  }

  if (e->row<0 || e->col<0)
  { fprintf(stderr,"mod2sparse_delete: Trying to delete a header entry\n");
    exit(1);
  }

  e->left->right = e->right;
  e->right->left = e->left;

  e->up->down = e->down;
  e->down->up = e->up;

  e->left = m->next_free;
  m->next_free = e;
}


/* TEST WHETHER TWO SPARSE MATRICES ARE EQUAL.  Returns one if they are,
   zero if they aren't.  The matrices must have the same dimensions. */

int LDPC_mod2::mod2sparse_equal
( mod2sparse *m1,
  mod2sparse *m2
)
{
  mod2entry *e1, *e2;
  int i;

  if (mod2sparse_rows(m1)!=mod2sparse_rows(m2)
   || mod2sparse_cols(m1)!=mod2sparse_cols(m2))
  { fprintf(stderr,"mod2sparse_equal: Matrices have different dimensions\n");
    exit(1);
  }

  for (i = 0; i<mod2sparse_rows(m1); i++)
  {
    e1 = mod2sparse_first_in_row(m1,i);
    e2 = mod2sparse_first_in_row(m2,i);

    while (!mod2sparse_at_end(e1) && !mod2sparse_at_end(e2))
    {
      if (mod2sparse_col(e1)!=mod2sparse_col(e2))
      { return 0;
      }

      e1 = mod2sparse_next_in_row(e1);
      e2 = mod2sparse_next_in_row(e2);
    }

    if (!mod2sparse_at_end(e1) || !mod2sparse_at_end(e2))
    { return 0;
    }
  }

  return 1;
}


/* COMPUTE THE TRANSPOSE OF A SPARSE MOD2 MATRIX.  Stores the transpose
   of the first matrix in the second matrix (which must already exist, and
   have as many rows as the first matrix has columns, and as many
   columns as the first matrix has rows).

   The result matrix must not be the same as the operand. */

void LDPC_mod2::mod2sparse_transpose
( mod2sparse *m,	/* Matrix to compute transpose of (left unchanged) */
  mod2sparse *r		/* Result of transpose operation */
)
{
  mod2entry *e;
  int i;

  if (mod2sparse_rows(m)!=mod2sparse_cols(r)
   || mod2sparse_cols(m)!=mod2sparse_rows(r))
  { fprintf(stderr,
     "mod2sparse_transpose: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (r==m)
  { fprintf(stderr,
     "mod2sparse_transpose: Result matrix is the same as the operand\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2sparse_rows(m); i++)
  {
    e = mod2sparse_first_in_row(m,i);

    while (!mod2sparse_at_end(e))
    { mod2sparse_insert(r,mod2sparse_col(e),i);
      e = mod2sparse_next_in_row(e);
    }
  }
}


/* ADD TWO SPARSE MOD2 MATRICES.  Adds the first matrix by the second matrix,
   storing the result in the third matrix.  Neither of the first two matrices
   is changed by this operation.  The three matrices must have the same numbers
   of rows and columns.

   The result matrix must not be the same as either operand. */

void LDPC_mod2::mod2sparse_add
( mod2sparse *m1,	/* Left operand of add */
  mod2sparse *m2,	/* Right operand of add */
  mod2sparse *r		/* Place to store result of add */
)
{
  mod2entry *e1, *e2;
  int i;

  if (mod2sparse_rows(m1)!=mod2sparse_rows(r)
   || mod2sparse_cols(m1)!=mod2sparse_cols(r)
   || mod2sparse_rows(m2)!=mod2sparse_rows(r)
   || mod2sparse_cols(m2)!=mod2sparse_cols(r))
  { fprintf(stderr,"mod2sparse_add: Matrices have different dimensions\n");
    exit(1);
  }

  if (r==m1 || r==m2)
  { fprintf(stderr,
     "mod2sparse_add: Result matrix is the same as one of the operands\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2sparse_rows(r); i++)
  {
    e1 = mod2sparse_first_in_row(m1,i);
    e2 = mod2sparse_first_in_row(m2,i);

    while (!mod2sparse_at_end(e1) && !mod2sparse_at_end(e2))
    {
      if (mod2sparse_col(e1)==mod2sparse_col(e2))
      { e1 = mod2sparse_next_in_row(e1);
        e2 = mod2sparse_next_in_row(e2);
      }

      else if (mod2sparse_col(e1)<mod2sparse_col(e2))
      { mod2sparse_insert(r,i,mod2sparse_col(e1));
        e1 = mod2sparse_next_in_row(e1);
      }

      else
      { mod2sparse_insert(r,i,mod2sparse_col(e2));
        e2 = mod2sparse_next_in_row(e2);
      }
    }

    while (!mod2sparse_at_end(e1))
    { mod2sparse_insert(r,i,mod2sparse_col(e1));
      e1 = mod2sparse_next_in_row(e1);
    }

    while (!mod2sparse_at_end(e2))
    { mod2sparse_insert(r,i,mod2sparse_col(e2));
      e2 = mod2sparse_next_in_row(e2);
    }
  }
}


/* MULTIPLY TWO SPARSE MOD2 MATRICES.  Multiplies the first matrix by the
   second matrix, storing the result in the third matrix.  Neither of the first
   two matrices is changed by this operation.  The previous contents of the
   third matrix are destroyed, and its entries release for general reuse.  The
   three matrices must have compatible numbers of rows and columns.

   The result matrix must not be the same as either operand. */

void LDPC_mod2::mod2sparse_multiply
( mod2sparse *m1, 	/* Left operand of multiply */
  mod2sparse *m2,	/* Right operand of multiply */
  mod2sparse *r		/* Place to store result of multiply */
)
{
  mod2entry *e1, *e2;
  int i, j, b;

  if (mod2sparse_cols(m1)!=mod2sparse_rows(m2)
   || mod2sparse_rows(m1)!=mod2sparse_rows(r)
   || mod2sparse_cols(m2)!=mod2sparse_cols(r))
  { fprintf (stderr,
      "mod2sparse_multiply: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (r==m1 || r==m2)
  { fprintf(stderr,
     "mod2sparse_multiply: Result matrix is the same as one of the operands\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2sparse_rows(m1); i++)
  {
    if (mod2sparse_at_end(mod2sparse_first_in_row(m1,i)))
    { continue;
    }

    for (j = 0; j<mod2sparse_cols(m2); j++)
    {
      b = 0;

      e1 = mod2sparse_first_in_row(m1,i);
      e2 = mod2sparse_first_in_col(m2,j);

      while (!mod2sparse_at_end(e1) && !mod2sparse_at_end(e2))
      {
        if (mod2sparse_col(e1)==mod2sparse_row(e2))
        { b ^= 1;
          e1 = mod2sparse_next_in_row(e1);
          e2 = mod2sparse_next_in_col(e2);
        }

        else if (mod2sparse_col(e1)<mod2sparse_row(e2))
        { e1 = mod2sparse_next_in_row(e1);
        }

        else
        { e2 = mod2sparse_next_in_col(e2);
        }
      }

      if (b)
      { mod2sparse_insert(r,i,j);
      }
    }
  }
}


/* MULTIPLY VECTOR BY SPARSE MATRIX.  Multiplies a column vector, unpacked
   one bit per byte, by a sparse matrix.  The result is stored in another
   unpacked vector, which must not overlap the input vector. */

void LDPC_mod2::mod2sparse_mulvec
( mod2sparse *m,	/* The sparse matrix, with M rows and N columns */
  char *u,		/* The input vector, N long */
  char *v		/* Place to store the result, M long */
)
{
  mod2entry *e;
  int M, N;
  int i, j;

  M = mod2sparse_rows(m);
  N = mod2sparse_cols(m);

  for (i = 0; i<M; i++) v[i] = 0;

  for (j = 0; j<N; j++)
  { if (u[j]==1)
    { for (e = mod2sparse_first_in_col(m,j);
           !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e))
      { v[mod2sparse_row(e)] ^= 1;
      }
    }
  }
}


/* COUNT ENTRIES IN A ROW. */

int LDPC_mod2::mod2sparse_count_row
( mod2sparse *m,
  int row
)
{
  mod2entry *e;
  int count;

  if (row<0 || row>=mod2sparse_rows(m))
  { fprintf(stderr,"mod2sparse_count_row: row index out of bounds\n");
    exit(1);
  }

  count = 0;

  for (e = mod2sparse_first_in_row(m,row);
       !mod2sparse_at_end(e);
       e = mod2sparse_next_in_row(e))
  { count += 1;
  }

  return count;
}


/* COUNT ENTRIES IN A COLUMN. */

int LDPC_mod2::mod2sparse_count_col
( mod2sparse *m,
  int col
)
{
  mod2entry *e;
  int count;

  if (col<0 || col>=mod2sparse_cols(m))
  { fprintf(stderr,"mod2sparse_count_col: column index out of bounds\n");
    exit(1);
  }

  count = 0;

  for (e = mod2sparse_first_in_col(m,col);
       !mod2sparse_at_end(e);
       e = mod2sparse_next_in_col(e))
  { count += 1;
  }

  return count;
}


/* ADD TO A ROW. */

void LDPC_mod2::mod2sparse_add_row
( mod2sparse *m1,	/* Matrix containing row to add to */
  int row1,		/* Index in this matrix of row to add to */
  mod2sparse *m2,	/* Matrix containing row to add from */
  int row2		/* Index in this matrix of row to add from */
)
{
  mod2entry *f1, *f2, *ft;

  if (mod2sparse_cols(m1)<mod2sparse_cols(m2))
  { fprintf (stderr,
     "mod2sparse_add_row: row added to is shorter than row added from\n");
    exit(1);
  }

  if (row1<0 || row1>=mod2sparse_rows(m1)
   || row2<0 || row2>=mod2sparse_rows(m2))
  { fprintf (stderr,"mod2sparse_add_row: row index out of range\n");
    exit(1);
  }

  f1 = mod2sparse_first_in_row(m1,row1);
  f2 = mod2sparse_first_in_row(m2,row2);

  while (!mod2sparse_at_end(f1) && !mod2sparse_at_end(f2))
  { if (mod2sparse_col(f1)>mod2sparse_col(f2))
    { mod2sparse_insert(m1,row1,mod2sparse_col(f2));
      f2 = mod2sparse_next_in_row(f2);
    }
    else
    { ft = mod2sparse_next_in_row(f1);
      if (mod2sparse_col(f1)==mod2sparse_col(f2))
      { mod2sparse_delete(m1,f1);
        f2 = mod2sparse_next_in_row(f2);
      }
      f1 = ft;
    }
  }

  while (!mod2sparse_at_end(f2))
  { mod2sparse_insert(m1,row1,mod2sparse_col(f2));
    f2 = mod2sparse_next_in_row(f2);
  }
}


/* ADD TO A COLUMN. */

void LDPC_mod2::mod2sparse_add_col
( mod2sparse *m1,	/* Matrix containing column to add to */
  int col1,		/* Index in this matrix of column to add to */
  mod2sparse *m2,	/* Matrix containing column to add from */
  int col2		/* Index in this matrix of column to add from */
)
{
  mod2entry *f1, *f2, *ft;

  if (mod2sparse_rows(m1)<mod2sparse_rows(m2))
  { fprintf (stderr,
     "mod2sparse_add_col: Column added to is shorter than column added from\n");
    exit(1);
  }

  if (col1<0 || col1>=mod2sparse_cols(m1)
   || col2<0 || col2>=mod2sparse_cols(m2))
  { fprintf (stderr,"mod2sparse_add_col: Column index out of range\n");
    exit(1);
  }

  f1 = mod2sparse_first_in_col(m1,col1);
  f2 = mod2sparse_first_in_col(m2,col2);

  while (!mod2sparse_at_end(f1) && !mod2sparse_at_end(f2))
  { if (mod2sparse_row(f1)>mod2sparse_row(f2))
    { mod2sparse_insert(m1,mod2sparse_row(f2),col1);
      f2 = mod2sparse_next_in_col(f2);
    }
    else
    { ft = mod2sparse_next_in_col(f1);
      if (mod2sparse_row(f1)==mod2sparse_row(f2))
      { mod2sparse_delete(m1,f1);
        f2 = mod2sparse_next_in_col(f2);
      }
      f1 = ft;
    }
  }

  while (!mod2sparse_at_end(f2))
  { mod2sparse_insert(m1,mod2sparse_row(f2),col1);
    f2 = mod2sparse_next_in_col(f2);
  }
}


/* FIND AN LU DECOMPOSITION OF A SPARSE MATRIX.  Takes as input a matrix, A,
   having M rows and N columns, and an integer K.  Finds an LU decomposition
   of a K by K sub-matrix of A.  The decomposition is stored in the matrix L,
   with M rows and K columns, and the matrix U, with K rows and N columns.
   The product of L and U will be equal to the K by K submatrix of A obtained
   by taking only rows and columns that are given in the first K elements of
   the 'rows' and 'cols' arrays, which are set by this procedure, with this
   sub-matrix distributed over the original M rows and N columns.  Furthermore,
   the ordering of the row and column indexes in these arrays will be set so
   that if the rows of L and the columns of U were rearranged in this order, L
   would be lower triangular, with zeros in rows past row K, and U would be
   upper triangular, with zeros in columns past column K.  The 'rows' array
   is M long, and the 'cols' array is N long.  The elements in both arrays
   after the first K contain the indexes of the rows and columns not selected
   to be part of the sub-matrix of A, in arbitrary order.

   If A is not of rank K or more, L will contain fewer than K non-zero columns,
   and U will contain an equal number of non-zero rows.  The entries in the
   'rows' and 'cols' arrays for the extra zero rows or columns will be
   arbitrary (but valid).  The number of extra zero columns is returned as
   the value of this procedure (hence a return value of zero indicates that
   a K by K sub-matrix of full rank was found).

   The matrix A is not altered.  The previous contents of L and U are cleared.
*/

int LDPC_mod2::mod2sparse_decomp
( mod2sparse *A,	/* Input matrix, M by N */
  int K,		/* Size of sub-matrix to find LU decomposition of */
  mod2sparse *L,	/* Matrix in which L is stored, M by K */
  mod2sparse *U,	/* Matrix in which U is stored, K by N */
  int *rows,		/* Array where row indexes are stored, M long */
  int *cols,		/* Array where column indexes are stored, N long */
  mod2sparse_strategy strategy, /* Strategy to follow in picking rows/columns */
  int abandon_number,	/* Number of columns to abandon at some point */
  int abandon_when	/* When to abandon these columns */
)
{
  int *rinv, *cinv, *acnt, *rcnt;
  mod2sparse *B;
  int M, N;

  mod2entry *e, *f, *fn, *e2, *e3;
  int i, j, k, cc, cc2, cc3, cr2, mxr, pr, fnd;
  int found, nnf;

  M = mod2sparse_rows(A);
  N = mod2sparse_cols(A);

  if (mod2sparse_cols(L)!=K || mod2sparse_rows(L)!=M
   || mod2sparse_cols(U)!=N || mod2sparse_rows(U)!=K)
  { fprintf (stderr,
      "mod2sparse_decomp: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (abandon_number>N-K)
  { fprintf(stderr,"Trying to abandon more columns than allowed\n");
    exit(1);
  }

  rinv = (int *)calloc (M, sizeof *rinv);
  cinv = (int *)calloc (N, sizeof *cinv);

  if (abandon_number>0)
  { acnt = (int *)calloc (M+1, sizeof *acnt);
  }

  if (strategy==Mod2sparse_minprod)
  { rcnt = (int *)calloc (M, sizeof *rcnt);
  }

  mod2sparse_clear(L);
  mod2sparse_clear(U);

  /* Copy A to B.  B will be modified, then discarded. */

  B = mod2sparse_allocate(M,N);
  mod2sparse_copy(A,B);

  /* Count 1s in rows of B, if using minprod strategy. */

  if (strategy==Mod2sparse_minprod)
  { for (i = 0; i<M; i++)
    { rcnt[i] = mod2sparse_count_row(B,i);
    }
  }

  /* Set up initial row and column choices. */

  for (i = 0; i<M; i++) rows[i] = rinv[i] = i;
  for (j = 0; j<N; j++) cols[j] = cinv[j] = j;

  /* Find L and U one column at a time. */

  nnf = 0;

  for (i = 0; i<K; i++)
  {
    /* Choose the next row and column of B. */

    switch (strategy)
    {
      case Mod2sparse_first:
      {
        found = 0;

        for (k = i; k<N; k++)
        { e = mod2sparse_first_in_col(B,cols[k]);
          while (!mod2sparse_at_end(e))
          { if (rinv[mod2sparse_row(e)]>=i)
            { found = 1;
              goto out_first;
            }
            e = mod2sparse_next_in_col(e);
          }
        }

      out_first:
        break;
      }

      case Mod2sparse_mincol:
      {
        found = 0;

        for (j = i; j<N; j++)
        { cc2 = mod2sparse_count_col(B,cols[j]);
          if (!found || cc2<cc)
          { e2 = mod2sparse_first_in_col(B,cols[j]);
            while (!mod2sparse_at_end(e2))
            { if (rinv[mod2sparse_row(e2)]>=i)
              { found = 1;
                cc = cc2;
                e = e2;
                k = j;
                break;
              }
              e2 = mod2sparse_next_in_col(e2);
            }
          }
        }

        break;
      }

      case Mod2sparse_minprod:
      {
        found = 0;

        for (j = i; j<N; j++)
        { cc2 = mod2sparse_count_col(B,cols[j]);
          e2 = mod2sparse_first_in_col(B,cols[j]);
          while (!mod2sparse_at_end(e2))
          { if (rinv[mod2sparse_row(e2)]>=i)
            { cr2 = rcnt[mod2sparse_row(e2)];
              if (!found || cc2==1 || (cc2-1)*(cr2-1)<pr)
              { found = 1;
                pr = cc2==1 ? 0 : (cc2-1)*(cr2-1);
                e = e2;
                k = j;
              }
            }
            e2 = mod2sparse_next_in_col(e2);
          }
        }

        break;
      }

      default:
      { fprintf(stderr,"mod2sparse_decomp: Unknown stategy\n");
        exit(1);
      }
    }

    if (!found)
    { nnf += 1;
    }

    /* Update 'rows' and 'cols'.  Looks at 'k' and 'e' found above. */

    if (found)
    {
      if (cinv[mod2sparse_col(e)]!=k) abort();

      cols[k] = cols[i];
      cols[i] = mod2sparse_col(e);

      cinv[cols[k]] = k;
      cinv[cols[i]] = i;

      k = rinv[mod2sparse_row(e)];

      if (k<i) abort();

      rows[k] = rows[i];
      rows[i] = mod2sparse_row(e);

      rinv[rows[k]] = k;
      rinv[rows[i]] = i;
    }

    /* Update L, U, and B. */

    f = mod2sparse_first_in_col(B,cols[i]);

    while (!mod2sparse_at_end(f))
    {
      fn = mod2sparse_next_in_col(f);
      k = mod2sparse_row(f);

      if (rinv[k]>i)
      { mod2sparse_add_row(B,k,B,mod2sparse_row(e));
        if (strategy==Mod2sparse_minprod)
        { rcnt[k] = mod2sparse_count_row(B,k);
        }
        mod2sparse_insert(L,k,i);
      }
      else if (rinv[k]<i)
      { mod2sparse_insert(U,rinv[k],cols[i]);
      }
      else
      { mod2sparse_insert(L,k,i);
        mod2sparse_insert(U,i,cols[i]);
      }

      f = fn;
    }

    /* Get rid of all entries in the current column of B, just to save space. */

    for (;;)
    { f = mod2sparse_first_in_col(B,cols[i]);
      if (mod2sparse_at_end(f)) break;
      mod2sparse_delete(B,f);
    }

    /* Abandon columns of B with lots of entries if it's time for that. */

    if (abandon_number>0 && i==abandon_when)
    {
      for (k = 0; k<M+1; k++)
      { acnt[k] = 0;
      }
      for (j = 0; j<N; j++)
      { k = mod2sparse_count_col(B,j);
        acnt[k] += 1;
      }

      cc = abandon_number;
      k = M;
      while (acnt[k]<cc)
      { cc -= acnt[k];
        k -= 1;
        if (k<0) abort();
      }

      cc2 = 0;
      for (j = 0; j<N; j++)
      { cc3 = mod2sparse_count_col(B,j);
        if (cc3>k || cc3==k && cc>0)
        { if (cc3==k) cc -= 1;
          for (;;)
          { f = mod2sparse_first_in_col(B,j);
            if (mod2sparse_at_end(f)) break;
            mod2sparse_delete(B,f);
          }
          cc2 += 1;
        }
      }

      if (cc2!=abandon_number) abort();

      if (strategy==Mod2sparse_minprod)
      { for (j = 0; j<M; j++)
        { rcnt[j] = mod2sparse_count_row(B,j);
        }
      }
    }
  }

  /* Get rid of all entries in the rows of L past row K, after reordering. */

  for (i = K; i<M; i++)
  { for (;;)
    { f = mod2sparse_first_in_row(L,rows[i]);
      if (mod2sparse_at_end(f)) break;
      mod2sparse_delete(L,f);
    }
  }

  mod2sparse_free(B);
  free(rinv);
  free(cinv);
  if (strategy==Mod2sparse_minprod) free(rcnt);
  if (abandon_number>0) free(acnt);

  return nnf;
}


/* SOLVE A LOWER-TRIANGULAR SYSTEM BY FORWARD SUBSTITUTION.  Solves Ly=x
   for y by forward substitution, based on L being lower triangular after
   its rows are reordered according to the given index array.  The vectors
   x and y are stored unpacked, one bit per character.  If L is M by K,
   then x should be M long, but only the K bits indexed by 'rows' are looked
   at.  The solution vector, y, must be K long.  Only K rows of L are used,
   as also determined by the K indexes in the 'rows' argument.  If 'rows' is
   null, the first K rows of L and the first K elements of x are used.

   If the matrix L does not have 1s on its diagonal (after row rearrangement),
   there may be no solution, depending on what x is.  If no solution exists,
   this procedure returns zero, otherwise it returns one.  Any arbitrary
   bits in the solution are set to zero.

   If L is not lower-triangular, a message is displayed on standard error and
   the program is terminated.
*/

int LDPC_mod2::mod2sparse_forward_sub
( mod2sparse *L,	/* Matrix that is lower triangular after reordering */
  int *rows,		/* Array of indexes (from 0) of rows for new order */
  char *x,		/* Vector on right of equation, also reordered */
  char *y		/* Place to store solution */
)
{
  int K, i, j, ii, b, d;
  mod2entry *e;

  K = mod2sparse_cols(L);

  /* Make sure that L is lower-triangular, after row re-ordering. */

  for (i = 0; i<K; i++)
  { ii = rows ? rows[i] : i;
    e = mod2sparse_last_in_row(L,ii);
    if (!mod2sparse_at_end(e) && mod2sparse_col(e)>i)
    { fprintf(stderr,
        "mod2sparse_forward_sub: Matrix is not lower-triangular\n");
      exit(1);
    }
  }

  /* Solve system by forward substitution. */

  for (i = 0; i<K; i++)
  {
    ii = rows ? rows[i] : i;

    /* Look at bits in this row, forming inner product with partial
       solution, and seeing if the diagonal is 1. */

    d = 0;
    b = 0;

    for (e = mod2sparse_first_in_row(L,ii);
         !mod2sparse_at_end(e);
         e = mod2sparse_next_in_row(e))
    {
      j = mod2sparse_col(e);

      if (j==i)
      { d = 1;
      }
      else
      { b ^= y[j];
      }
    }

    /* Check for no solution if the diagonal isn't 1. */

    if (!d && b!=x[ii])
    { return 0;
    }

    /* Set bit of solution, zero if arbitrary. */

    y[i] = b^x[ii];
  }

  return 1;
}


/* SOLVE AN UPPER-TRIANGULAR SYSTEM BY BACKWARD SUBSTITUTION.  Solves Uz=y
   for z by backward substitution, based on U being upper triangular after
   its colums are reordered according to the given index array.  The vectors
   y and z are stored unpacked, one bit per character.  If U is K by N,
   then the solution vector, z, should be N long, but only the K bits indexed
   by 'cols' are set.  The vector y must be K long.  Only K columns of U are
   used, as also determined by the K indexes in the 'cols' argument.  The other
   columns of U must be zero (this is not checked, but is necessary for the
   method used to work).  If 'cols' is null, the first K columns of U and the
   first K elements of z are used.

   If the matrix U does not have 1s on its diagonal (after column rearrangement)
   there may be no solution, depending on what y is.  If no solution exists,
   this procedure returns zero, otherwise it returns one.  Any arbitrary
   bits in the solution are set to zero.

   If U is not upper-triangular, a message is displayed on standard error and
   the program is terminated.
*/

int LDPC_mod2::mod2sparse_backward_sub
( mod2sparse *U,	/* Matrix that is upper triangular after reordering */
  int *cols,		/* Array of indexes (from 0) of columns for new order */
  char *y,		/* Vector on right of equation */
  char *z		/* Place to store solution, also reordered */
)
{
  int K, i, j, ii, b, d;
  mod2entry *e;

  K = mod2sparse_rows(U);

  /* Make sure that U is upper-triangular, after column re-ordering. */

  for (i = 0; i<K; i++)
  { ii = cols ? cols[i] : i;
    e = mod2sparse_last_in_col(U,ii);
    if (!mod2sparse_at_end(e) && mod2sparse_row(e)>i)
    { fprintf(stderr,
        "mod2sparse_backward_sub: Matrix is not upper-triangular\n");
      exit(1);
    }
  }

  /* Solve system by backward substitution. */

  for (i = K-1; i>=0; i--)
  {
    ii = cols ? cols[i] : i;

    /* Look at bits in this row, forming inner product with partial
       solution, and seeing if the diagonal is 1. */

    d = 0;
    b = 0;

    for (e = mod2sparse_first_in_row(U,i);
         !mod2sparse_at_end(e);
         e = mod2sparse_next_in_row(e))
    {
      j = mod2sparse_col(e);

      if (j==ii)
      { d = 1;
      }
      else
      { b ^= z[j];
      }
    }

    /* Check for no solution if the diagonal isn't 1. */

    if (!d && b!=y[i])
    { return 0;
    }

    /* Set bit of solution, zero if arbitrary. */

    z[ii] = b^y[i];
  }

  return 1;
}


/* READ AN INTEGER ONE BYTE AT A TIME.  Four bytes are read, ordered from
   low to high order.  These are considered to represent a signed integer,
   in two's complement form.  The value returned is this integer, converted
   to whatever a C "int" is.  The conversion should work as long as an "int"
   is at least four bytes, even if it's not in two's complement representation
   (except for the largest two's complement negative integer).

   If an error or eof is encountered, zero is returned.  The caller can
   check for these events using feof and ferror.

   The file read from should have been opened as "binary".
*/

int LDPC_mod2::intio_read
( FILE *f   /* File to read from */
)
{
  unsigned char b[4];
  int top;
  int i;

  for (i = 0; i<4; i++)
  {
	  if (fread(&b[i],1,1,f) != 1)
		return 0;
  }

  top = b[3]>127 ? (int)b[3] - 256 : b[3];

  return (top<<24) + (b[2]<<16) + (b[1]<<8) + b[0];
}


/* WRITE AN INTEGER ONE BYTE AT A TIME.  Four bytes are written, ordered from
   low to high order.  These are considered to represent a signed integer,
   in two's complement form.  This should work as long as the integer passed
   can be represented in four bytes, even if a C "int" is longer than this.

   The file written to should have been opened as "binary".
*/

void LDPC_mod2::intio_write
( FILE *f,  /* File to write to */
  int v     /* Value to write to file */
)
{
  unsigned char b;
  int i;

  for (i = 0; i<3; i++)
  { b = v&0xff;
    fwrite(&b,1,1,f);
    v >>= 8;
  }

  b = v>0 ? v : v+256;
  fwrite(&b,1,1,f);
}



void *LDPC_mod2::mod2sparse_remove
( mod2sparse *m,
  int row,
  int col
)
{ // added by dgq on 20100907am: remove one entry.

  mod2entry *re, *ce, *ne;

  re = mod2sparse_last_in_row(m,row); //the last element in row-th row.

  if (!mod2sparse_at_end(re) && mod2sparse_col(re)==col) //if ((e)->row>0) => this is valid entry. //!mod2sparse_at_end(re) == re is valid entry.
  {//if this entry already exist, then return;
	  mod2sparse_delete(m,re);
	  return 0; //
  }


  //find the entry on the right of insert entry:
  if (mod2sparse_at_end(re) || mod2sparse_col(re)<col) //if ((e)->row<0), this entry is not valid entry, and it should be the H->rows[i], which is not entry, but a indicator/pointer.
  {//if there is no element in current chain, then insert it.
	  //re = re->right; //the result is the root entry:
	  printf("No such entry.\n");
	  exit(0);
	  return 0;
  }
  else
  {//if there are already entries in the chain, then find the position to insert the entry:
    re = mod2sparse_first_in_row(m,row);
    for (;;)
    {
      if (!mod2sparse_at_end(re) && mod2sparse_col(re)==col)
      {
		  mod2sparse_delete(m,re);
		  return 0;
      }
      if (mod2sparse_at_end(re) || mod2sparse_col(re)>col) //find the first element with col larger than the input col. The new element should be added before this found element.
      {
		  printf("No such entry.\n");
		  exit(0);
		  return 0;
		  //break;
      }

	  //if this element's col value smaller than expected, then go to check the next element to see whether its col larger than expected.
	  //this shows that ->right ->right is the order of col increasing, as common sense.
      re = mod2sparse_next_in_row(re);
    }
  }

  //if reach here, it means there is not entry at such position.
	//
  printf("No such entry.\n");
  exit(0);

}
