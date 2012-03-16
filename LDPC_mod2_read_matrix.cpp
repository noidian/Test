#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "LDPC_mod2.h"







mod2sparse * LDPC_mod2::read_binary_matrix_H_to_mod2sparse(char *filename, int n_col, int m_row)
{
//read in the 0,1 matrix file, and then turn it into mod2sparse format:
	int i,j;

	mod2sparse *H; //* Representation of a sparse matrix */
	mod2entry *e;
	int N1,M1;
	FILE *fp;

	N1=n_col;
	M1=m_row;
	N=n_col;
	M=m_row;

	//1: allocate memory space first:
	H = (mod2sparse *)calloc(1, sizeof(*H));
	H->n_rows = M1;
	H->n_cols = N1;

	H->rows = (mod2entry *) calloc(H->n_rows,sizeof(mod2entry)); //the first element in each row.
	H->cols = (mod2entry *) calloc(H->n_cols,sizeof(mod2entry)); //the first element in each column.

	H->blocks = 0;
	H->next_free = 0;

	//initial the indicator in all rows and all columns:
	//only N+M indicator exist at the beginning.
	//NOTE: this is actually not a entry, just a indicator, and its row and col value are always -1.
	for(i=0; i<M1; i++)
	{
		e= & H->rows[i];
		e->left = e->right = e->up = e-> down = e;
		e->row = e->col =-1;
	}
	for(i=0; i<N1; i++)
	{
		e= & H->cols[i];
		e->left = e->right = e->up = e-> down = e;
		e->row = e->col =-1;
	}


	//---read file and insert entries:
	int it;
	fp=fopen(filename,"r");
	if(fp==0)
	{
		printf("filename is wrong. exit.\n");
		exit(0);
	}
	for(j=0;j<M1;j++)
	{
		for(i=0;i<N1;i++)
		{
			fscanf(fp,"%1d",&it);
			if(it==1)
			{//insert [j,i] entry to H:
				mod2sparse_insert(H,j,i);
			}
		}//for i.
	}//for j.
	fclose(fp);


return H;

}

mod2sparse * LDPC_mod2::read_H_positions_in_rows_to_mod2sparse(char *filename, int n_col, int m_row, int row_weight)
{
//read in the 0,1 matrix file, and then turn it into mod2sparse format:
	int i,j;

	mod2sparse *H; //* Representation of a sparse matrix */
	mod2entry *e;
	int N1,M1;
	FILE *fp;

	N1=n_col;
	M1=m_row;
	N=n_col;
	M=m_row;

	//1: allocate memory space first:
	H = (mod2sparse *)calloc(1, sizeof(*H));
	H->n_rows = M1;
	H->n_cols = N1;


	H->rows = (mod2entry *) calloc(H->n_rows,sizeof(mod2entry)); //the first element in each row.
	H->cols = (mod2entry *) calloc(H->n_cols,sizeof(mod2entry)); //the first element in each column.

	H->blocks = 0;
	H->next_free = 0;

	//initial the indicator in all rows and all columns:
	//only N+M indicator exist at the beginning.
	//NOTE: this is actually not a entry, just a indicator, and its row and col value are always -1.
	for(i=0; i<M1; i++)
	{
		e= & H->rows[i];
		e->left = e->right = e->up = e-> down = e;
		e->row = e->col =-1;
	}
	for(i=0; i<N1; i++)
	{
		e= & H->cols[i];
		e->left = e->right = e->up = e-> down = e;
		e->row = e->col =-1;
	}


	//---read file and insert entries:
	int it;
	fp=fopen(filename,"r");
	if(fp==0)
	{
		printf("filename is wrong. exit.\n");
		exit(0);
	}
	for(j=0;j<M1;j++)
	{
		for(i=0;i<row_weight;i++)
		{
			fscanf(fp,"%d",&it);
			{//insert [j,i] entry to H:
				mod2sparse_insert(H,j,it);
			}
		}//for i.
	}//for j.
	fclose(fp);


return H;

}


mod2sparse * LDPC_mod2::read_RS_LDPC_base_matrix_to_sparse_matrix(
				char * input_filename,
				int b,
				char * output_filename,
				int row,
				int col
				)
{
// input: the large basic matrix, which is b x b size, with each element in GF(b).
	mod2sparse *H; //* Representation of a sparse matrix */
	mod2entry *e;
	int N1,M1;
	int i,j;
	FILE *fp;

	int it;
	int base_row1,base_col1,k;
	int **base_h;
	int **selected_base_h;

	N1=col*b;
	M1=row*b;

	//---------------------------------------------------------------------------------------------
	//1: allocate memory space first:
	H = (mod2sparse *)calloc(1, sizeof(*H));
	H->n_rows = M1;
	H->n_cols = N1;

	H->rows = (mod2entry *) calloc(H->n_rows,sizeof(mod2entry)); //the first element in each row.
	H->cols = (mod2entry *) calloc(H->n_cols,sizeof(mod2entry)); //the first element in each column.

	H->blocks = 0;
	H->next_free = 0;


	//---------------------------------------------------------------------------------------------
	//2: initial the indicator in all rows and all columns:
	//only N+M indicator exist at the beginning.
	//NOTE: this is actually not a entry, just a indicator, and its row and col value are always -1.
	for(i=0; i<M1; i++)
	{
		e= & H->rows[i];
		e->left = e->right = e->up = e-> down = e;
		e->row = e->col =-1;
	}
	for(i=0; i<N1; i++)
	{
		e= & H->cols[i];
		e->left = e->right = e->up = e-> down = e;
		e->row = e->col =-1;
	}


	//---------------------------------------------------------------------------------------------
	//3: read QC_BASE_MATRIX:
	base_h = (int **)calloc(b, sizeof(int *));
	for(i=0;i<b;i++)
		base_h[i] = (int *)calloc(b, sizeof(int ));

	selected_base_h=(int **)calloc(row, sizeof(int *));
	for(i=0;i<row;i++)
		selected_base_h[i] = (int *)calloc(col, sizeof(int ));
	//
	fp=fopen(input_filename,"r");
	for(j=0;j<b;j++)
	{
		for(i=0;i<b;i++)
		{
			fscanf(fp," %d",&it);
			base_h[j][i]=it;
		}//for i.
	}//for j.
	fclose(fp);

	//get the selected_non-binary-H:
	base_row1 = col;
	base_col1 = 0;
	for(j=0;j<row;j++)
	for(i=0;i<col;i++)
	{
		selected_base_h[j][i] = base_h[j+base_row1][i+base_col1];
	}//for i.

//	show_int_matrix(selected_base_h,row,col,"./selected_base_h.txt");

	for(j=0;j<row;j++)
	{
		for(k=0;k<b;k++)
		{//each row in j-th circulant-row:
			for(i=0;i<col;i++)
			{//each COC in the j-th circulant-row:
				it=selected_base_h[j][i];
				//[it+k]
				mod2sparse_insert(H,  (j*b+k),  (i*b+ (it-1+k+b)%b));
			}//for i.
		}//for k.
	}//for j.

//	show_row_entry_matrix(matrix_H,row*b,col*b,output_filename);

	//free:
	for(i=0;i<b;i++)
		free(base_h[i]);
	free(base_h);

	for(i=0;i<row;i++)
		free(selected_base_h[i]);
	free(selected_base_h);

return H;
}


