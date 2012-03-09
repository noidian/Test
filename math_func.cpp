#include "math_func.h"


long double norm(long double *val, int numberofvalues)
{
    long double net_val;
    register int i;

    net_val=0;
    for (i=0;i<numberofvalues;i++)
        net_val+=pow(val[i],2);

    net_val=sqrt(net_val);

return net_val;
}


long double norm(int *val, int numberofvalues)
{
    long double net_val;
    register int i;

    net_val=0;
    for (i=0;i<numberofvalues;i++)
        net_val+=pow(val[i],2);

    net_val=sqrt(net_val);

return net_val;
}


long double norm(int **val, int rows, int cols)
{
    long double net_val;
    register int i,j;

    net_val=0;
    for (i=0;i<rows;i++)
        for (j=0;j<cols;j++)
            net_val+=pow(val[i][j],2);

    net_val=sqrt(net_val);

return net_val;
}

long double norm(long double **val, int rows, int cols)
{
    long double net_val;
    register int i,j;

    net_val=0;
    for (i=0;i<rows;i++)
        for (j=0;j<cols;j++)
            net_val+=pow(val[i][j],2);

    net_val=sqrt(net_val);

return net_val;
}



int array_divide (long double *val,long double div_val, int numberofvalues)
{
   register int i;

   for (i=0 ;i<numberofvalues;i++)
        val[i]=val[i]/div_val;

return OK;
}


int array_addition (long double *output,long double *input1,long double *input2, int numberofvalues)
{
   register int i;
   for (i=0 ;i<numberofvalues;i++)
        output[i]=input1[i]+input2[i];

return OK;
}


int array_divide (long double **val,long double div_val, int rows,int numberofvalues)
{
   register int i,j;

    for (i=0;i<rows;i++)
        for (j=0 ;j<numberofvalues;j++)
            val[i][j]=val[i][j]/div_val;

return OK;
}


int array_multiply (long double *val,long double mult_val, int numberofvalues)
{
   register int i;

   for (i=0 ;i<numberofvalues;i++)
        val[i]=val[i]*mult_val;

return OK;
}



int  x_corr(int *x,int *y,int cshift, int numberofvals, long double *x_corr_op)
{
    register int i,j;

    for (i=0;i<(2*(numberofvals-1));i++)
    {

        if (i<=numberofvals-1 && i>=(numberofvals-1-cshift)) //e.g. for 10000, 9987 is turning point
        {
            for (j=0;j<=i;j++)
                x_corr_op[i]=x_corr_op[i]+(x[(numberofvals-1)-i+j]*y[j]);
        }
        else if ((i>numberofvals-1)&&(i<=(numberofvals+cshift-1)))
        {
            for (j=0;j<=(2*(numberofvals-1))-i;j++)
                x_corr_op[i]=x_corr_op[i]+(y[i-(numberofvals-1)+j]*x[j]);
        }
    }


    return OK;
}


int  x_corr(int *x,long double *y,int cshift, int numberofvals, long double *x_corr_op)
{
    register int i,j;

    for (i=0;i<(2*(numberofvals-1));i++)
    {

        if (i<=numberofvals-1 && i>=(numberofvals-1-cshift)) //e.g. for 10000, 9987 is turning point
        {
            for (j=0;j<=i;j++)
                x_corr_op[i]=x_corr_op[i]+(x[(numberofvals-1)-i+j]*y[j]);
        }
        else if ((i>numberofvals-1)&&(i<=(numberofvals+cshift-1)))
        {
            for (j=0;j<=(2*(numberofvals-1))-i;j++)
                x_corr_op[i]=x_corr_op[i]+(y[i-(numberofvals-1)+j]*x[j]);
        }
    }


    return OK;
}


long double  min(long double input1,long double input2)
{
    long double output;

    if (input1>input2)
        output=input2;
    else
        output=input1;

    return output;
}



int  x_corr(long double *x,long double *y,int cshift, int numberofvals, long double *x_corr_op)
{
    register int i,j;

    for (i=0;i<(2*(numberofvals-1));i++)
    {

        if (i<=numberofvals-1 && i>=(numberofvals-1-cshift)) //e.g. for 10000, 9987 is turning point
        {
            for (j=0;j<=i;j++)
                x_corr_op[i]=x_corr_op[i]+(x[(numberofvals-1)-i+j]*y[j]);
        }
        else if ((i>numberofvals-1)&&(i<=(numberofvals+cshift-1)))
        {
            for (j=0;j<=(2*(numberofvals-1))-i;j++)
                x_corr_op[i]=x_corr_op[i]+(y[i-(numberofvals-1)+j]*x[j]);
        }
    }
    return OK;
}

int x_corr2(long double **input1,int rows1, int cols1,long double **input2, int rows2, int cols2, int rshift,int cshift,int Nb_samples, long double **output)
{

    int x_corr_cols,x_corr_rows;

    int row_shift,column_shift,row_shift_val,column_shift_val;

    register int i,j;

    x_corr_cols=cols1+cols2-1;
    x_corr_rows=rows1+rows2-1;

    int temp;
    //x_corr_op=output
    //Keep 1st matrix static and other one to convolve
    //Notchecked for full width works for required values!!
  for (row_shift=-(rshift-1);row_shift<=(rshift-1);row_shift++)
  {
    row_shift_val=abs(row_shift);
    for (column_shift=-cshift;column_shift<=cshift;column_shift++)
    {
        column_shift_val=abs(column_shift);
        for (i=row_shift_val;i<=rows1;i++)
        {
            for (j=column_shift_val;j<=cols1;j++)
            {
                if (j==cols1)
                    temp=0;
                if ((row_shift>=0) && (column_shift>=0))
                {
                    output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i-row_shift][j-column_shift]*input2[i][j];
                }
                if ((row_shift < 0) && (column_shift>=0))
                {
                     output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i][j-column_shift]*input2[i+row_shift][j];
                }
                if ((row_shift >=0) && (column_shift<0))
                {
                     output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i-row_shift][j]*input2[i][j+column_shift];
                }
                if ((row_shift <0) && (column_shift<0))
                {
                     output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i][j]*input2[i+row_shift][j+column_shift];
                }
            }
        }
    }
  }

    return OK;
}

int x_corr2(int **input1,int rows1, int cols1,long double **input2, int rows2, int cols2, int rshift,int cshift,int Nb_samples, long double **output)
{

    int x_corr_cols,x_corr_rows;

    int row_shift,column_shift,row_shift_val,column_shift_val;

    register int i,j;

    x_corr_cols=cols1+cols2-1;
    x_corr_rows=rows1+rows2-1;

    //x_corr_op=output
    //Keep 1st matrix static and other one to convolve
    //Notchecked for full width works for required values!!
  for (row_shift=-(rshift-1);row_shift<=(rshift-1);row_shift++)
  {
    row_shift_val=abs(row_shift);
    for (column_shift=-cshift;column_shift<=cshift;column_shift++)
    {
        column_shift_val=abs(column_shift);
        for (i=row_shift_val;i<=rows1;i++)
        {
            for (j=column_shift_val;j<=cols1;j++)
            {

                if ((row_shift>=0) && (column_shift>=0))
                {
                    output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i-row_shift][j-column_shift]*input2[i][j];
                }
                if ((row_shift < 0) && (column_shift>=0))
                {
                     output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i][j-column_shift]*input2[i+row_shift][j];
                }
                if ((row_shift >=0) && (column_shift<0))
                {
                     output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i-row_shift][j]*input2[i][j+column_shift];
                }
                if ((row_shift <0) && (column_shift<0))
                {
                     output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i][j]*input2[i+row_shift][j+column_shift];
                }
            }
        }
    }
  }

    return OK;
}


int x_corr2(int **input1,int rows1, int cols1,int **input2, int rows2, int cols2, int rshift,int cshift,int Nb_samples, long double **output)
{

    int x_corr_cols,x_corr_rows;

    int row_shift,column_shift,row_shift_val,column_shift_val;

    register int i,j;

    x_corr_cols=cols1+cols2-1;
    x_corr_rows=rows1+rows2-1;
    //x_corr_op=output
    //Keep 1st matrix static and other one to convolve
    //Notchecked for full width works for required values!!
  for (row_shift=-(rshift-1);row_shift<=(rshift-1);row_shift++)
  {
    row_shift_val=abs(row_shift);
    for (column_shift=-cshift;column_shift<=cshift;column_shift++)
    {
        column_shift_val=abs(column_shift);
        for (i=row_shift_val;i<=rows1;i++)
        {
            for (j=column_shift_val;j<=cols1;j++)
            {
                if ((row_shift>=0) && (column_shift>=0))
                {
                    output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i-row_shift][j-column_shift]*input2[i][j];
                }
                if ((row_shift < 0) && (column_shift>=0))
                {
                     output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i][j-column_shift]*input2[i+row_shift][j];
                }
                if ((row_shift >=0) && (column_shift<0))
                {
                     output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i-row_shift][j]*input2[i][j+column_shift];
                }
                if ((row_shift <0) && (column_shift<0))
                {
                     output[rows1+row_shift][cshift+column_shift]=output[rows1+row_shift][cshift+column_shift]+input1[i][j]*input2[i+row_shift][j+column_shift];
                }
            }
        }
    }
  }

    return OK;
}


int matrix_transpose(long double **input_transpose,int rows,int cols, long double **input)

{
    register int i,j;

       for (i=0;i<cols;i++)
    {
        for (j=0;j<rows;j++)
        {
            input_transpose[i][j]=input[j][i];
        }
    }


    return OK;

}

//MxN with NxQ matrix
int matrix_multiply (long double **output,long double **input1,int row1, int col1,long double **input2,int row2,int col2)
{
    register int i,j,k;

       for (i=0;i<row1;i++)
    {
        for (j=0;j<col2;j++)
        {

            for (k=0;k<col1;k++)
            {
                output[i][j]=output[i][j]+(input1[i][k] *input2[k][j]);
            }
        }
    }
    return OK;
}



//1xN with NxQ matrix
int matrix_multiply (long double *output,long double *input1,int col1,long double **input2,int row2,int col2)
{
    register int i,j;

       for (i=0;i<col2;i++)
    {
        for (j=0;j<col1;j++)
        {
           output[i]=output[i]+(input1[j]*input2[j][i]);
        }
    }
    return OK;
}

//NxQ matrix with Qx1
int matrix_multiply (long double *output,long double **input1,int row1,int col1,long double *input2,int row2)
{
    register int i,j;

       for (i=0;i<row1;i++)
    {
        for (j=0;j<col1;j++)
        {
           output[i]=output[i]+(input1[i][j]*input2[j]);
        }
    }
    return OK;
}





int  matrix_subtraction(long double **output,long double **input1,long double **input2,int row,int col)
{
   register int i,j;

       for (i=0;i<row;i++)
    {
        for (j=0;j<col;j++)
        {
                output[i][j]=(input1[i][j] - input2[i][j]);
        }
    }
    return OK;
}

int matrix_inverse(long double **output,int row,int col,long double **input)
{
    register int i,j;

    gsl_matrix *m;
    m=gsl_matrix_calloc (row,col);

    for (i=0;i<row;i++)
    {
        for (j=0;j<col;j++)
        {
           gsl_matrix_set(m,i,j,input[i][j]);
        }
    }

    gsl_matrix *inv;
    gsl_permutation *p;
    int *signum;

    inv=gsl_matrix_calloc (row,col);
    p=gsl_permutation_calloc(row);
    signum=(int *)malloc( sizeof(int)*row);

    gsl_linalg_LU_decomp(m,p,signum);
    gsl_linalg_LU_invert(m,p,inv);


    for (i=0;i<row;i++)
    {
        for (j=0;j<col;j++)
        {
            output[i][j]=(long double) gsl_matrix_get(inv,i,j);;
        }
    }

    free(signum);
    gsl_matrix_free (m);
    gsl_matrix_free (inv);
    gsl_permutation_free(p);

    return OK;
}
