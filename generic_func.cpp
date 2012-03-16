

#include "generic_func.h"


int fileopen(char *filename,int *input,int numberofvalues)
{
	register int i;

	FILE *fp;
	fp=fopen(filename,"r");
	for(i=0;i<numberofvalues;i++)
		fscanf(fp,"%d",&input[i]);
	fclose(fp);

	return OK;
}


int fileopen(char *filename,char *input,int numberofvalues)
{
	register int i;

	FILE *fp;
	fp=fopen(filename,"r");
	for(i=0;i<numberofvalues;i++)
	{
	    fscanf(fp,"%c",&input[i]);
	    if (input[i]=='\n')             //accounts for the input data stored in next line
            i--;
        if (input[i]=='-')             //remap -1 to 0
        {
           fscanf(fp,"%c",&input[i]);
           input[i]='\000';
        }
        if (input[i]=='1')
        {
            input[i]='\001';
        }
	}

	fclose(fp);

	return OK;
}


int fileopen(char *filename,long double *input,int numberofvalues)
{
	register int i;

	FILE *fp;
	fp=fopen(filename,"r");
	for(i=0;i<numberofvalues;i++)
		fscanf(fp,"%Lf",&input[i]);
	fclose(fp);

	return OK;
}

int array_append(int *output,int *input1,int *input2, int append_length,int end_val)
{
    register int i;
    for (i=0;i<end_val;i++)
     {
        if(i<append_length)
            output[i]=input1[i];
        else if(i<end_val-append_length)
            output[i]=input2[i-append_length];
        else
            output[i]=input1[i-(end_val-append_length)];
     }

return OK;
}


int array_resize(int *input,int start_val,int end_val)
{
    register int i;
    for (i=start_val;i<end_val;i++)
        input[i-start_val]=input[i];
    input[i-start_val]=NULL;
return OK;
}

int array_resize(long double *input,int start_val,int end_val)
{
    register int i;
    for (i=start_val;i<end_val;i++)
        input[i-start_val]=input[i];
    input[i-start_val]=NULL;
return OK;
}


int shiftback_wd (char *cwd,int length,char *opwd,int shiftback_number)
{
    //j must be > shiftback_number
    register int i,j;
    int slash_location[10],stop_val;
    j=0;
    for (i=0;i<length;i++)
    {
        if (cwd[i]=='/')
        {  slash_location[j]=i;
           j++;
        }
    }

    stop_val=slash_location[(j)-shiftback_number];

    for (i=0;i<stop_val;i++)
        opwd[i]=cwd[i];
    opwd[i]='\0';

return OK;
}

int dot2underscore(char *str, int length)
{
    register int i;
    for (i=0;i<length;i++)
    {
        if (str[i]=='.')
        {
            str[i]='_';
        }
    }
return OK;
}
int dot2null(char *str, int length)
{
    register int i;
    for (i=0;i<length;i++)
    {
        if (str[i]=='.')
        {
            str[i]='\0';
        }
    }
return OK;
}


//calloc issues!
long double **matrix_calloc(int row, int col)
{
    register int i;
    long double **input;


    input=(long double **)malloc(sizeof(long double* ) * row);
    for( i = 0; i < row; i++ )
        input[i]=(long double *)malloc( sizeof(long double ) * col );

    var_init(input,row,col);
    return (input);
}


long double ***matrix_calloc(int x, int y,int z)
{
    register int i,j;
    long double ***input;


    input=(long double ***)malloc(sizeof(long double** ) * x);
    for (i= 0;i < x; i++)
    {
        input[i]=(long double **)malloc( sizeof(long double ) * y );
        for( j = 0; j < y; j++ )
        {
            input[i][j]=(long double *)malloc( sizeof(long double ) * z );
        }
    }

    var_init(input,x,y,z);
    return (input);
}



int **matrix_calloc_int(int row, int col)
{
    register int i;
    int **input;

    input=(int **)malloc(sizeof(int* ) * row);
    for( i = 0; i < row; i++ )
        input[i]=(int *)malloc( sizeof(int) * col );

    var_init(input,row,col);
    return (input);
}

int array_copy(long double *val1,long double *val, int start_val, int end_val)
{
   register int i,j;
    j=0;
   for (i=start_val ;i<end_val;i++)
        {
            val1[j]=val[i];
            j++;
        }

return OK;
}


int array_copy (long double *val1,long double *val, int numberofvalues)
{
   register int i;

   for (i=0 ;i<numberofvalues;i++)
        val1[i]=val[i];

return OK;
}


int twoD_array_copy_int (int **val1,int *val, int row,int col)
{
   register int i;

   for (i=0 ;i<col;i++)
        val1[row][i]=val[i];

return OK;
}

int twoD_array_copy_ldouble (long double **val1,long double *val, int row,int col)
{
   register int i;

   for (i=0 ;i<col;i++)
        val1[row][i]=val[i];

return OK;
}


int vector_copy (int *val1,int *val, int numberofvalues)
{
   register int i;

   for (i=0 ;i<numberofvalues;i++)
        val1[i]=val[i];

return OK;
}


long double *vector_calloc(int col)
{
    long double *input;

    input=(long double *)malloc( sizeof(long double) * col);

    var_init(input,col);

    return input;
}


int *vector_calloc_int(int col)
{
    int *input;

    input=(int *)malloc( sizeof(int) * col);

    var_init(input,col);

    return input;
}

char *vector_calloc_char(int col)
{
    char *input;

    input=(char *)malloc( sizeof(char) * col);

    var_init(input,col);

    return input;
}

int matrix_free(long double **input,int row)
{
    register int i;


    for(i=0; i<row; i++ )
        free(input[i]);
    free(input);

    return OK;
}

int matrix_free(int **input,int row)
{
    register int i;


    for(i=0; i<row; i++ )
        free(input[i]);
    free(input);

    return OK;
}


int matrix_free(long double ***input,int row, int col)
{
    register int i,j;


    for(i=0; i<row; i++ )
        for (j=0; j<col;j++)
        free(input[i][j]);
    free(input);

    return OK;
}
int filestore_ber(long double ber,long double ldpc_ber, char *filename)
{

	FILE *fp;
	fp=fopen(filename,"aw");
	fprintf(fp,"%Le \t%Le \n",ber,ldpc_ber);
	fclose(fp);

	return OK;
}


//int filestore_ber(long double SNR,long double ber, char *filename)
//{
//
//	FILE *fp;
//	fp=fopen(filename,"aw");
//	fprintf(fp,"%Lf \t%Le \n",SNR,ber);
//	fclose(fp);
//
//	return OK;
//}


int filestore_prob(long double **input, char *filename, int rows)
{

	FILE *fp;
	fp=fopen(filename,"aw");
	for (int i=0;i<rows;i++)
	{
	    fprintf(fp,"%Le \t %Le \n",input[i][0],input[i][1]);
	}

	fclose(fp);

	return OK;
}


int filestore_ber(long double ber, char *filename)
{

	FILE *fp;
	fp=fopen(filename,"aw");
	fprintf(fp,"%Le \n",ber);
	fclose(fp);

	return OK;
}

