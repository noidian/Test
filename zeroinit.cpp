
#include "init.h"

int var_init(int *var,int k)
{
	register int i;
	for(i=0;i<k;i++)
	{
		var[i]=0;
	}

	return OK;
}

int var_init(char *var,int k)
{
	register int i;
	for(i=0;i<k;i++)
	{
		var[i]=0;
	}

	return OK;
}



int var_init(long double *var,int k)
{
	register int i;
	for(i=0;i<k;i++)
	{
		var[i]=0;
	}
	return OK;
}

int var_init(long double **var,int row, int col)
{
	register int i,j;
	for(i=0;i<row;i++)
        for (j=0;j<col;j++)
            var[i][j]=0;
return OK;
}

int var_init(long double ***var,int x, int y, int z)
{
	register int i,j,k;
	for(i=0;i<x;i++)
        for (j=0;j<y;j++)
            for (k=0;k<z;k++)
                var[i][j][k]=0;
return OK;
}


int var_init(int **var,int row,int col)
{
	register int i,j;
	for(i=0;i<row;i++)
        for (j=0;j<col;j++)
            var[i][j]=0;
return OK;
}
