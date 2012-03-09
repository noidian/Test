#include "comm_func.h"
#include "math_func.h"



int convolution(long double *signal1,long double *signal2,int length1,int length2,long double *output)
{
    register int i,j;

    for (i=0;i<length1+length2;i++)
    {
            for (j=0;j<length2;j++)
            {
              if((j>=0 && i-j>=0)&&i<length1) //takes care of intial bits
              {
                  output[i]=output[i]+ (signal2[j]*signal1[i-j]);
              }
              else if (i>=length1 && j>(i-length1))     //last bits
                    {
                       output[i]=output[i]+ (signal2[j]*signal1[i-j]);
                    }
            }
    }
return OK;
}


int convolution(int *signal1,long double *signal2,int length1,int length2,long double *output)
{
    register int i,j;

    for (i=0;i<length1+length2;i++)
    {
            for (j=0;j<length2;j++)
            {
              if((j>=0 && i-j>=0)&&i<length1) //takes care of intial bits
              {
                  output[i]=output[i]+ (signal2[j]*signal1[i-j]);
              }
              else if (i>=length1 && j>(i-length1))     //last bits
                    {
                       output[i]=output[i]+ (signal2[j]*signal1[i-j]);
                    }
            }
    }
return OK;
}

int input_ITI_adder(long double *output,int *adj_OD,long double *filter_n1,int *main,long double *filter_0, int *adj_ID, long double *filter_1,int Nb_samples,int Nb_filter)
{
  long double *data_main,*data_adj1_OD,*data_adj1_ID;

    /**********Adj1 OD track******************************/
    //Input Filter output = sample length + filter length
    data_adj1_OD = vector_calloc(Nb_samples+Nb_filter);

    convolution(adj_OD,filter_n1,Nb_samples,Nb_filter,data_adj1_OD);

    //truncate initial samples and end values from convolution process
    array_resize(data_adj1_OD,Nb_filter-1,Nb_samples);
    /**********Adj1 OD track******************************/

    /**********Main track******************************/
    //Input Filter output = sample length + filter length
    data_main = vector_calloc(Nb_samples+Nb_filter);

    convolution(main,filter_0,Nb_samples,Nb_filter,data_main);

    //truncate initial samples and end values from convolution process
    array_resize(data_main,Nb_filter-1,Nb_samples);
    /**********Main track******************************/

    /**********Adj1 ID track******************************/
    //Input Filter output = sample length + filter length
    data_adj1_ID = vector_calloc(Nb_samples+Nb_filter);

    convolution(adj_ID,filter_1,Nb_samples,Nb_filter,data_adj1_ID);

    //truncate initial samples and end values from convolution process
    array_resize(data_adj1_ID,Nb_filter-1,Nb_samples);
    /**********Adj1 ID track******************************/

   //Length of sample is sample_length-(Filter_coeff-1);
    Nb_samples=Nb_samples-(Nb_filter-1);

    array_addition(output,data_main,data_adj1_OD,Nb_samples);
    array_addition(output,output,data_adj1_ID,Nb_samples);


    free(data_main);
    free(data_adj1_ID);
    free(data_adj1_OD);

    return OK;
}

//Notation same as 2D equalizer
int lagrange_err_min_1D(long double *equalizer,int no_equalizer_coeff,long double *target,int no_of_target_vals,int *a,int numberofvalues,long double *matched_main)
{
    int i,j;
    int track_count;
    int P1,P2,K1,K2,N1,N2,M1,M2,L,K;

    int *a_merge;
    long double *y_merge;
    long double normalize;

    //matrix parameters
    int main_rows_A,main_cols_A,sub_rows_A,sub_cols_A;
    int correlation_values,column_count,cshift;

    long double *A0,**A;

    int main_rows_T,main_cols_T,sub_rows_T,sub_cols_T;

    long double *T0,**T,**T_transpose;


    int main_rows_R,main_cols_R,sub_rows_R,sub_cols_R;

    long double *R0,**R,**R_inverse;


    long double **temp,**temp1,**temp2,**temp_prod;
    long double *temp3;

    long double *E_i,lambda;




    //Equalizer coefficients going from -K1 to K2.
    if (no_equalizer_coeff%2==0)
    {
        K1=(no_equalizer_coeff/2)-1;
        K2=(no_equalizer_coeff/2);
    }
    else
    {
        K1=(no_equalizer_coeff/2);
        K2=(no_equalizer_coeff/2);
    }

    //Tracks to consider for GPR,  P1 adjacent tracks towards ID side, P2 is adjacent tracks towards OD side
    //1D GPR no adjacent tracks
    P1=0;
    P2=0;


    //number of Target values if [0 1 2] then corresponds to -1 to 1
    if (no_of_target_vals%2==0)
    {
        M1=(no_of_target_vals/2)-1;
        M2=(no_of_target_vals/2);
    }
    else
    {
        M1=(no_of_target_vals/2);
        M2=(no_of_target_vals/2);
    }

    //Number of tracks considered for target computation
    N1=0;
    N2=0;

    //Target Length
    L=M1+M2+1;

    //Equalizer Length
    K=K1+K2+1;

    //Tabulate input a into a common matrix (notation same as 2D)


    a_merge=a;

 /****Note on Normalization****/
 //Normalizing zero lag to '1' is biased and provides poor results
 //Normalizing by computing the mean and variance is not optimal either
//Normalizing using track length is rather sub optimal
//Normalizing using equal energy norms is to a small extent sub optimal
//Normalize using corresponding norm seems to have better performance

//anamoly if normalize R matrix with normalization factor of T and normalize A with the same


    y_merge=matched_main;

//Test
normalize=norm(a_merge,numberofvalues)*norm(y_merge,numberofvalues);
    //normalize=norm(a_merge,numberofvalues)*norm(a_merge,numberofvalues);

    track_count=1;          //relevant only for 2D

    main_rows_A=abs((N1+N2))+1;
    sub_rows_A=abs((M2+M1))+1;

    main_cols_A=abs((-N1-N2))+1;
    sub_cols_A=abs((-M2-M1))+1;


    correlation_values=2*(numberofvalues-1); //includes 0
    //A0 = new long double [correlation_values];

    //Can use calloc but not supported on some compilers
    A0=vector_calloc(correlation_values);

    if (main_cols_A*sub_cols_A>main_rows_A*sub_rows_A)
    {
            cshift=(main_cols_A*sub_cols_A)+1;
    }
    else
    {
            cshift=(main_rows_A*sub_rows_A)+1;
    }


    x_corr(a_merge, a_merge, cshift, numberofvalues, A0);
    array_divide(A0,normalize,correlation_values);

    //A matrix formation
    A=matrix_calloc((main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A));

    for (i=0;i<(main_rows_A*sub_rows_A);i++)
    {
        for (j=0;j<(main_cols_A*sub_cols_A);j++)
        {
            column_count=(numberofvalues-1)-j+i;
            A[i][j]=A0[column_count];
        }
    }





//    normalize=norm(y_merge,numberofvalues)*norm(a_merge,numberofvalues);

    main_rows_T=(P2+N1)-(-P1+N1)+1;
    sub_rows_T=(K2+M1)-(-K1+M1)+1;

    main_cols_T=abs((-P1-N2)-(-P1+N1))+1;
    sub_cols_T=abs((-K1-M2)-(-K1+M1))+1;


    //Resuse A0 memory
    T0 = A0;
    var_init(T0,correlation_values);

   if (main_cols_T*sub_cols_T>main_rows_T*sub_rows_T)
    {
            cshift=(main_cols_T*sub_cols_T)+1;
    }
    else
    {
            cshift=(main_rows_T*sub_rows_T)+1;
    }

    x_corr(a_merge, y_merge, cshift, numberofvalues, T0);
    array_divide(T0,normalize,correlation_values);


     //T matrix formation
    T=matrix_calloc((main_rows_T*sub_rows_T ),(main_cols_T*sub_cols_T));

    for (i=0;i<(main_rows_T*sub_rows_T);i++)
    {
        for (j=0;j<(main_cols_T*sub_cols_T);j++)
        {
            column_count=(numberofvalues-1)-K1+M1-(j%L)+(i%K);
            T[i][j]=T0[column_count];
        }
    }

    //R Matrix formulation


//Test
//normalize=norm(a_merge,numberofvalues)*norm(y_merge,numberofvalues);

    //normalize=norm(y_merge,numberofvalues)*norm(y_merge,numberofvalues);

    main_rows_R=(P2+P1)+1;
    sub_rows_R=(K2+K1)+1;

    main_cols_R=abs((-P1-P2))+1;
    sub_cols_R=abs((-K2-K1))+1;

    //Resuse A0 memory
    R0 = T0;
    var_init(R0,correlation_values);

    if (main_cols_R*sub_cols_R>main_rows_R*sub_rows_R)
    {
            cshift=(main_cols_R*sub_cols_R)+1;
    }
    else
    {
            cshift=(main_rows_R*sub_rows_R)+1;
    }

    x_corr(y_merge, y_merge, cshift, numberofvalues, R0);
    array_divide(R0,normalize,correlation_values);


     //R matrix formation
     R=matrix_calloc((main_rows_R*sub_rows_R ),(main_cols_R*sub_cols_R));

    for (i=0;i<(main_rows_R*sub_rows_R);i++)
    {
        for (j=0;j<(main_cols_R*sub_cols_R);j++)
        {
            column_count=(numberofvalues-1)-(j%K)+(i%K);
            R[i][j]=R0[column_count];
        }
    }


    //R_inverse
    R_inverse=matrix_calloc((main_rows_R*sub_rows_R ),(main_cols_R*sub_cols_R));

    matrix_inverse (R_inverse,(main_rows_R*sub_rows_R),(main_cols_R*sub_cols_R),R);

     //T' matrix formation
    T_transpose=matrix_calloc((main_cols_T*sub_cols_T),(main_rows_T*sub_rows_T) );


    matrix_transpose(T_transpose,main_rows_T*sub_rows_T,main_cols_T*sub_cols_T,T);

    //Compute T'*R
    temp=matrix_calloc((main_cols_T*sub_cols_T),(main_cols_R*sub_cols_R) );

    matrix_multiply(temp,T_transpose,(main_cols_T*sub_cols_T),(main_rows_T*sub_rows_T),R_inverse,(main_rows_R*sub_rows_R ),(main_cols_R*sub_cols_R));

    //Compute temp*T
    temp1=matrix_calloc((main_cols_T*sub_cols_T),(main_cols_T*sub_cols_T) );

    matrix_multiply(temp1,temp,(main_cols_T*sub_cols_T),(main_cols_R*sub_cols_R),T,(main_rows_T*sub_rows_T ),(main_cols_T*sub_cols_T));

    //Compute A- temp1
    temp2=matrix_calloc((main_cols_A*sub_cols_A),(main_cols_A*sub_cols_A) );

    matrix_subtraction(temp2,A,temp1,(main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A));

    //Compute inv(A-T'R_inv*T)
     //temp_prod matrix formation
    temp_prod=matrix_calloc((main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A) );
    matrix_inverse (temp_prod,(main_rows_A*sub_rows_A),(main_cols_A*sub_cols_A),temp2);

    //free memory for temporary variables
    matrix_free(temp,(main_cols_T*sub_cols_T));
    matrix_free(temp1,(main_cols_T*sub_cols_T));
    matrix_free(temp2,(main_cols_T*sub_cols_T));


    //Defines contraint matrix "i"
    E_i=vector_calloc(L);
    E_i[int(L/2)]=1;

    //i'*temp_prod
    temp3=vector_calloc((main_rows_A*sub_rows_A ));
    matrix_multiply(temp3,E_i,L,temp_prod,(main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A));

    lambda=0;
    //temp1=(temp*i), inv(temp1)
    for (i=0;i<L;i++)
        lambda+=temp3[i]*E_i[i];

    lambda=1/lambda;

    //temp3=E_i
    var_init(temp3,L);
    temp3[(int(L/2))]=1;

    array_multiply(temp3,lambda,L);
    //G Matrix = target


    matrix_multiply(target,temp_prod,(main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A),temp3,L);

    //Free temporary storage
    free(temp3);

    //equalizer=F matrix
    temp3=vector_calloc((main_rows_T*sub_rows_T));

    matrix_multiply(temp3,T,(main_rows_T*sub_rows_T),(main_cols_T*sub_cols_T),target,L);

    matrix_multiply(equalizer,R_inverse,(main_rows_R*sub_rows_R),(main_cols_R*sub_cols_R),temp3,(main_rows_T*sub_rows_T));

    //Free temporary storage
    free(temp3);

    free(A0);

    matrix_free(A,(main_rows_A*sub_rows_A));

    matrix_free(T,(main_rows_T*sub_rows_T));

    matrix_free(R,(main_rows_R*sub_rows_R));

    matrix_free(R_inverse,(main_rows_R*sub_rows_R));

    matrix_free(T_transpose,(main_cols_T*sub_cols_T ));


return OK;
}


int lagrange_err_min_2D(long double *f_0,long double *f_n1,long double *f_1,int Nb_eq_0,int Nb_eq_n1,int Nb_eq_1,long double *G,int Nb_target,int **a_merge,int Nb_samples,long double *matched_main,long double *matched_adj1_OD,long double *matched_adj1_ID)
{
        int i,j;
    int track_count;
    int P1,P2,K1,K2,N1,N2,M1,M2,L,K;

    long double **y_merge;
    long double normalize;

    //matrix parameters
    int main_rows_A,main_cols_A,sub_rows_A,sub_cols_A;
    int correlation_values,corr_rows,column_count,cshift,rshift;

    int start_row,c_row_count,start_val;

    long double **A0,**A;

    int main_rows_T,main_cols_T,sub_rows_T,sub_cols_T;

    long double **T0,**T,**T_transpose;


    int main_rows_R,main_cols_R,sub_rows_R,sub_cols_R;

    long double **R0,**R,**R_inverse;


    long double **temp,**temp1,**temp2,**temp_prod;
    long double *temp3;

    long double *E_i,lambda;

    long double *target,*equalizer;

//Assuming all tracks have equal coeeficients

    //Equalizer coefficients going from -K1 to K2.
    if (Nb_eq_0%2==0)
    {
        K1=(Nb_eq_0/2)-1;
        K2=(Nb_eq_0/2);
    }
    else
    {
        K1=(Nb_eq_0/2);
        K2=(Nb_eq_0/2);
    }

    //Tracks to consider for GPR,  P1 adjacent tracks towards ID side, P2 is adjacent tracks towards OD side
    P1=1;
    P2=1;


    //number of Target values if [0 1 2] then corresponds to -1 to 1
    if (Nb_target%2==0)
    {
        M1=(Nb_target/2)-1;
        M2=(Nb_target/2);
    }
    else
    {
        M1=(Nb_target/2);
        M2=(Nb_target/2);
    }

    //Number of tracks considered for target computation
    N1=0;
    N2=0;

    //Target Length
    L=M1+M2+1;

    //Equalizer Length
    K=K1+K2+1;

    //Nb of tracks
    track_count=P2-(-P1)+1;





 /****Note on Normalization****/
 //Normalizing zero lag to '1' is biased and provides poor results
 //Normalizing by computing the mean and variance is not optimal either
//Normalizing using track length is rather sub optimal
//Normalizing using equal energy norms is to a small extent sub optimal
//Normalize using corresponding norm seems to have better performance

//anamoly if normalize R matrix with normalization factor of T and normalize A with the same


    //Tabulate input a into a common matrix
    y_merge=matrix_calloc(3,Nb_samples);
    twoD_array_copy_ldouble(y_merge,matched_adj1_OD,0,Nb_samples);
    twoD_array_copy_ldouble(y_merge,matched_main,1,Nb_samples);
    twoD_array_copy_ldouble(y_merge,matched_adj1_ID,2,Nb_samples);

//Test
    normalize=norm(a_merge,3,Nb_samples)*norm(y_merge,3,Nb_samples); //use a and y to normalize
    //normalize=norm(a_merge,3,Nb_samples)*norm(a_merge,3,Nb_samples);


    main_rows_A=abs((N1+N2))+1;
    sub_rows_A=abs((M2+M1))+1;

    main_cols_A=abs((-N1-N2))+1;
    sub_cols_A=abs((-M2-M1))+1;


    correlation_values=2*(Nb_samples-1); //includes 0 so no need for addtional mid point

    //Can use calloc but not supported on some compilers



    if (main_cols_A*sub_cols_A>main_rows_A*main_rows_A)
    {
            cshift=(main_cols_A*sub_cols_A)+1;
    }
    else
    {
            cshift=(main_rows_A*sub_rows_A)+1;
    }

    rshift=3;

      //Assuming the number of equalizer coefficients is greater than the number of targets....
    A0=matrix_calloc(5,(2*cshift)+1);

    //rows start from "0"so 0-2 corresponds to 3
    //Nb_sample-1 since we start from 0
    x_corr2(a_merge,2,Nb_samples-1, a_merge,2,Nb_samples-1,rshift,cshift, Nb_samples-1, A0);

    corr_rows=5; //rows1=2 rows2=2 i.e. 0-2 0-2 corresponds to 5
    //normalize total of 2wice + "0th compoenent
    array_divide(A0,normalize,corr_rows,(2*cshift)+1);

    //A matrix formation
    A=matrix_calloc((main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A));

    for (i=0;i<(main_rows_A*sub_rows_A);i++)
    {
        start_row=(track_count-1)+floor(i/L);       //track_count-1 since tracks start from 0
        for (j=0;j<(main_cols_A*sub_cols_A);j++)
        {
            column_count=(cshift)-(j%L)+(i%L);
            c_row_count=start_row-floor(j/L);
                    A[i][j]=A0[c_row_count][column_count];
        }
    }



//    normalize=norm(y_merge,3,Nb_samples)*norm(a_merge,3,Nb_samples);

    main_rows_T=(P2+N1)-(-P1+N1)+1;
    sub_rows_T=(K2+M1)-(-K1+M1)+1;

    main_cols_T=abs((-P1-N2)-(-P1+N1))+1;
    sub_cols_T=abs((-K1-M2)-(-K1+M1))+1;




    if (main_cols_T*sub_cols_T>main_rows_T*main_rows_T)
    {
            cshift=(main_cols_T*sub_cols_T)+1;
    }
    else
    {
            cshift=(main_rows_T*sub_rows_T)+1;
    }
    rshift=3;

    T0 = matrix_calloc(5,(2*cshift)+1);

        //rows start from "0"so 0-2 corresponds to 3
    //Nb_sample-1 since we start from 0

    x_corr2(a_merge,2,Nb_samples-1, y_merge,2,Nb_samples-1,rshift,cshift, Nb_samples-1, T0);


    corr_rows=5; //rows1=2 rows2=2 i.e. 0-2 0-2 corresponds to 5

    array_divide(T0,normalize,corr_rows,(2*cshift)+1);

    start_val=-P1+N1;
     //T matrix formation
    T=matrix_calloc((main_rows_T*sub_rows_T ),(main_cols_T*sub_cols_T));

    for (i=0;i<(main_rows_T*sub_rows_T);i++)
    {
        start_row=(track_count-1)+floor(i/K)+start_val;
        for (j=0;j<(main_cols_T*sub_cols_T);j++)
        {
            c_row_count=start_row-floor(j/L);
            column_count=cshift-K1+M1-(j%L)+(i%K);
            T[i][j]=T0[c_row_count][column_count];
        }
    }


    //R Matrix formulation

    //Test
   // normalize=norm(a_merge,3,Nb_samples)*norm(y_merge,3,Nb_samples);

    //normalize=norm(y_merge,3,Nb_samples)*norm(y_merge,3,Nb_samples);


    main_rows_R=(P2+P1)+1;
    sub_rows_R=(K2+K1)+1;

    main_cols_R=abs((-P1-P2))+1;
    sub_cols_R=abs((-K2-K1))+1;

    if (main_cols_R*sub_cols_R>main_rows_R*main_rows_R)
    {
            cshift=(main_cols_R*sub_cols_R)+1;
    }
    else
    {
            cshift=(main_rows_R*sub_rows_R)+1;
    }
    rshift=3;


    R0 = matrix_calloc(5,(2*cshift)+1);

       //rows start from "0"so 0-2 corresponds to 3
    //Nb_sample-1 since we start from 0

    x_corr2(y_merge,2,Nb_samples-1, y_merge,2,Nb_samples-1,rshift,cshift, Nb_samples-1, R0);

    //free memory
    matrix_free(y_merge,3);


    corr_rows=5; //rows1=2 rows2=2 i.e. 0-2 0-2 corresponds to 5

    array_divide(R0,normalize,corr_rows,(2*cshift)+1);

     //R matrix formation
     R=matrix_calloc((main_rows_R*sub_rows_R ),(main_cols_R*sub_cols_R));

    for (i=0;i<(main_rows_R*sub_rows_R);i++)
    {
        start_row=(track_count-1)+floor(i/K);
        for (j=0;j<(main_cols_R*sub_cols_R);j++)
        {
            c_row_count=start_row-floor(j/K);
            column_count=cshift-(j%K)+(i%K);
            R[i][j]=R0[c_row_count][column_count];
        }
    }

    //R_inverse
    R_inverse=matrix_calloc((main_rows_R*sub_rows_R ),(main_cols_R*sub_cols_R));

    matrix_inverse (R_inverse,(main_rows_R*sub_rows_R),(main_cols_R*sub_cols_R),R);

     //T' matrix formation
    T_transpose=matrix_calloc((main_cols_T*sub_cols_T),(main_rows_T*sub_rows_T) );


    matrix_transpose(T_transpose,main_rows_T*sub_rows_T,main_cols_T*sub_cols_T,T);

    //Compute T'*R_inverse
    temp=matrix_calloc((main_cols_T*sub_cols_T),(main_cols_R*sub_cols_R) );

    matrix_multiply(temp,T_transpose,(main_cols_T*sub_cols_T),(main_rows_T*sub_rows_T),R_inverse,(main_rows_R*sub_rows_R ),(main_cols_R*sub_cols_R));

    //Compute temp*T
    temp1=matrix_calloc((main_cols_T*sub_cols_T),(main_cols_T*sub_cols_T) );

    matrix_multiply(temp1,temp,(main_cols_T*sub_cols_T),(main_cols_R*sub_cols_R),T,(main_rows_T*sub_rows_T ),(main_cols_T*sub_cols_T));

    //Compute A- temp1
    temp2=matrix_calloc((main_cols_A*sub_cols_A),(main_cols_A*sub_cols_A) );

    matrix_subtraction(temp2,A,temp1,(main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A));

    //Compute inv(A-T'R_inv*T)
     //temp_prod matrix formation
    temp_prod=matrix_calloc((main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A) );
    matrix_inverse (temp_prod,(main_rows_A*sub_rows_A),(main_cols_A*sub_cols_A),temp2);

    //free memory for temporary variables
    matrix_free(temp,(main_cols_T*sub_cols_T));
    matrix_free(temp1,(main_cols_T*sub_cols_T));
    matrix_free(temp2,(main_cols_A*sub_cols_A));


    //Defines contraint matrix "i"
    E_i=vector_calloc(L);
    E_i[int(L/2)]=1;

    //i'*temp_prod
    temp3=vector_calloc((main_rows_A*sub_rows_A ));
    matrix_multiply(temp3,E_i,L,temp_prod,(main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A));

    lambda=0;
    //temp1=(temp*i), inv(temp1)
    for (i=0;i<L;i++)
        lambda+=temp3[i]*E_i[i];

    lambda=1/lambda;

    //temp3=E_i
    var_init(temp3,L);
    temp3[(int(L/2))]=1;

    array_multiply(temp3,lambda,L);
    //G Matrix = target

    //Works when all threee equalizers are same width
    target=vector_calloc(L);
    equalizer=vector_calloc(track_count*K);

    matrix_multiply(target,temp_prod,(main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A),temp3,L);

    //Free temporary storage
    free(temp3);

    //equalizer=F matrix
    temp3=vector_calloc((main_rows_T*sub_rows_T));

    matrix_multiply(temp3,T,(main_rows_T*sub_rows_T),(main_cols_T*sub_cols_T),target,L);

    matrix_multiply(equalizer,R_inverse,(main_rows_R*sub_rows_R),(main_cols_R*sub_cols_R),temp3,(main_rows_T*sub_rows_T));

    //Resassign each equalizer to respective tracks

    array_copy(f_n1,equalizer,0,(Nb_eq_n1));
    array_copy(f_0,equalizer,Nb_eq_n1,(Nb_eq_n1+Nb_eq_0));
    array_copy(f_1,equalizer,(Nb_eq_n1+Nb_eq_0),(Nb_eq_n1+Nb_eq_0+Nb_eq_1));
    array_copy(G,target,3);

    //Free temporary storage

    free(temp3);

    free(target);
    free(equalizer);

    matrix_free(temp_prod,(main_rows_A*sub_rows_A));

    matrix_free(A0,corr_rows);
    matrix_free(T0,corr_rows);
    matrix_free(R0,corr_rows);

    free(E_i);

    matrix_free(A,(main_rows_A*sub_rows_A));

    matrix_free(T,(main_rows_T*sub_rows_T));

    matrix_free(R,(main_rows_R*sub_rows_R));

    matrix_free(R_inverse,(main_rows_R*sub_rows_R));

    matrix_free(T_transpose,(main_cols_T*sub_cols_T ));

    return OK;

}


int lagrange_err_min_onesided(long double *f_0,long double *f_n1,int Nb_eq_0,int Nb_eq_n1,long double *G,int Nb_target,int **a_merge,int Nb_samples,long double *matched_main,long double *matched_adj1_OD)
{
        int i,j;
    int track_count;
    int P1,P2,K1,K2,N1,N2,M1,M2,L,K;

    long double **y_merge;
    long double normalize;

    //matrix parameters
    int main_rows_A,main_cols_A,sub_rows_A,sub_cols_A;
    int correlation_values,corr_rows,column_count,cshift,rshift;

    int start_row,c_row_count,start_val;

    long double **A0,**A;

    int main_rows_T,main_cols_T,sub_rows_T,sub_cols_T;

    long double **T0,**T,**T_transpose;


    int main_rows_R,main_cols_R,sub_rows_R,sub_cols_R;

    long double **R0,**R,**R_inverse;


    long double **temp,**temp1,**temp2,**temp_prod;
    long double *temp3;

    long double *E_i,lambda;

    long double *target,*equalizer;

//Assuming all tracks have equal coeeficients

    //Equalizer coefficients going from -K1 to K2.
    if (Nb_eq_0%2==0)
    {
        K1=(Nb_eq_0/2)-1;
        K2=(Nb_eq_0/2);
    }
    else
    {
        K1=(Nb_eq_0/2);
        K2=(Nb_eq_0/2);
    }

    //Tracks to consider for GPR,
    //For Onesided only P1 is taken P2=0
    P1=1;
    P2=0;


    //number of Target values if [0 1 2] then corresponds to -1 to 1
    if (Nb_target%2==0)
    {
        M1=(Nb_target/2)-1;
        M2=(Nb_target/2);
    }
    else
    {
        M1=(Nb_target/2);
        M2=(Nb_target/2);
    }

    //Number of tracks considered for target computation
    N1=0;
    N2=0;

    //Target Length
    L=M1+M2+1;

    //Equalizer Length
    K=K1+K2+1;

    //Nb of tracks
    track_count=P2-(-P1)+1;





 /****Note on Normalization****/
 //Normalizing zero lag to '1' is biased and provides poor results
 //Normalizing by computing the mean and variance is not optimal either
//Normalizing using track length is rather sub optimal
//Normalizing using equal energy norms is to a small extent sub optimal
//Normalize using corresponding norm seems to have better performance

//anamoly if normalize R matrix with normalization factor of T and normalize A with the same

   //Tabulate input a into a common matrix
    y_merge=matrix_calloc(track_count,Nb_samples);
    twoD_array_copy_ldouble(y_merge,matched_adj1_OD,0,Nb_samples);
    twoD_array_copy_ldouble(y_merge,matched_main,1,Nb_samples);


    normalize=norm(a_merge,track_count,Nb_samples)*norm(y_merge,track_count,Nb_samples);


    main_rows_A=abs((N1+N2))+1;
    sub_rows_A=abs((M2+M1))+1;

    main_cols_A=abs((-N1-N2))+1;
    sub_cols_A=abs((-M2-M1))+1;


    correlation_values=2*(Nb_samples-1); //includes 0 so no need for addtional mid point

    //Can use calloc but not supported on some compilers



    if (main_cols_A*sub_cols_A>main_rows_A*main_rows_A)
    {
            cshift=(main_cols_A*sub_cols_A)+1;
    }
    else
    {
            cshift=(main_rows_A*sub_rows_A)+1;
    }

    rshift=track_count;

      //Assuming the number of equalizer coefficients is greater than the number of targets....
    A0=matrix_calloc((2*track_count)-1,(2*cshift)+1);

    //rows start from "0"so 0-2 corresponds to 3
    //Nb_sample-1 since we start from 0
    x_corr2(a_merge,1,Nb_samples-1, a_merge,1,Nb_samples-1,rshift,cshift, Nb_samples-1, A0);

    corr_rows=3; //rows1=1 rows2=1 i.e. 0-1 0-1 corresponds to 3
    //normalize total of 2wice + "0th compoenent
    array_divide(A0,normalize,corr_rows,(2*cshift)+1);

    //A matrix formation
    A=matrix_calloc((main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A));

    for (i=0;i<(main_rows_A*sub_rows_A);i++)
    {
        start_row=(track_count-1)+floor(i/L);       //track_count-1 since tracks start from 0
        for (j=0;j<(main_cols_A*sub_cols_A);j++)
        {
            column_count=(cshift)-(j%L)+(i%L);
            c_row_count=start_row-floor(j/L);
                    A[i][j]=A0[c_row_count][column_count];
        }
    }




    //normalize=norm(y_merge,track_count,Nb_samples)*norm(a_merge,track_count,Nb_samples);

    main_rows_T=(P2+N1)-(-P1+N1)+1;
    sub_rows_T=(K2+M1)-(-K1+M1)+1;

    main_cols_T=abs((-P1-N2)-(-P1+N1))+1;
    sub_cols_T=abs((-K1-M2)-(-K1+M1))+1;




    if (main_cols_T*sub_cols_T>main_rows_T*main_rows_T)
    {
            cshift=(main_cols_T*sub_cols_T)+1;
    }
    else
    {
            cshift=(main_rows_T*sub_rows_T)+1;
    }
    rshift=track_count;

    T0 = matrix_calloc((2*track_count)-1,(2*cshift)+1);

        //rows start from "0"so 0-2 corresponds to 3
    //Nb_sample-1 since we start from 0

    x_corr2(a_merge,1,Nb_samples-1, y_merge,1,Nb_samples-1,rshift,cshift, Nb_samples-1, T0);


    corr_rows=(2*track_count)-1; //rows1=1 rows2=1 i.e. 0-1 0-1 corresponds to 4-1=3

    array_divide(T0,normalize,corr_rows,(2*cshift)+1);

    start_val=-P1+N1;
     //T matrix formation
    T=matrix_calloc((main_rows_T*sub_rows_T ),(main_cols_T*sub_cols_T));

    for (i=0;i<(main_rows_T*sub_rows_T);i++)
    {
        start_row=(track_count-1)+floor(i/K)+start_val;
        for (j=0;j<(main_cols_T*sub_cols_T);j++)
        {
            c_row_count=start_row-floor(j/L);
            column_count=cshift-K1+M1-(j%L)+(i%K);
            T[i][j]=T0[c_row_count][column_count];
        }
    }


    //R Matrix formulation

   // normalize=norm(y_merge,track_count,Nb_samples)*norm(a_merge,track_count,Nb_samples);


    main_rows_R=(P2+P1)+1;
    sub_rows_R=(K2+K1)+1;

    main_cols_R=abs((-P1-P2))+1;
    sub_cols_R=abs((-K2-K1))+1;

    if (main_cols_R*sub_cols_R>main_rows_R*main_rows_R)
    {
            cshift=(main_cols_R*sub_cols_R)+1;
    }
    else
    {
            cshift=(main_rows_R*sub_rows_R)+1;
    }
    rshift=track_count;


    R0 = matrix_calloc((2*track_count)-1,(2*cshift)+1);

       //rows start from "0"so 0-1 corresponds to 2
    //Nb_sample-1 since we start from 0

    x_corr2(y_merge,1,Nb_samples-1, y_merge,1,Nb_samples-1,rshift,cshift, Nb_samples-1, R0);

    //free memory
    matrix_free(y_merge,2);


    corr_rows=(2*track_count)-1; //rows1=2 rows2=2 i.e. 0-2 0-2 corresponds to 5

    array_divide(R0,normalize,corr_rows,(2*cshift)+1);

     //R matrix formation
     R=matrix_calloc((main_rows_R*sub_rows_R ),(main_cols_R*sub_cols_R));

    for (i=0;i<(main_rows_R*sub_rows_R);i++)
    {
        start_row=(track_count-1)+floor(i/K);
        for (j=0;j<(main_cols_R*sub_cols_R);j++)
        {
            c_row_count=start_row-floor(j/K);
            column_count=cshift-(j%K)+(i%K);
            R[i][j]=R0[c_row_count][column_count];
        }
    }

    //R_inverse
    R_inverse=matrix_calloc((main_rows_R*sub_rows_R ),(main_cols_R*sub_cols_R));

    matrix_inverse (R_inverse,(main_rows_R*sub_rows_R),(main_cols_R*sub_cols_R),R);

     //T' matrix formation
    T_transpose=matrix_calloc((main_cols_T*sub_cols_T),(main_rows_T*sub_rows_T) );


    matrix_transpose(T_transpose,main_rows_T*sub_rows_T,main_cols_T*sub_cols_T,T);

    //Compute T'*R_inverse
    temp=matrix_calloc((main_cols_T*sub_cols_T),(main_cols_R*sub_cols_R) );

    matrix_multiply(temp,T_transpose,(main_cols_T*sub_cols_T),(main_rows_T*sub_rows_T),R_inverse,(main_rows_R*sub_rows_R ),(main_cols_R*sub_cols_R));

    //Compute temp*T
    temp1=matrix_calloc((main_cols_T*sub_cols_T),(main_cols_T*sub_cols_T) );

    matrix_multiply(temp1,temp,(main_cols_T*sub_cols_T),(main_cols_R*sub_cols_R),T,(main_rows_T*sub_rows_T ),(main_cols_T*sub_cols_T));

    //Compute A- temp1
    temp2=matrix_calloc((main_cols_A*sub_cols_A),(main_cols_A*sub_cols_A) );

    matrix_subtraction(temp2,A,temp1,(main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A));

    //Compute inv(A-T'R_inv*T)
     //temp_prod matrix formation
    temp_prod=matrix_calloc((main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A) );
    matrix_inverse (temp_prod,(main_rows_A*sub_rows_A),(main_cols_A*sub_cols_A),temp2);

    //free memory for temporary variables
    matrix_free(temp,(main_cols_T*sub_cols_T));
    matrix_free(temp1,(main_cols_T*sub_cols_T));
    matrix_free(temp2,(main_cols_A*sub_cols_A));


    //Defines contraint matrix "i"
    E_i=vector_calloc(L);
    E_i[int(L/2)]=1;

    //i'*temp_prod
    temp3=vector_calloc((main_rows_A*sub_rows_A ));
    matrix_multiply(temp3,E_i,L,temp_prod,(main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A));

    lambda=0;
    //temp1=(temp*i), inv(temp1)
    for (i=0;i<L;i++)
        lambda+=temp3[i]*E_i[i];

    lambda=1/lambda;

    //temp3=E_i
    var_init(temp3,L);
    temp3[(int(L/2))]=1;

    array_multiply(temp3,lambda,L);
    //G Matrix = target

    //Works when all threee equalizers are same width
    target=vector_calloc(L);
    equalizer=vector_calloc(track_count*K);

    matrix_multiply(target,temp_prod,(main_rows_A*sub_rows_A ),(main_cols_A*sub_cols_A),temp3,L);

    //Free temporary storage
    free(temp3);

    //equalizer=F matrix
    temp3=vector_calloc((main_rows_T*sub_rows_T));

    matrix_multiply(temp3,T,(main_rows_T*sub_rows_T),(main_cols_T*sub_cols_T),target,L);

    matrix_multiply(equalizer,R_inverse,(main_rows_R*sub_rows_R),(main_cols_R*sub_cols_R),temp3,(main_rows_T*sub_rows_T));

    //Resassign each equalizer to respective tracks

    array_copy(f_n1,equalizer,0,(Nb_eq_n1));
    array_copy(f_0,equalizer,Nb_eq_n1,(Nb_eq_n1+Nb_eq_0));
    array_copy(G,target,3);

    //Free temporary storage

    free(temp3);

    free(target);
    free(equalizer);

    matrix_free(temp_prod,(main_rows_A*sub_rows_A));

    matrix_free(A0,corr_rows);
    matrix_free(T0,corr_rows);
    matrix_free(R0,corr_rows);

    free(E_i);

    matrix_free(A,(main_rows_A*sub_rows_A));

    matrix_free(T,(main_rows_T*sub_rows_T));

    matrix_free(R,(main_rows_R*sub_rows_R));

    matrix_free(R_inverse,(main_rows_R*sub_rows_R));

    matrix_free(T_transpose,(main_cols_T*sub_cols_T ));

    return OK;

}



//General functions used by the communication functions

int awgn(long double *output,long double *input,int Nb_samples,float SNR)
{
    long double sigma,es_over_n0,es,mean,sn_ratio,noise;
    register int i;
    long double u,r;
 //N.B. Es/N0 = Eb/N0 x Rate x log(2)M where Es/N0 is the SNR Rate is FEC rate of coder and Log(2)M for M=8

  es_over_n0=SNR;  //Asuming the input is the Es_over_N0

  es = 1;  //iess normalized to 1

  sn_ratio =  pow(10,( es_over_n0 / 10) );

  sigma =  sqrt (es / ( sn_ratio ) ); // Std deviation fo 8psk

  sigma/=sqrt(2); //for N0/2 for each component and each bit

  mean=0;


  for(i=0;i<Nb_samples;i++)
   {
   /* generate a uniformly distributed random number u between 0 and 1 - 1E-6*/
    u = (long double)rand() / RAND_MAX;
    if (u == 1.0) u = 0.999999999;

    /* generate a Rayleigh-distributed random number r using u */
    r = (sigma)* sqrt( 2.0 * log( 1.0 / (1.0 - u) ) );

    /* generate another uniformly-distributed random number u as before*/
    u = (long double)rand() /RAND_MAX;
    if (u == 1.0) u = 0.999999999;

    /* generate and return a Gaussian-distributed random number using r and u */
    noise=mean + (r * cos(2 * PI * u));

    output[i]=input[i]+noise;
    }

return OK;
}





int int2bin (const int number, int table[], const int length_table)
{

	int cpt;
	int decimal;

	decimal=number;

	//intilisation of tabel to zero
	for (cpt=0;cpt<length_table;cpt++)
	{
		table[cpt]=0;
	}


	for (cpt=0;cpt<length_table;cpt++)
	{
		if (decimal>=(int)pow(2,length_table-cpt-1))
		{
			table[cpt]=1;
			decimal-=(int)pow(2,length_table-cpt-1);
		}
	}

	return 1;
}


int shift_reg (int reg[], const int length_reg, const int new_value)


{
	int cpt;


	for (cpt=length_reg-1;cpt>0;cpt--)
	{
		reg[cpt]=reg[cpt-1];
	}
	reg[0]=new_value;
	return 1;
}


int bin2int (int *number, const int table[], const int length_table)


{
	int cpt;
	int decimal;

	decimal=0;
	for (cpt=0;cpt<length_table;cpt++)
	{
		decimal+=(int)(table[cpt]*pow(2,(length_table-cpt-1)));
	}
	*number=decimal;
	return 1;
}


int encode_bit (const int Polynom1, const int Polynom2,const int shift_register[], const int length_register, double *out_poly1, double *out_poly2)


{
	int cpt,temp1,temp2;
	int *poly1,*poly2;

	poly1= (int*)malloc(length_register*sizeof(int));
	poly2= (int*)malloc(length_register*sizeof(int));


	int2bin (Polynom1, poly1, length_register);
	int2bin (Polynom2, poly2, length_register);


	// init the outputs

	temp1=0;
	temp2=0;


	// compute the xor

	for (cpt=0;cpt<length_register;cpt++)

	{
		temp1^= ( poly1[cpt] && shift_register[cpt]) ;
		temp2^= ( poly2[cpt] && shift_register[cpt]) ;


	}

	*out_poly1=(double)((1-2*temp1));
	*out_poly2=(double)((1-2*temp2));


	free(poly1);
	free(poly2);


	return 1;
}




int gene_viterbi_matrix (const int Polynom1, const int Polynom2,int **inputmat, int **nextstate, int **outmat)

{

	int cpt_line, cpt_column,cpt;
	int next_decimal,out;
	//int current_state[Length_Poly-1],reg[Length_Poly];
	int*current_state,*reg;
	double out1,out2;
    int Length_Poly=3;

	current_state = (int*)malloc((Length_Poly-1)*sizeof(int));
	reg = (int*)malloc((Length_Poly)*sizeof(int));


	out1=0.;
	out2=0.;
	next_decimal=0;

	// initilisation of the transition matrix
	// 2 for impossible transition

	for ( cpt_line=0; cpt_line<(int)(pow(2,Length_Poly-1)); cpt_line++)
	{
		for ( cpt_column=0; cpt_column<(int)(pow(2,Length_Poly-1)); cpt_column++)
		{
			inputmat[cpt_line][cpt_column]=2;
		}
	}

	// compute output matrix and nextstate matrix

	for ( cpt_line=0; cpt_line<(int)(pow(2,Length_Poly-1)); cpt_line++)
	{
		for ( cpt_column=0; cpt_column<2; cpt_column++)
		{
			//compute the actual state in binary
			int2bin (cpt_line,current_state,Length_Poly-1);
			//compute the next state in binary
			shift_reg (current_state,Length_Poly-1, cpt_column);
			//compute the next state in integer
			bin2int (&next_decimal,current_state,Length_Poly-1);
			// stockage
			inputmat[cpt_line][next_decimal]=cpt_column;
			nextstate[cpt_line][cpt_column]=next_decimal;
			//compute the actual state in binary
			int2bin (cpt_line,current_state,Length_Poly-1);

			// register initialisation for the computation of the outputs encoder
			for (cpt=1; cpt<Length_Poly; cpt++)
			{
				reg[cpt]=current_state[cpt-1];
			}
			reg[0]=cpt_column;
			//computation of the outputs encoder
			encode_bit (Polynom1,Polynom2,reg,Length_Poly,&out1,&out2);

			// transform the outputs in decimal
			if ((out1==1) && (out2==1)) {out=0;}
			if ((out1==-1) && (out2==1)) {out=2;}
			if ((out1==1) && (out2==-1)) {out=1;}
			if ((out1==-1) && (out2==-1)) {out=3;}
            outmat[cpt_line][cpt_column]=out;
		}
	}
	free(current_state);
	free(reg);
	return OK;
}




int va_detection(int MemoryDepth, int Nbstates, long double *target, int target_length, long double *input, int *output,int ip_length)
{


     int CS[4][2]={ { 1,1 }, { 2, 2 }, { 3,3 },{ 4, 4}};
     int NS[4][2]={{0,1},{2,3},{0,1},{2,3}};
     int IP[4][4]={{-1,1,0,0},{0,0,-1,1},{-1,1,0,0},{0,0,-1,1}};
     int OP[4][2]={{0,1},{1,0},{1,0},{0,1}};
     int state_values[4][2]={{-1,-1},{-1,1},{1,-1},{1,1}};

	register int i,j,k;

	int error_temp,buffer_position,min_metric_val,min_metric_state;
	int **state_history,*state_sequence;

	long double x_temp,x,max_weight,branch_metric;
	long double **error_accumulation;

    error_accumulation=matrix_calloc(Nbstates,2);
	state_sequence = vector_calloc_int(MemoryDepth+1);
    state_history=matrix_calloc_int(Nbstates,MemoryDepth+1);

	/****************************************************************************/
	/*								Initialisation		 						*/
	/****************************************************************************/
	// give a huge value for the maximum metric
	max_weight=(long double)(ip_length*10000);
	// initialisation of the metrics
	for (i=0; i<Nbstates; i++)
	{
		error_accumulation[i][0]=0;
		error_accumulation[i][1]=max_weight;
	} // end of  for (i=0; i<NBSTATE; i++)

	/****************************************************************************/
	/*							Decode the data									*/
	/****************************************************************************/
	for (k=0; k<ip_length; k++)
	{
		// compute the circulary column indice
		buffer_position = ( (k + 1) % ( MemoryDepth + 1 ));
		// compute the metric for each state
		for ( i=0; i<Nbstates; i++)
		{
			// compute the metric for each transition
			for ( j=0; j<2; j++)
			{
				branch_metric = 0.;
				if (j==0)
                {
                        branch_metric=(-1*target[0])+(target[1]*state_values[i][1])+(target[2]*state_values[i][0]);
                }
                else
                {
                        branch_metric=(1*target[0])+(target[1]*state_values[i][1])+(target[2]*state_values[i][0]);
                }
				// give the output encoder in function of the state and the transition
            	branch_metric=pow((input[k]-branch_metric),2);
							// choose the survivor path
				if (error_accumulation[NS[i][j]][1]>(error_accumulation[i][0]+branch_metric))

				{
					error_accumulation[NS[i][j]][1]=(error_accumulation[i][0]+branch_metric);
					state_history[NS[i][j]][buffer_position]=i;
				}

			}// end of  for ( j=0; j<2; j++)

		}
		// update the trellis metric
		for (i=0; i < Nbstates ; i++)
		{
			error_accumulation[i][0]=error_accumulation[i][1];
			error_accumulation[i][1]=max_weight;
		}
		// Reconstrcution of the data when the state_history is full
		if (k>=MemoryDepth-1)
		{
			// initilisation of state_sequence !!! need not need ???
			var_init(state_sequence,MemoryDepth);
			// find the minimum metric
			// up and down searh
			x=max_weight;
			for (i=0; i <(Nbstates/2); i++ )
			{
				if (error_accumulation[i][0]<error_accumulation[Nbstates-1-i][0])
				{
					x_temp=error_accumulation[i][0];
					error_temp=i;
				}
				else
				{
					x_temp=error_accumulation[Nbstates-1-i][0];
					error_temp=Nbstates-1-i;

				}
				if (x_temp<x)
				{
					min_metric_state=error_temp;
					x=x_temp;
				} // if (x_temp<x)
			} // end of  for (i=0; i <(NBSTATE/2); i++ )
			// initialisation of the start point for the data reconstruction
			state_sequence[MemoryDepth]=min_metric_state;
			// constitution of the path with the minimum metric
			for (i=MemoryDepth; i > 0 ; i--)
			{
				min_metric_val= i + buffer_position - ( MemoryDepth /*+ 1*/ );
				if ( min_metric_val < 0 )
				{
					min_metric_val+=(MemoryDepth+1);
				} // end of  if ( min_metric_val < 0 )
				state_sequence[i-1]=state_history[state_sequence[i]][min_metric_val];
			} // end of  for (i= MEMORY_DEPTH; i < 1; i--);
			output[k-MemoryDepth+1]=IP[state_sequence[0]][state_sequence[1]];
		} // end of  if (k>=MEMORY_DEPTH)
	} // end of  for (k=0; k<nb_bits; k++)
	/****************************************************************************/
	/*						empty the state_sequence							*/
	/****************************************************************************/
	for (i=1; i < (MemoryDepth ); i++)
		{
		output[ip_length-MemoryDepth+i]=IP[state_sequence[i]][state_sequence[i+1]];
		} // end of  for (i=1; i < MEMORY_DEPTH-LENGTHPOLY; i++)
	x=max_weight;
	x_temp=0;
	// compute the survivor path metric and the mean of the non survivor
	for (i=0; i<Nbstates; i++)
	{
		x_temp+=error_accumulation[i][0];
		if (error_accumulation[i][0]<x) { x=error_accumulation[i][0];}
	}
	x_temp-=x;
	free(state_sequence);
	matrix_free(error_accumulation,Nbstates);
	matrix_free(state_history,Nbstates);
/*	for(i =0; i<NbStates; i++)
		free(state_history[i]);
	free(state_history);
*/
	return OK;
}

int sova_detection(int MemoryDepth, int Nbstates, long double *target, int target_length, long double *input, int *output,long double *llr,long int *llr_length,int ip_length)
{
    //Memory depth is 5*contraint length
    //SOVA depth is 2*contraint length

    //Here we do an overkill to avoid any limitation i.e. depth >10 times constraint length

     int CS[4][2]={ { 1,1 }, { 2, 2 }, { 3,3 },{ 4, 4}};
     int NS[4][2]={{0,1},{2,3},{0,1},{2,3}};
     int IP[4][4]={{-1,1,0,0},{0,0,-1,1},{-1,1,0,0},{0,0,-1,1}};
     int OP[4][2]={{0,1},{1,0},{1,0},{0,1}};
     int comp_state[4][2]={{2,2},{3,3},{0,0},{1,1}}; //state zero corresponds to the "1" state matlab

     int state_values[4][2]={{-1,-1},{-1,1},{1,-1},{1,1}};

	register int i,j,k;

	int bit,corr_col,l,reqd_col,m;

	int error_temp,buffer_position,min_metric_val,min_metric_state;
	int **state_history,**complement_state_history,*state_sequence,*comp_state_sequence;

	long double x_temp,x,max_weight,branch_metric;
	long double **error_accumulation,**metr_diff;

    error_accumulation=matrix_calloc(Nbstates,2);
    metr_diff=matrix_calloc(Nbstates,MemoryDepth+1);

	state_sequence = vector_calloc_int(MemoryDepth+1);
	comp_state_sequence = vector_calloc_int(MemoryDepth+1);

    state_history=matrix_calloc_int(Nbstates,MemoryDepth+1);
    complement_state_history=matrix_calloc_int(Nbstates,MemoryDepth+1);

	/****************************************************************************/
	/*								Initialisation		 						*/
	/****************************************************************************/
	// give a huge value for the maximum metric
	max_weight=(long double)(ip_length*10000);
	// initialisation of the metrics
	for (i=0; i<Nbstates; i++)
	{
		error_accumulation[i][0]=0;
		error_accumulation[i][1]=max_weight;
	} // end of  for (i=0; i<NBSTATE; i++)

	/****************************************************************************/
	/*							Decode the data									*/
	/****************************************************************************/
	for (k=0; k<ip_length; k++)
	{
		// compute the circulary column indice
		buffer_position = ( (k) % ( MemoryDepth ));
		// compute the metric for each state
		for ( i=0; i<Nbstates; i++)
		{
			// compute the metric for each transition
			for ( j=0; j<2; j++)
			{
				branch_metric = 0.;
				if (j==0)
                {
                        branch_metric=(-1*target[0])+(target[1]*state_values[i][1])+(target[2]*state_values[i][0]);
                }
                else
                {
                        branch_metric=(1*target[0])+(target[1]*state_values[i][1])+(target[2]*state_values[i][0]);
                }
				// give the output encoder in function of the state and the transition
            	branch_metric=pow((input[k]-branch_metric),2);
            	//Store metric difference once both path branch metric are computed
            	if (error_accumulation[NS[i][j]][1]!=max_weight)
            	{
                    metr_diff[NS[i][j]][buffer_position]=(error_accumulation[NS[i][j]][1]-(error_accumulation[i][0] + branch_metric));
                    if (metr_diff[NS[i][j]][buffer_position] < 0)
                        {
                            metr_diff[NS[i][j]][buffer_position]=-metr_diff[NS[i][j]][buffer_position];
                        }
            	}

                // choose the survivor path
				if (error_accumulation[NS[i][j]][1]>(error_accumulation[i][0]+branch_metric))
				{
					error_accumulation[NS[i][j]][1]=(error_accumulation[i][0]+branch_metric);
					state_history[NS[i][j]][buffer_position]=i;
					complement_state_history[NS[i][j]][buffer_position]=comp_state[i][j];  //stores complement state for each possible state
				}

			}// end of  for ( j=0; j<2; j++)

		}
		// update the trellis metric
		for (i=0; i < Nbstates ; i++)
		{
			error_accumulation[i][0]=error_accumulation[i][1];
			error_accumulation[i][1]=max_weight;
		}
		// Reconstrcution of the data when the state_history is full
		if (k>=MemoryDepth-1)
		{
			// initilisation of state_sequence !!! need not need ???
			var_init(state_sequence,MemoryDepth);
			// find the minimum metric
			// up and down searh
			x=max_weight;
			for (i=0; i <(Nbstates/2); i++ )
			{
				if (error_accumulation[i][0]<error_accumulation[Nbstates-1-i][0])
				{
					x_temp=error_accumulation[i][0];
					error_temp=i;
				}
				else
				{
					x_temp=error_accumulation[Nbstates-1-i][0];
					error_temp=Nbstates-1-i;

				}
				if (x_temp<x)
				{
					min_metric_state=error_temp;
					x=x_temp;
				} // if (x_temp<x)
			} // end of  for (i=0; i <(NBSTATE/2); i++ )
			// initialisation of the start point for the data reconstruction
			state_sequence[MemoryDepth-1]=state_history[min_metric_state][buffer_position];

			//complement start point
			comp_state_sequence[MemoryDepth-1]=complement_state_history[min_metric_state][buffer_position];

			// constitution of the path with the minimum metric
			for (i=MemoryDepth-1-1; i >= 0 ; i--)
			{
				min_metric_val= i + buffer_position - ( MemoryDepth-1 /*+ 1*/ );
				if ( min_metric_val < 0 )
				{
					min_metric_val+=(MemoryDepth);
				} // end of  if ( min_metric_val < 0 )
				state_sequence[i]=state_history[state_sequence[i+1]][min_metric_val];
				comp_state_sequence[i]=state_history[comp_state_sequence[i+1]][min_metric_val];
			} // end of  for (i= MEMORY_DEPTH; i < 1; i--);
			output[k-MemoryDepth+1]=IP[state_sequence[0]][state_sequence[1]];

			//SOVA traceback and llr computation
            //state history works on a 1 to memory depth basis
			llr[k-MemoryDepth+1]=max_weight;
            bit=IP[comp_state_sequence[0]][comp_state_sequence[1]];

            if (output[k-MemoryDepth+1]!=bit)
                llr[k-MemoryDepth+1]=min(llr[k-MemoryDepth+1], metr_diff[min_metric_state][buffer_position]);

            for (i=1;i<=MemoryDepth-1-1;i++) //mem_depth-sova_depth and 1 less transition as already covered above
            //Takes the next survivor path point and takes the complement path for this point.
            {
                corr_col=buffer_position-i;

                if (corr_col < 0)
                {
                     corr_col=corr_col+(MemoryDepth);
                }

                l=MemoryDepth-i-1; //-1 for C coding
                comp_state_sequence[l]=complement_state_history[state_sequence[l+1]][corr_col];

                for (j=1;j<=l;j++) //corresponds to 1 less transition to trace back
                {
                        m=l-j;
                    reqd_col=corr_col - j;
                    if (reqd_col < 0)
                    {
                    reqd_col=reqd_col+MemoryDepth;
                    }
                comp_state_sequence[m]=state_history[comp_state_sequence[m+1]][reqd_col];
                }

                bit=IP[comp_state_sequence[0]][comp_state_sequence[1]];
                if (output[k-MemoryDepth+1]!=bit)
                {
                    llr[k-MemoryDepth+1]=min(llr[k-MemoryDepth+1], metr_diff[state_sequence[l+1]][corr_col]);
                }
            }

            llr[k-MemoryDepth+1]=llr[k-MemoryDepth+1]*output[k-MemoryDepth+1];

		} // end of  if (k>=MEMORY_DEPTH)
	} // end of  for (k=0; k<nb_bits; k++)

	*llr_length=k-MemoryDepth+1+1;   //+1 as net samples is from 1->k

	//Last bits in the buffer, the llr is not RELIABLE

	/****************************************************************************/
	/*						empty the state_sequence							*/
	/****************************************************************************/
	for (i=1; i < (MemoryDepth ); i++)
		{
		output[ip_length-MemoryDepth+i]=IP[state_sequence[i]][state_sequence[i+1]];
		} // end of  for (i=1; i < MEMORY_DEPTH-LENGTHPOLY; i++)
	x=max_weight;
	x_temp=0;
	// compute the survivor path metric and the mean of the non survivor
	for (i=0; i<Nbstates; i++)
	{
		x_temp+=error_accumulation[i][0];
		if (error_accumulation[i][0]<x) { x=error_accumulation[i][0];}
	}
	x_temp-=x;
	free(state_sequence);
	matrix_free(error_accumulation,Nbstates);
	matrix_free(state_history,Nbstates);
/*	for(i =0; i<NbStates; i++)
		free(state_history[i]);
	free(state_history);
*/
	return OK;
}



//int old_map_detection(int MemoryDepth, int Nbstates, long double *target, int target_length, long double *input, int *output,long double *llr,long int *llr_length,int ip_length,float SNR)
//{
//
//    float es_over_n0,es,sn_ratio;
//    long double sigma,noise_sigma,variance,channel_reliability,No;
//
//
//    int K,m,states;
//    //MAP compute sigma for branch metric computation
//    es_over_n0=SNR;  //Asuming the input is the Es_over_N0
//    es = 1;  //iess normalized to  1
//    sn_ratio =  pow(10,( es_over_n0 / 10));
//    No =  es/(sn_ratio);
//
//    //variance=pow(noise_sigma,2);
//
//    //channel_reliability=variance/2;
//
//    //input_length corresponds to block_s in matlab
//    K=3; //Constraint Length
//    m=K-1;
//    states=pow(2,m); //Number of states
//
//
//
//    int op_level=2;
//
//
//    //Here we do an overkill to avoid any limitation i.e. depth >10 times constraint length
//
//     int CS[4][2]={ { 1,1 }, { 2, 2 }, { 3,3 },{ 4, 4}};
//     int NS[4][2]={{0,1},{2,3},{0,1},{2,3}};
//     int PS[4][2]={{0,2},{2,0},{3,1},{1,3}};
//     int OP[4][2]={{-1,1},{1,-1},{1,-1},{-1,1}};
//
//     int comp_state[4][2]={{2,2},{3,3},{0,0},{1,1}}; //state zero corresponds to the "1" state matlab
//
//     int state_values[4][2]={{-1,-1},{-1,1},{1,-1},{1,1}};
//
//     long double **prob,***gamma, **alpha , **beta, *A_norm, *B_norm;
//
//
//    //Probability computation
//    prob=matrix_calloc(ip_length,2);
//
//    //forward recursion init
//    gamma=matrix_calloc(MemoryDepth,states,op_level);
//    alpha=matrix_calloc(MemoryDepth,states);
//    A_norm=vector_calloc(MemoryDepth);
//    //backward recursion
//    //gamma_b=zeros(mem_depth,no_of_states,op_level);
//
//
//
//    register int h,i,j,k;
//    int depth_history,depth_history_beta,prev_val,col_val,next_val,bit_col,temp;
//    long double branch_metric,Nx;
//
//    for (i=0;i<ip_length;i++)
//    {
//    /****************************************************************************/
//	/*						Forward Recursion       							*/
//	/****************************************************************************/
//         //depth is from 0- Memory Depth
//         depth_history= (i)%(MemoryDepth);
//         //initialize alpha
//        for (j=0;j<states;j++)
//        {
//           alpha[depth_history][j]=0;
//        }
//        A_norm[depth_history]=0;
//        //forward recursion
//        for (j=0;j<states;j++)
//        {
//            for (k=0;k<op_level;k++)
//                {
//                  if (i==0)
//                  {
//                    alpha[0][0]=1;
//                  }
//                  else
//                  {
//                    if (depth_history-1<0)
//                       prev_val=MemoryDepth-1;
//                   else
//                       prev_val=depth_history-1;
//
//                  alpha[depth_history][j]+=alpha[prev_val][PS[j][k]]*gamma[prev_val][PS[j][k]][k];
//                  }
//                  branch_metric = 0.;
//                  if (k==0)
//                  {
//                    branch_metric=(-1*target[0])+(target[1]*state_values[j][1])+(target[2]*state_values[j][0]);
//                  }
//                  else
//                  {
//                        branch_metric=(1*target[0])+(target[1]*state_values[j][1])+(target[2]*state_values[j][0]);
//                  }
//                  // give the output encoder in function of the state and the transition
//
//                  branch_metric=pow((input[i]-branch_metric),2);
//                  //gamma[depth_history][j][k]=exp(-(1/(2*variance))*branch_metric);
//                  gamma[depth_history][j][k]=exp(-(1/No)*branch_metric);
//                }
//             A_norm[depth_history]+=alpha[depth_history][j];
//        }
//
//       for (j=0;j<states;j++)
//           alpha[depth_history][j]/=A_norm[depth_history];
//
//    if (i>=MemoryDepth-1)
//    {
//        //depth_history_beta=(depth_history+1)%MemoryDepth;
//        //backward recursion Initialize
//        beta=matrix_calloc(MemoryDepth,states);
//        //initialize
//        //for (j=0;j<states;j++)
//        //    beta[depth_history][j]=(long double)(1/(long double)states);
//
//        B_norm=vector_calloc(MemoryDepth);
//        //Initialize beta with alpha
//        for (j=0;j<states;j++)
//            beta[depth_history][j]=alpha[depth_history][j];
//        for (h=0;h<MemoryDepth-1;h++)
//        {
//            col_val=depth_history-h-1;
//            if(col_val<0)
//                col_val=col_val+MemoryDepth;
//            next_val=col_val+1;
//            if (next_val>MemoryDepth-1)
//                next_val=next_val-MemoryDepth;
//
//            for (j=0;j<states;j++)
//            {
//                for (k=0;k<op_level;k++)
//                {
//                  beta[col_val][j]+=beta[next_val][NS[j][k]]*gamma[next_val][j][k];
//                }
//                B_norm[col_val]+=beta[col_val][j];
//            }
//            for (j=0;j<states;j++)
//               {
//                 beta[col_val][j]/=B_norm[col_val];
//               }
//
//        }
//
//        bit_col=(depth_history)+1;
//        if (bit_col>MemoryDepth-1)
//            bit_col=0;
//        for (j=0;j<states;j++)
//        {
//             prob[i-MemoryDepth+1][0]+=alpha[bit_col][j]*gamma[bit_col][j][0]*beta[bit_col][NS[j][0]];
//             prob[i-MemoryDepth+1][1]+=alpha[bit_col][j]*gamma[bit_col][j][1]*beta[bit_col][NS[j][1]];
//        }
//
//        Nx=prob[i-MemoryDepth+1][0]+prob[i-MemoryDepth+1][1];
//
//        //Normalize Probability
//        prob[i-MemoryDepth+1][0]=prob[i-MemoryDepth+1][0]/Nx;
//        prob[i-MemoryDepth+1][1]=prob[i-MemoryDepth+1][1]/Nx;
//
//
//
//        if((prob[i-MemoryDepth+1][0])>0.5)
//        {
//             output[i-MemoryDepth+1]=-1;
//        }
//        else
//        {
//            output[i-MemoryDepth+1]=1;
//        }
//
//        if (i!=ip_length-1)         //last set of betas required to empty out
//        {
//            matrix_free(beta,MemoryDepth);
//        }
//
//
//        free(B_norm);
//
//    }
//    }
//    *llr_length=i-MemoryDepth+1+1;
//
//    for (i=1;i<MemoryDepth-1;i++)
//    {
//        bit_col=bit_col+1;
//        if (bit_col>MemoryDepth-1)
//            bit_col=0;
//
//
//        for (j=0;j<states;j++)
//            {
//                prob[ip_length-MemoryDepth+i][0]+=(alpha[bit_col][j]*gamma[bit_col][j][0]*beta[bit_col][NS[j][0]]);
//                prob[ip_length-MemoryDepth+i][1]+=(alpha[bit_col][j]*gamma[bit_col][j][1]*beta[bit_col][NS[j][1]]);
//            }
//
//
//
//        Nx=prob[ip_length-MemoryDepth+i][0]+prob[ip_length-MemoryDepth+i][1];
//        prob[ip_length-MemoryDepth+i][0]=prob[ip_length-MemoryDepth+i][0]/Nx;
//        prob[ip_length-MemoryDepth+i][1]=prob[ip_length-MemoryDepth+i][1]/Nx;
//    if((prob[ip_length-MemoryDepth+i][0])>0.5)
//        output[ip_length-MemoryDepth+i]=-1;
//    else
//        output[ip_length-MemoryDepth+i]=1;
//
//
//    }
//
//
//    matrix_free(prob,ip_length);
//    matrix_free(gamma,MemoryDepth,states);
//    matrix_free(alpha,MemoryDepth);
//    free(A_norm);
///*	for(i =0; i<NbStates; i++)
//		free(state_history[i]);
//	free(state_history);
//*/
//	return OK;
//}



int map_detection(int MemoryDepth, int Nbstates, long double *target, int target_length, long double *input, int *output,long double *llr,long int *llr_length,int ip_length,float SNR)
{

    float es_over_n0,es,sn_ratio;
    long double sigma,noise_sigma,variance,channel_reliability,No;


    int K,m,states;
    //MAP compute sigma for branch metric computation
    es_over_n0=SNR;  //Asuming the input is the Es_over_N0
    es = 1;  //iess normalized to  1
    sn_ratio =  pow(10,( es_over_n0 / 10));
    sigma = sqrt (es/sn_ratio);
    sigma = sigma/sqrt(2);
    variance=pow(sigma,2);

    //channel_reliability=variance/2;

    //input_length corresponds to block_s in matlab
    K=3; //Constraint Length
    m=K-1;
    states=pow(2,m); //Number of states



    int op_level=2;


    //Here we do an overkill to avoid any limitation i.e. depth >10 times constraint length

     int CS[4][2]={ { 0,0 }, { 1, 1 }, { 2,2 },{ 3, 3}};
     int NS[4][2]={{0,1},{2,3},{0,1},{2,3}};
     int PS[4][2]={{0,2},{2,0},{3,1},{1,3}};
     int OP[4][2]={{-1,1},{-1,1},{-1,1},{-1,1}};

     int comp_state[4][2]={{2,2},{3,3},{0,0},{1,1}}; //state zero corresponds to the "1" state matlab

     int state_values[4][2]={{-1,-1},{-1,1},{1,-1},{1,1}};

     long double **prob,***gamma, **alpha , **beta, *A_norm, *B_norm;
     register int h,i,j,k;

    //Probability computation
    prob=matrix_calloc(ip_length,2);

    //forward recursion init
    gamma=matrix_calloc(MemoryDepth,states,op_level);
    alpha=matrix_calloc(MemoryDepth+1,states);
    A_norm=vector_calloc(MemoryDepth+1);
    //backward recursion
    //gamma_partial is not required as apriori probability is "0"


    //Initialize alpha
    alpha[0][0]=0;
    for (j=1;j<states;j++)
    {
      alpha[0][j]=-1.5e-200;
    }


    int depth_history,alpha_depth_history,next_alpha,prev_val,col_val,next_val,bit_col,alpha_bit_col,temp;
    long double branch_metric,Nx;

    for (i=0;i<ip_length;i++)
    {
    /****************************************************************************/
	/*						Forward Recursion       							*/
	/****************************************************************************/
         //depth is from 0- Memory Depth
         depth_history= (i)%(MemoryDepth);
         alpha_depth_history= (i)%(MemoryDepth+1);



         if((alpha_depth_history+1)==(MemoryDepth+1))
        {
             next_alpha=0;
        }
        else
         {
             next_alpha=alpha_depth_history+1;
         }

         //initialize alpha
        for (j=0;j<states;j++)
        {
           alpha[next_alpha][j]=0;
        }
        A_norm[next_alpha]=0;


        //forward recursion
        for (j=0;j<states;j++)
        {
            for (k=0;k<op_level;k++)
                {
                  branch_metric = 0.;
                  if (k==0)
                  {
                    branch_metric=(-1*target[0])+(target[1]*state_values[j][1])+(target[2]*state_values[j][0]);
                  }
                  else
                  {
                    branch_metric=(1*target[0])+(target[1]*state_values[j][1])+(target[2]*state_values[j][0]);
                  }
                  // give the output encoder in function of the state and the transition

                  branch_metric=pow((input[i]-branch_metric),2);
                  //gamma[depth_history][j][k]=exp(-(1/(2*variance))*branch_metric);
                  gamma[depth_history][j][k]=exp(-(1/variance)*branch_metric);

                  alpha[next_alpha][NS[j][k]]+=alpha[alpha_depth_history][j]*gamma[depth_history][j][k];
                }

        }
       for (j=0;j<states;j++)
            A_norm[next_alpha]+=alpha[next_alpha][j];
       for (j=0;j<states;j++)
           alpha[next_alpha][j]/=A_norm[next_alpha];

    if (i>=MemoryDepth-1)
    {
        //backward recursion Initialize
        beta=matrix_calloc(MemoryDepth,states);
        //initialize
        for (j=0;j<states;j++)
            beta[depth_history][j]=(long double)(1/(long double)states);

        B_norm=vector_calloc(MemoryDepth);
        //Initialize beta with alpha
        //for (j=0;j<states;j++)
        //    beta[depth_history][j]=alpha[depth_history][j];

        for (h=0;h<MemoryDepth-1;h++)
        {
            col_val=depth_history-h-1;
            if(col_val<0)
                col_val=col_val+MemoryDepth;
            next_val=col_val+1;
            if (next_val>MemoryDepth-1)
                next_val=next_val-MemoryDepth;

            for (j=0;j<states;j++)
            {
                for (k=0;k<op_level;k++)
                {
                  beta[col_val][j]+=beta[next_val][NS[j][k]]*gamma[next_val][j][k];
                }
             B_norm[col_val]+=beta[col_val][j];
            }
            for (j=0;j<states;j++)
               {
                 beta[col_val][j]/=B_norm[col_val];
               }

        }

        bit_col=(depth_history)+1;
        if (bit_col>MemoryDepth-1)
            bit_col=0;

        alpha_bit_col=next_alpha+1;

        if(alpha_bit_col>MemoryDepth)
            alpha_bit_col=0;



        for (j=0;j<states;j++)
        {
             prob[i-MemoryDepth+1][0]+=alpha[alpha_bit_col][j]*gamma[bit_col][j][0]*beta[bit_col][NS[j][0]];
             prob[i-MemoryDepth+1][1]+=alpha[alpha_bit_col][j]*gamma[bit_col][j][1]*beta[bit_col][NS[j][1]];
        }

        Nx=prob[i-MemoryDepth+1][0]+prob[i-MemoryDepth+1][1];

        //Normalize Probability
        prob[i-MemoryDepth+1][0]=prob[i-MemoryDepth+1][0]/Nx;
        prob[i-MemoryDepth+1][1]=prob[i-MemoryDepth+1][1]/Nx;

        if (prob[i-MemoryDepth+1][0]==0)
            prob[i-MemoryDepth+1][0]=1e-200;
        if (prob[i-MemoryDepth+1][1]==0)
            prob[i-MemoryDepth+1][1]=1e-200;


        llr[i-MemoryDepth+1]=prob[i-MemoryDepth+1][1]/prob[i-MemoryDepth+1][0];

        llr[i-MemoryDepth+1]=log(llr[i-MemoryDepth+1]);

        if((prob[i-MemoryDepth+1][0])>0.5)
        {
             output[i-MemoryDepth+1]=-1;
        }
        else
        {
            output[i-MemoryDepth+1]=1;
        }

        if (i!=ip_length-1)         //last set of betas required to empty out
        {
            matrix_free(beta,MemoryDepth);
        }


        free(B_norm);

    }
    }
    *llr_length=i-MemoryDepth+1+1;          //Do not consider last bits as they are not as accurate as the other bits
                                            //through the window and could potentially increase ber

    for (i=1;i<MemoryDepth-1;i++)
    {
        bit_col=bit_col+1;
        if (bit_col>MemoryDepth-1)
            bit_col=0;

        alpha_bit_col=next_alpha+1;

        if(alpha_bit_col>MemoryDepth)
            alpha_bit_col=0;



        for (j=0;j<states;j++)
            {
                prob[ip_length-MemoryDepth+i][0]+=(alpha[alpha_bit_col][j]*gamma[bit_col][j][0]*beta[bit_col][NS[j][0]]);
                prob[ip_length-MemoryDepth+i][1]+=(alpha[alpha_bit_col][j]*gamma[bit_col][j][1]*beta[bit_col][NS[j][1]]);
            }



        Nx=prob[ip_length-MemoryDepth+i][0]+prob[ip_length-MemoryDepth+i][1];
        prob[ip_length-MemoryDepth+i][0]=prob[ip_length-MemoryDepth+i][0]/Nx;
        prob[ip_length-MemoryDepth+i][1]=prob[ip_length-MemoryDepth+i][1]/Nx;

    if((prob[ip_length-MemoryDepth+i][0])>0.5)
        output[ip_length-MemoryDepth+i]=-1;
    else
        output[ip_length-MemoryDepth+i]=1;


    }


    matrix_free(prob,ip_length);
    matrix_free(gamma,MemoryDepth,states);
    matrix_free(alpha,MemoryDepth);
    free(A_norm);
/*	for(i =0; i<NbStates; i++)
		free(state_history[i]);
	free(state_history);
*/
	return OK;
}


long double ber_compute(int *input1,int *input2,int Nb_values)
{

    register int i;http://www.nytimes.com/
    long double ber_count;

    ber_count=0;
    //1 bit is swallowed by viterbi!!:-(
    for (i=0;i<Nb_values-1;i++)
    {
            if (input1[i]!=input2[i])
            {
                    ber_count=ber_count+1;
            }
    }
    ber_count=ber_count/Nb_values;
 return ber_count;
}





long double oneD_equalize_oneD_SOVA(int *tx_bits,int current_sample_length, long double *rx_bits,int Nb_eq,int Nb_target,int *ber_bits,float SNR)
{


  long double *equalizer,*target;

  int track_count;

  long double *equalized_op;

  long double ber_val;

    int MemoryDepth;
    int NbStates;
    int *SOVA_out;

    long double *llr;
    long int *llr_length;

    /**********Equalizer computation*********/
    track_count=1;          // Assuming 1 track on each side influence
    /*************************************/

    /**********SOVA Parameters*********/
    //Currently har coded for 4 states
    //MemoryDepth=200;
    MemoryDepth=40;
    NbStates=4;
    /*************************************/


/********************EQUALIZER*************************/

    equalizer= vector_calloc(Nb_eq);
    target= vector_calloc(Nb_target);

//    //Compute Equalizer and Target Coeff
    lagrange_err_min_1D(equalizer,Nb_eq,target,Nb_target,tx_bits,current_sample_length,rx_bits);

    //Equalize
    equalized_op=vector_calloc(current_sample_length+Nb_eq);
    convolution(rx_bits,equalizer,current_sample_length,Nb_eq,equalized_op);

    //truncate initial samples and end values from convolution process
    array_resize(equalized_op,Nb_eq-1,current_sample_length);

    //Length of sample is sample_length-(Filter_coeff-1);
    current_sample_length=current_sample_length-(Nb_eq-1);

//Under investigation!!! how to generate states for variable target size
//int MemoryDepth=200, NbStates=4;
    /*int Polynom1=3;
    int Polynom2=1;
    int Nb_states;
    int **inputmat,**nextstate,**outmat;

    Nb_states=pow(2,no_of_target_vals-1);

    inputmat=matrix_calloc_int(Nb_states,Nb_states);
    nextstate=matrix_calloc_int(Nb_states,2);
    outmat=matrix_calloc_int(Nb_states,2);

    gene_viterbi_matrix (Polynom1,Polynom2,inputmat,nextstate,outmat);
*/

    SOVA_out=vector_calloc_int(current_sample_length);


    llr=vector_calloc(current_sample_length);

    /************Include for VA***********************/
    va_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,current_sample_length);
    /************Include for VA***********************/

    llr_length=(long int *)malloc( sizeof(long int));

    /************Include for SOVA***********************/
    //sova_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,llr,llr_length,current_sample_length);

    //current_sample_length=*llr_length;

    /************Include for SOVA***********************/


    /************Include for MAP***********************/
    //map_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,llr,llr_length,current_sample_length,SNR);

    //current_sample_length=*llr_length;

    /************Include for MAP***********************/



    ber_val=ber_compute(ber_bits,SOVA_out,current_sample_length);


   free (equalized_op);

   free (llr_length);
   free (llr);

   free (equalizer);
   free (target);

    free (SOVA_out);

    return ber_val;
}

long double twoD_equalize_oneD_SOVA(int **tx_merge,int current_sample_length, long double *rx_main,long double *rx_adj1_OD,long double *rx_adj1_ID,int Nb_eq_main,int Nb_eq_adj1_OD,int Nb_eq_adj1_ID,int Nb_target,int *ber_bits,float SNR)
{

  long double *equalizer_main,*equalizer_adj1_OD,*equalizer_adj1_ID,*target;

  int track_count;

  long double *equalized_main,*equalized_adj1_OD,*equalized_adj1_ID,*equalized_op;

  long double ber_val;

    int MemoryDepth;
    int NbStates;
    int *SOVA_out;

    long double *llr;
    long int *llr_length;

    /**********Equalizer computation*********/
    track_count=1;          // Assuming 1 track on each side influence
    /*************************************/

    /**********SOVA Parameters*********/
    //Currently har coded for 4 states
    MemoryDepth=40;
    NbStates=4;
    /*************************************/


/********************EQUALIZER*************************/

    equalizer_main= vector_calloc(Nb_eq_main);
    equalizer_adj1_OD= vector_calloc(Nb_eq_adj1_OD);
    equalizer_adj1_ID= vector_calloc(Nb_eq_adj1_ID);

    target= vector_calloc(Nb_target);

    //Compute Equalizer and Target Coeff
    lagrange_err_min_2D(equalizer_main,equalizer_adj1_OD,equalizer_adj1_ID,Nb_eq_main,Nb_eq_adj1_OD,Nb_eq_adj1_ID,target,Nb_target,tx_merge,current_sample_length,rx_main,rx_adj1_OD,rx_adj1_ID);

    //Equalize Main
    equalized_main=vector_calloc(current_sample_length+Nb_eq_main);
    convolution(rx_main,equalizer_main,current_sample_length,Nb_eq_main,equalized_main);

    //Equalize Adj1 OD
    equalized_adj1_OD=vector_calloc(current_sample_length+Nb_eq_adj1_OD);
    convolution(rx_adj1_OD,equalizer_adj1_OD,current_sample_length,Nb_eq_adj1_OD,equalized_adj1_OD);

    //Equalize Adj1 ID
    equalized_adj1_ID=vector_calloc(current_sample_length+Nb_eq_adj1_ID);
    convolution(rx_adj1_ID,equalizer_adj1_ID,current_sample_length,Nb_eq_adj1_ID,equalized_adj1_ID);


    //truncate initial samples and end values from convolution process
    array_resize(equalized_main,Nb_eq_main-1,current_sample_length);

    //truncate initial samples and end values from convolution process
    array_resize(equalized_adj1_OD,Nb_eq_adj1_OD-1,current_sample_length);

//truncate initial samples and end values from convolution process
    array_resize(equalized_adj1_ID,Nb_eq_adj1_ID-1,current_sample_length);

    //!!!!!!!!!!!!!!!!!Have to see how to handle variable equalizer length??!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Length of sample is sample_length-(Filter_coeff-1);
    current_sample_length=current_sample_length-(Nb_eq_main-1);



//Under investigation!!! how to generate states for variable target size
//int MemoryDepth=200, NbStates=4;
    /*int Polynom1=3;
    int Polynom2=1;
    int Nb_states;
    int **inputmat,**nextstate,**outmat;

    Nb_states=pow(2,no_of_target_vals-1);

    inputmat=matrix_calloc_int(Nb_states,Nb_states);
    nextstate=matrix_calloc_int(Nb_states,2);
    outmat=matrix_calloc_int(Nb_states,2);

    gene_viterbi_matrix (Polynom1,Polynom2,inputmat,nextstate,outmat);
*/

    //Add Equalized outputs
    equalized_op=vector_calloc(current_sample_length);

    array_addition (equalized_op,equalized_main,equalized_adj1_OD,current_sample_length);
    array_addition (equalized_op,equalized_op,equalized_adj1_ID,current_sample_length);


    SOVA_out=vector_calloc_int(current_sample_length);


    llr=vector_calloc(current_sample_length);

    /************Include for VA***********************/
     va_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,current_sample_length);
    /************Include for VA***********************/

    llr_length=(long int *)malloc( sizeof(long int));

    /************Include for SOVA***********************/
    //sova_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,llr,llr_length,current_sample_length);

    //current_sample_length=*llr_length;

    /************Include for SOVA***********************/


    /************Include for MAP***********************/
    //map_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,llr,llr_length,current_sample_length,SNR);

    //current_sample_length=*llr_length;

    /************Include for MAP***********************/


    /************OLD CODE***********************/

    //va_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,current_sample_length);

    /************OLD CODE***********************/

    ber_val=ber_compute(ber_bits,SOVA_out,current_sample_length);


   free (equalized_op);

   free (equalizer_main);
   free (equalizer_adj1_ID);
   free (equalizer_adj1_OD);


   free (equalized_main);
   free (equalized_adj1_ID);
   free (equalized_adj1_OD);


   free (target);

    free (SOVA_out);
    free (llr_length);
   free (llr);

    return ber_val;
}

long double onesided_equalize_oneD_SOVA(int **tx_merge,int current_sample_length, long double *rx_main,long double *rx_adj1_OD,int Nb_eq_main,int Nb_eq_adj1_OD,int Nb_target,int *ber_bits,float SNR)
{

  long double *equalizer_main,*equalizer_adj1_OD,*target;

  int track_count;

  long double *equalized_main,*equalized_adj1_OD,*equalized_op;

  long double ber_val;

    int MemoryDepth;
    int NbStates;
    int *SOVA_out;


    long double *llr;
    long int *llr_length;


    /**********Equalizer computation*********/
    track_count=1;          // Assuming 1 track on each side influence
    /*************************************/

    /**********SOVA Parameters*********/
    //Currently har coded for 4 states
    MemoryDepth=40;
    NbStates=4;
    /*************************************/


/********************EQUALIZER*************************/

    equalizer_main= vector_calloc(Nb_eq_main);
    equalizer_adj1_OD= vector_calloc(Nb_eq_adj1_OD);

    target= vector_calloc(Nb_target);

    //Compute Equalizer and Target Coeff
    lagrange_err_min_onesided(equalizer_main,equalizer_adj1_OD,Nb_eq_main,Nb_eq_adj1_OD,target,Nb_target,tx_merge,current_sample_length,rx_main,rx_adj1_OD);

    //Equalize Main
    equalized_main=vector_calloc(current_sample_length+Nb_eq_main);
    convolution(rx_main,equalizer_main,current_sample_length,Nb_eq_main,equalized_main);

    //Equalize Adj1 OD
    equalized_adj1_OD=vector_calloc(current_sample_length+Nb_eq_adj1_OD);
    convolution(rx_adj1_OD,equalizer_adj1_OD,current_sample_length,Nb_eq_adj1_OD,equalized_adj1_OD);

      //truncate initial samples and end values from convolution process
    array_resize(equalized_main,Nb_eq_main-1,current_sample_length);

    //truncate initial samples and end values from convolution process
    array_resize(equalized_adj1_OD,Nb_eq_adj1_OD-1,current_sample_length);


    //!!!!!!!!!!!!!!!!!Have to see how to handle variable equalizer length??!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Length of sample is sample_length-(Filter_coeff-1);
    current_sample_length=current_sample_length-(Nb_eq_main-1);



//Under investigation!!! how to generate states for variable target size
//int MemoryDepth=200, NbStates=4;
    /*int Polynom1=3;
    int Polynom2=1;
    int Nb_states;
    int **inputmat,**nextstate,**outmat;

    Nb_states=pow(2,no_of_target_vals-1);

    inputmat=matrix_calloc_int(Nb_states,Nb_states);
    nextstate=matrix_calloc_int(Nb_states,2);
    outmat=matrix_calloc_int(Nb_states,2);

    gene_viterbi_matrix (Polynom1,Polynom2,inputmat,nextstate,outmat);
*/

    //Add Equalized outputs
    equalized_op=vector_calloc(current_sample_length);

    array_addition (equalized_op,equalized_main,equalized_adj1_OD,current_sample_length);


    SOVA_out=vector_calloc_int(current_sample_length);

    llr=vector_calloc(current_sample_length);

    /************Include for VA***********************/
    va_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,current_sample_length);
    /************Include for VA***********************/

    llr_length=(long int *)malloc( sizeof(long int));

    /************Include for SOVA***********************/
    //sova_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,llr,llr_length,current_sample_length);

    //current_sample_length=*llr_length;

    /************Include for SOVA***********************/


    /************Include for MAP***********************/
    //map_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,llr,llr_length,current_sample_length,SNR);

    //current_sample_length=*llr_length;

    /************Include for MAP***********************/


    /************OLD CODE***********************/
    //va_detection(MemoryDepth,NbStates,target,Nb_target,equalized_op, SOVA_out,current_sample_length);
    /************OLD CODE***********************/

    ber_val=ber_compute(ber_bits,SOVA_out,current_sample_length);


   free (equalized_op);

   free (equalizer_main);
   free (equalizer_adj1_OD);


   free (equalized_main);
   free (equalized_adj1_OD);


   free (target);

    free (SOVA_out);

       free (llr_length);
   free (llr);

    return ber_val;
}











