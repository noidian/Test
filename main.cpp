#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <sstream>

#include "init.h"   /*initialization function deifnitions*/
#include "impulse.h"   /* TDMR Parameters*/
#include "generic_func.h" //Contains generic fucntions for files etc.
#include "math_func.h" //Contains generic math functions
#include "comm_func.h"  //Contains Communication System Fuunctions

//Receiver model alone...Transmitter to be checked later.

#include "LDPC_mod2.h"

//TMR Model
#define TMR_Main 0 //Percentage
#define TMR_OD 0 //Percentage
#define TMR_ID 0 //Percentage
#define TMR_Main_Enable 1 //0 is disabled , 1 is Enabled
#define TMR_OD_Enable 0 //0 is disabled , 1 is Enabled
#define TMR_ID_Enable 0 //0 is disabled , 1 is Enabled

//#define ITI 25
//check impulse.h for P_cross to modify
/****Specify equalizer*********************/
//#define oneD_EQ 0// 0 implies false 1 implies selected
//#define twoD_EQ 0// 0 implies false 1 implies selected
//#define onesided_EQ 1//1 implies it is selected only with full flow i.e. not with fixed input
//#define onesided_offset 10 // is (1sided_offset/100)*track_pitch ITI = 20 => 15% offset, 27 => 5% offset
/****Specify simulation type*********************/
#define fixed_input 0 //0 implies using full program

/****Path to save final BER result*********************/
//#define save_path "ber_SMR_1Sided_EQ_11tap_1DVA_ITI_25_full_offset_10_TMR_0_LDPC_9_10.txt"
/****Path to save final BER result*********************/

//MATLAB precision is a little better meaning the BER comes out better in comparison to the C program. but both are identical for high SNR's
//i.e. 1D Equalizer above 22 both are equal below there is small difference creeps in...(matrix computations needs to be better i.e. precision)

/*SOVA metric normalization is not implemented since for current test lengths there is no overflow for long double*/

using namespace std;


int ITI,onesided_offset;
char save_path[300];
bool oneD_EQ,twoD_EQ,onesided_EQ;

int fix_input();

int full_flow();

int smr_parameters();

template <class T>
bool from_string(T& t,const std::string& s,std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

int main()
{

    smr_parameters();
    if (fixed_input==1)
        {
            fix_input();
        }
    else
        {
            full_flow();
        }

}

int smr_parameters()
{
    ifstream smr_param;

    smr_param.open("smr_parameters.cfg");

    string temp;        //temporary string to parse sim_param
    string check ("="); //sub string to match with
    size_t found_equal;

    string str_ITI,str_onesided_offset,str_save_path,str_oneD_EQ,str_twoD_EQ,str_onesided_EQ;

    int line_count=0;
    int str_length;



    while (smr_param.good())
    {
        getline(smr_param,temp);
        found_equal=temp.find(check);
        switch(line_count)
        {
            case 0: str_ITI=temp.substr(found_equal+1,temp.length());
                    from_string<int>(ITI, std::string(str_ITI), std::dec);
                    break;
            case 1: str_oneD_EQ=temp.substr(found_equal+1,temp.length());
                    from_string<bool>(oneD_EQ, std::string(str_oneD_EQ), std::dec);
                    break;
            case 2: str_onesided_EQ=temp.substr(found_equal+1,temp.length());
                    from_string<bool>(onesided_EQ, std::string(str_onesided_EQ), std::dec);
                    break;
            case 3: str_twoD_EQ=temp.substr(found_equal+1,temp.length());
                    from_string<bool>(twoD_EQ, std::string(str_twoD_EQ), std::dec);
                    break;
            case 4: str_onesided_offset=temp.substr(found_equal+1,temp.length());
                    from_string<int>(onesided_offset, std::string(str_onesided_offset), std::dec);
                    break;
            case 5: str_save_path=temp.substr(found_equal+1,temp.length());
                    strcpy(save_path,str_save_path.c_str());
                    break;
        }
        line_count++;
    }

    smr_param.close();

    return OK;
}


int full_flow()
{
    int *a;         //input data pointer for main track
    int *b,*c;         //input data pointer for adj ID,OD tracks
    int *d,*e;         //input data pointer for adj2 ID,OD tracks
    int **a_merge;
    int current_sample_length; //sample length
    float SNR;
    char *fname_input_main,*fname_input_adj1_OD,*fname_input_adj1_ID,*fname_input_adj2_ID,*fname_input_adj2_OD, *current_wd,str_SNR[5];
    long double *h_main,*h_adj; //Main track response

    /****Input Filter*****/
    long double *data_main_net,*data_adj1_OD_net,*data_adj1_ID_net;

    int mid_filter_val;
    long double *noised_main,*matched_main,*noised_adj1_OD,*matched_adj1_OD,*noised_adj1_ID,*matched_adj1_ID;


    /****Trackwise dimension specifications*******/
    //Main Track
    long double *main_h_n1,*main_h_1,*main_h;
    long double *main_PW50_cross_n1,*main_Ty_n1;
    long double *main_PW50_cross_1,*main_Ty_1;
    long double *main_PW50_rpm_along;
    long double main_track_pitch;

    //TMR parameter for main track impulse response
    float TMR_offset;

    //OD Track
    long double *od_h_n1,*od_h_1,*od_h;
    long double *od_PW50_cross_n1,*od_Ty_n1;
    long double *od_PW50_cross_1,*od_Ty_1;
    long double *od_PW50_rpm_along;
    long double od_track_pitch;


    //ID Track
    long double *id_h_n1,*id_h_1,*id_h;
    long double *id_PW50_cross_n1,*id_Ty_n1;
    long double *id_PW50_cross_1,*id_Ty_1;
    long double *id_PW50_rpm_along;
    long double id_track_pitch;

    //LLR Allocation
    long double *llr;
    long int *llr_length;

    llr_length=(long int *)malloc( sizeof(long int));


    //LDPC class
    LDPC_mod2 ldpc;

    int *a_encoded;         //input data pointer for main track

    char *a_uncoded;
    int *b_uncoded,*c_uncoded;
    int *d_uncoded,*e_uncoded;

    /**************LDPC Variables *************/
    int i;
    char *a_op_ber;
    int res,sum_block_error=0;

    /**********BER Computation*********/
    int *a_det;
    long double ber=0.0;
    long double ldpc_ber=0.0;

    char *ber_filename;
    /*************************************/

    ber_filename= new char [300];



    /**********Equalizer computation*********/
    int no_equalizer_coeff;
    int no_of_target_vals;
    /*************************************/

    long double *PW50_rpm_along;
    PW50_rpm_along=(long double *)calloc(1, sizeof (long double ));

   /*****************Main Track paramters variables********************/
    main_PW50_cross_1=(long double *)calloc(1, sizeof (long double ));
    main_PW50_cross_n1=(long double *)calloc(1, sizeof (long double ));
    main_Ty_n1=(long double *)calloc(1, sizeof (long double ));
    main_Ty_1=(long double *)calloc(1, sizeof (long double ));
    main_PW50_rpm_along=(long double *)calloc(1, sizeof (long double ));


    /*****************Main Track paramters variables********************/
    od_PW50_cross_1=(long double *)calloc(1, sizeof (long double ));
    od_PW50_cross_n1=(long double *)calloc(1, sizeof (long double ));
    od_Ty_n1=(long double *)calloc(1, sizeof (long double ));
    od_Ty_1=(long double *)calloc(1, sizeof (long double ));
    od_PW50_rpm_along=(long double *)calloc(1, sizeof (long double ));


     /*****************Main Track paramters variables********************/
    id_PW50_cross_1=(long double *)calloc(1, sizeof (long double ));
    id_PW50_cross_n1=(long double *)calloc(1, sizeof (long double ));
    id_Ty_n1=(long double *)calloc(1, sizeof (long double ));
    id_Ty_1=(long double *)calloc(1, sizeof (long double ));
    id_PW50_rpm_along=(long double *)calloc(1, sizeof (long double ));


    //Get input file name
   fname_input_main = new char [200];           //Increase depending on working directory length
   fname_input_adj1_OD = new char [200];           //Increase depending on working directory length
   fname_input_adj1_ID = new char [200];           //Increase depending on working directory length
   fname_input_adj2_OD = new char [200];
   fname_input_adj2_ID = new char [200];

   current_wd = new char [200];

   /*************LDPC Initialization***************/
    ldpc.initial_matrix_file();
    //read in H and G file:
    ldpc.read_pchk(ldpc.pchk_file);
    ldpc.read_gen(ldpc.gen_file,0,0);

    //initial decoding space:
    ldpc.initial_decoding_space();
   /*************LDPC Initialization***************/

    ldpc.received_block=(double *)calloc(ldpc.N,sizeof *ldpc.received_block);
    ldpc.iter_distr=(int *)calloc(ldpc.max_iter,sizeof(int));
    //Given a sample length then for a code rate use the maximum possible input data bits..

    int Nb_code_blocks,current_enc_sample_length;

    //Reverse computation
    float net_sample_length;
    net_sample_length=(sample_length)*((float)ldpc.info_length/(float)ldpc.N);
    current_enc_sample_length=floor(net_sample_length);

    Nb_code_blocks=current_enc_sample_length/ldpc.info_length;

    current_sample_length=Nb_code_blocks*ldpc.N; //sample_length of output codeword
    current_enc_sample_length=Nb_code_blocks*ldpc.info_length; //sample_length of input codeword

    /***************Pad input codeblocks to ensure codeblockwise decoding***************/

    int *pad_values;
    int pad_length;

    no_equalizer_coeff=11;

    pad_length=(Filter_coeff-1)/2+(Filter_coeff-1)/2+(no_equalizer_coeff-1)/2+1; //input filter+output filter + equalizer  + 1 bit loss +1 ??

    pad_values   = (int *)calloc(pad_length,sizeof *pad_values);
    for (i=0;i<pad_length;i++)
          pad_values[i]=-1;

    /**********************************************************************************/


    getcwd(current_wd,200);

    //strcpy(fname_input_main,current_wd);
    shiftback_wd(current_wd,200,fname_input_main,1);          //Only for GCC compiler Linux PCs
    strcat(fname_input_main,"/C_Input_Files/input_main.txt");   //input main data file

//    current_sample_length=sample_length; //sample_length defined in init.h

    /******Main Track*************/
    a_uncoded = vector_calloc_char(current_enc_sample_length);
    //Copy input information bits to variable "a"
    fileopen(fname_input_main,a_uncoded,current_enc_sample_length);
    /******Main Track*************/

    shiftback_wd(current_wd,200,fname_input_adj1_OD,1);          //Only for GCC compiler Linux PCs
    strcat(fname_input_adj1_OD,"/C_Input_Files/input_adj1_OD.txt");   //input main data file
    /******Adjacent OD Track*************/
    c_uncoded = vector_calloc_int(current_sample_length+2*pad_length);
    //Copy input information bits to variable "c"
    fileopen(fname_input_adj1_OD,c_uncoded,current_sample_length+2*pad_length);
    /******Adjacent OD Track*************/

    shiftback_wd(current_wd,200,fname_input_adj1_ID,1);          //Only for GCC compiler Linux PCs
    strcat(fname_input_adj1_ID,"/C_Input_Files/input_adj1_ID.txt");   //input main data file
   /******Adjacent ID Track*************/
   b_uncoded = vector_calloc_int(current_sample_length+2*pad_length);
    //Copy input information bits to variable "b"
   fileopen(fname_input_adj1_ID,b_uncoded,current_sample_length+2*pad_length);
    /******Adjacent ID Track*************/




    if (twoD_EQ==1)
    {
        shiftback_wd(current_wd,200,fname_input_adj2_OD,1);          //Only for GCC compiler Linux PCs
        strcat(fname_input_adj2_OD,"/C_Input_Files/input_adj2_OD.txt");   //input adj2_OD data file
        /******Adjacent OD Track*************/
        e_uncoded = vector_calloc_int(current_sample_length+2*pad_length);
        //Copy input information bits to variable "b"
        fileopen(fname_input_adj2_OD,e_uncoded,current_sample_length+2*pad_length);
        /******Adjacent OD Track*************/

        shiftback_wd(current_wd,200,fname_input_adj2_OD,1);          //Only for GCC compiler Linux PCs
        strcat(fname_input_adj2_ID,"/C_Input_Files/input_adj2_ID.txt");   //input adj2_OD data file
        /******Adjacent OD Track*************/
        d_uncoded = vector_calloc_int(current_sample_length+2*pad_length);
        //Copy input information bits to variable "b"
        fileopen(fname_input_adj2_OD,d,current_sample_length+2*pad_length);
       /******Adjacent OD Track*************/

    }
    else if(onesided_EQ==1)
    {
           shiftback_wd(current_wd,200,fname_input_adj2_OD,1);          //Only for GCC compiler Linux PCs
        strcat(fname_input_adj2_OD,"/C_Input_Files/input_adj2_OD.txt");   //input adj2_OD data file
        /******Adjacent OD Track*************/
        e_uncoded = vector_calloc_int(current_sample_length+2*pad_length);
        //Copy input information bits to variable "b"
        fileopen(fname_input_adj2_OD,e_uncoded,current_sample_length+2*pad_length);

    }


        //free memory
        delete [] fname_input_adj2_OD;
        delete [] fname_input_adj2_ID;
        delete [] fname_input_main;
        delete [] fname_input_adj1_OD;
        delete [] fname_input_adj1_ID;


/************************LDPC ENCODE**********************/
    //If LDPC is used:
    //Change above variables a,b,c to a_uncoded, b_uncoded, c_uncoded

    //Encode only main track
    //use a,b,c as ldpc output variables
    a_encoded = vector_calloc_int(current_sample_length);
    a_op_ber = vector_calloc_char(current_sample_length);
    //b = vector_calloc_int(current_sample_length);
    //c = vector_calloc_int(current_sample_length);

    b=b_uncoded;
    c=c_uncoded;
    d=d_uncoded;
    e=e_uncoded;

    char *encoded_data,*source_data;

	encoded_data = (char *)calloc(ldpc.N,sizeof *encoded_data);
	source_data   = (char *)calloc(ldpc.info_length,sizeof *source_data);

//for(i=0;i<ldpc.info_length;i++)
//            source_data[i]=0;//rand()%2;

    ldpc.encoded_data=(char *)calloc(ldpc.N,sizeof *ldpc.encoded_data);
    ldpc.source_data=(char *)calloc(ldpc.N,sizeof *ldpc.source_data);

    unsigned long long int sample_start_pt,coded_start_pt;
    unsigned long long int k,j;
    //encoding
    int chk;

    for (i=0;i<Nb_code_blocks;i++)
	{

        sample_start_pt=i*ldpc.info_length;
        //Send for encoding block by block


        ldpc.sparse_encode(a_uncoded,encoded_data,sample_start_pt);  //

        coded_start_pt=i*ldpc.N;
        k=0;

        chk = ldpc.check(ldpc.H, encoded_data, ldpc.parity_checks); //1 for correct; 0 for fails.
        if(chk==0)
			{
				printf("Encoding Wrong!\n");
			}
//			else
//				printf("Encoding check is correct.\n");

        for (j=coded_start_pt;j<((i+1)*ldpc.N);j++)
            {
                a_encoded[j]=(int)encoded_data[k];
                a_op_ber[j]=encoded_data[k];
                if (a_encoded[j]==0)
                 {
                    a_encoded[j]=-1;
                 }
                k++;
            }
	}

    //free corresponding memory
    free (encoded_data);
    free(source_data);



/*****************LDPC Encoded Data***********************/

/************************PAD encoded codeblock at start and end***********************/

    a= vector_calloc_int(current_sample_length+2*pad_length);

    array_append(a,pad_values,a_encoded,pad_length,current_sample_length+2*pad_length);

    //current_sample_length=current_sample_length+2*pad_length; //check if this is ok...??~~!!!!


//    int temp1=0;
//    int chk2[40];
//    long int loop1=0;
//    for (loop1=current_sample_length;loop1>current_sample_length-40;loop1--)
//        {
//            chk2[temp1]=a_encoded[loop1];
//            temp1++;
//
//        }
//

    free(a_encoded);




//
//    /***********Cross track PW 50 and bit width deifnition by R.Wood**********/
              //track_pitch=6.4;
//            //Define PW50_cross and Ty according to required ITI
//            //As=0.2774
//            //Max h[0] amp =0.8686 and h[1]=0.34
//            *PW50_cross_n1=4.8; //OD track
//            *Ty_n1=6.4;
//            *PW50_cross_1=4.8; //ID track
//            *Ty_1=6.4;           //was 5.7nm track pitch for 0.27 amplitude
//            //0.21
//    /********************RWood*******************************************/
    /********************Mine********************************************/
    /***********Cross track PW 50 and bit width deifnition**********/
            //track_pitch=5.6;
            //Define PW50_cross and Ty according to required ITI
            //0.1984
           // *PW50_cross_n1=4.2; //OD track
           // *Ty_n1=5.6;
          //  *PW50_cross_1=4.2; //ID track
           // *Ty_1=5.6;
            //0.21
    /********************Mine********************************************/
    //Track pitch parameters (typically should be equal)
    main_track_pitch=5.6;
    od_track_pitch=5.6;
    id_track_pitch=5.6;


    /***********Main Track Filter Initialization**********/
    main_h=vector_calloc(Filter_coeff);    //_c denotes cross track _a denotes along track
    main_h_n1=vector_calloc(Filter_coeff);
    main_h_1=vector_calloc(Filter_coeff);

    main_impulse_response_parameter_init(main_PW50_cross_n1,main_Ty_n1,main_PW50_cross_1,main_Ty_1,TMR_Main,main_track_pitch,main_h_n1,main_h,main_h_1);

    /***********OD Track********************************************/
    od_h=vector_calloc(Filter_coeff);  //OD filter init
    od_h_n1=vector_calloc(Filter_coeff);
    od_h_1=vector_calloc(Filter_coeff);

    od_impulse_response_parameter_init(od_PW50_cross_n1,od_Ty_n1,od_PW50_cross_1,od_Ty_1,TMR_OD,od_track_pitch,od_h_n1,od_h,od_h_1);

    /***********ID Track********************************************/
    id_h=vector_calloc(Filter_coeff);    //ID filter init
    id_h_n1=vector_calloc(Filter_coeff);
    id_h_1=vector_calloc(Filter_coeff);

    id_impulse_response_parameter_init(id_PW50_cross_n1,id_Ty_n1,id_PW50_cross_1,id_Ty_1,TMR_ID,id_track_pitch,id_h_n1,id_h,id_h_1);

    if (twoD_EQ==1)
    {

    data_adj1_OD_net = vector_calloc(current_sample_length+2*pad_length-(Filter_coeff-1));
    input_ITI_adder(data_adj1_OD_net,e,od_h_n1,c,od_h,a,od_h_1,current_sample_length+2*pad_length,Filter_coeff);

    data_adj1_ID_net = vector_calloc(current_sample_length+2*pad_length-(Filter_coeff-1));
    input_ITI_adder(data_adj1_ID_net,d,id_h_n1,b,id_h,a,id_h_1,current_sample_length+2*pad_length,Filter_coeff);

    data_main_net = vector_calloc(current_sample_length+2*pad_length-(Filter_coeff-1));
    input_ITI_adder(data_main_net,c,main_h_n1,a,main_h,b,main_h_1,current_sample_length+2*pad_length,Filter_coeff);

    //free memory
    free(d);
    free(e);
    }
    else if (oneD_EQ==1)
    {
        data_main_net = vector_calloc(current_sample_length+2*pad_length-(Filter_coeff-1));
        input_ITI_adder(data_main_net,c,main_h_n1,a,main_h,b,main_h_1,current_sample_length+2*pad_length,Filter_coeff);
    }

    else if (onesided_EQ==1)
    {

       data_adj1_OD_net = vector_calloc(current_sample_length+2*pad_length-(Filter_coeff-1));
       input_ITI_adder(data_adj1_OD_net,e,od_h_n1,c,od_h,a,od_h_1,current_sample_length+2*pad_length,Filter_coeff);

//     *******************Rwood with 20%*******************************************
//        *********Cross track PW 50 and bit width deifnition********
//            Define PW50_cross and Ty according to required ITI
//            0.1984
//            *PW50_cross_n1=4.8; //OD track       //0.6nm offset  for best performance
//            *Ty_n1=5.8;
//            *PW50_cross_1=4.8; //ID track
//            *Ty_1=7.0;
//            0.21
//        *********Cross track PW 50 and bit width deifnition********
//   ******************Rwood******************************************

     /********************Mine********************************************/
        /**********Cross track PW 50 and bit width deifnition*********/
            //Define PW50_cross and Ty according to required ITI
            //0.1984
           // track_pitch=5.0;        //offset track pitch
           // *PW50_cross_n1=4.2; //OD track       //0.6nm offset  for best performance
           // *Ty_n1=5.0;
           // *PW50_cross_1=4.2; //ID track
           // *Ty_1=6.2;
            //0.21
        /**********Cross track PW 50 and bit width deifnition*********/
   /*******************Mine*******************************************/

    if (ITI==20)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        //0.1984main_track_pitch
        *main_PW50_cross_n1=4.2; //OD track         //0.6nm offset  for best performance
        *main_Ty_n1=(1-(float)(onesided_offset/100))*main_track_pitch;
        *main_PW50_cross_1=4.2; //ID track
        *main_Ty_1=(1+(float)(onesided_offset/100))*main_track_pitch;
        //0.21
        /********************Mine********************************************/
    }
    else if (ITI==25)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=4.5; //OD track
        *main_Ty_n1=(1-(float)(onesided_offset/100))*main_track_pitch;
        *main_PW50_cross_1=4.5; //ID track
        *main_Ty_1=(1+(float)(onesided_offset/100))*main_track_pitch;
        //0.2444
        /********************Mine********************************************/
    }
    else if (ITI==27)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=4.8; //OD track
        *main_Ty_n1=(1-(float)(onesided_offset/100))*main_track_pitch;
        *main_PW50_cross_1=4.8; //ID track
        *main_Ty_1=(1+(float)(onesided_offset/100))*main_track_pitch;
        //0.2898
        /********************Mine********************************************/
    }
    else if (ITI==15)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=3.9; //OD track
        *main_Ty_n1=(1-(float)(onesided_offset/100))*main_track_pitch;
        *main_PW50_cross_1=3.9; //ID track
        *main_Ty_1=(1+(float)(onesided_offset/100))*main_track_pitch;
        //0.1532
        /********************Mine********************************************/
    }
    else if (ITI==10)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=3.5; //OD track
        *main_Ty_n1=(1-(float)(onesided_offset/100))*main_track_pitch;
        *main_PW50_cross_1=3.5; //ID track
        *main_Ty_1=(1+(float)(onesided_offset/100))*main_track_pitch;
        //0.0974
        /********************Mine********************************************/
    }
    else if (ITI==5)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=3.1; //OD track
        *main_Ty_n1=(1-(float)(onesided_offset/100))*main_track_pitch;
        *main_PW50_cross_1=3.1; //ID track
        *main_Ty_1=(1+(float)(onesided_offset/100))*main_track_pitch;
        //0.0513
        /********************Mine********************************************/
    }

    tmr_addition(TMR_Main,main_Ty_n1,main_Ty_1,main_track_pitch,*main_Ty_n1,*main_Ty_1);    //1.2 indicates the 0.6nm offset from center

    smr_impulse_response_2D(main_h_n1,main_h,main_h_1,PW50_along,Tx,*main_PW50_cross_n1,*main_Ty_n1,*main_PW50_cross_1,*main_Ty_1,Amplitude,Filter_coeff);
    data_main_net = vector_calloc(current_sample_length+2*pad_length-(Filter_coeff-1));
    input_ITI_adder(data_main_net,c,main_h_n1,a,main_h,b,main_h_1,current_sample_length+2*pad_length,Filter_coeff);

    free(e);

    }

    //Need to resize input matrix since the filtering at the input has eliminated Filter_Coeff-1 values
    mid_filter_val=(Filter_coeff-1)/2;  //works ok for Filter Coeff 7 and odd value

    array_resize(a,mid_filter_val,current_sample_length+2*pad_length-mid_filter_val);
    array_resize(b,mid_filter_val,current_sample_length+2*pad_length-mid_filter_val);
    array_resize(c,mid_filter_val,current_sample_length+2*pad_length-mid_filter_val);

   //Length of sample is sample_length-(Filter_coeff-1);
    current_sample_length=current_sample_length-(Filter_coeff-1);


    //Free memory for input filter
//    free (h);
//    free (h_n1);
//    free (h_1);


    //Free memory for Main input filter
    free (main_h);
    free (main_h_n1);
    free (main_h_1);

    free (od_h);
    free (od_h_n1);
    free (od_h_1);


    free (id_h);
    free (id_h_n1);
    free (id_h_1);




    //OD
    free(od_PW50_cross_1);
    free(od_PW50_cross_n1);
    free(od_Ty_n1);
    free(od_Ty_1);
    free(od_PW50_rpm_along);


    //ID
    free(id_PW50_cross_1);
    free(id_PW50_cross_n1);
    free(id_Ty_n1);
    free(id_Ty_1);
    free(id_PW50_rpm_along);





    *PW50_rpm_along=5.2; //NORMAL 5.2

  /********************Works for Mine + R.Wood********************************************/
    /***************************ISI**************************/
    /* ISI    PW50    RPM                                         */
    /* 1.1    5.72    5940                                      */
    /* 1.3    6.76    7020                                      */
    /* 1.4    7.28    7560                                      */
    /* 1.5    7.8     8100                                      */
    /* 1.6    8.32    8640                                      */
    /* 1.7    8.84    9180                                      */
    /* 1.9    9.88    10260                                      */
    /***************************ISI**************************/


    //TMR One sided
    if (onesided_EQ==1)
    {
        TMR_offset=(((float)onesided_offset/100)-((float)TMR_Main/100))*main_track_pitch;
    }
    else
    {
     //TMR normal
        TMR_offset=((float)TMR_Main/100)*main_track_pitch;
    }



      h_main=vector_calloc(Filter_coeff);    //_c denotes cross track _a denotes along track
      h_adj=vector_calloc(Filter_coeff);    //_c denotes cross track _a denotes along track
     smr_impulse_response_1D(h_main,*PW50_rpm_along,Tx,Amplitude,TMR_offset,Filter_coeff);
     smr_impulse_response_1D(h_adj,*PW50_rpm_along,Tx,Amplitude,0,Filter_coeff);


    /**********1D Equalization******************************/
    if (oneD_EQ==1)
    {       no_equalizer_coeff=11;
             no_of_target_vals=3;
    }
    /**********1D Equalization******************************/

    /********************2D Equalization********************/
    else if (twoD_EQ==1)
        {
           no_equalizer_coeff=11;
           no_of_target_vals=3;
        }
    else if (onesided_EQ==1)
    {
            no_equalizer_coeff=11;
           no_of_target_vals=3;

    }


    /********************1D EQ Matched filter correction********************/
    //The matched filter output is less so for doing the computation the input matrix is resized
    array_resize(a,mid_filter_val,current_sample_length+2*pad_length-mid_filter_val);


    /********************2D Equalization********************/

    if (twoD_EQ==1)
    {
     array_resize(b,mid_filter_val,current_sample_length+2*pad_length-mid_filter_val);
    array_resize(c,mid_filter_val,current_sample_length+2*pad_length-mid_filter_val);

    //Tabulate input a into a common matrix for post matched computation
    a_merge=matrix_calloc_int(3,current_sample_length+2*pad_length-(Filter_coeff-1));
    twoD_array_copy_int(a_merge,c,0,current_sample_length+2*pad_length-(Filter_coeff-1));
    twoD_array_copy_int(a_merge,a,1,current_sample_length+2*pad_length-(Filter_coeff-1));
    twoD_array_copy_int(a_merge,b,2,current_sample_length+2*pad_length-(Filter_coeff-1));
    }
    else if (onesided_EQ==1)
    {
        array_resize(c,mid_filter_val,current_sample_length+2*pad_length-mid_filter_val);
        //Tabulate input a into a common matrix for post matched computation
        a_merge=matrix_calloc_int(2,current_sample_length+2*pad_length-(Filter_coeff-1));
        twoD_array_copy_int(a_merge,c,0,current_sample_length+2*pad_length-(Filter_coeff-1));
        twoD_array_copy_int(a_merge,a,1,current_sample_length+2*pad_length-(Filter_coeff-1));


    }



    //free memory no need of b and c input data for 1D or 2D
    free (b);
    free (c);

//Note here the array of "input" is less than 6 in comparison to output to compensate for matching operation.

/*******************************Post detection input bits for ber *******************/


  /***************Separate comparison data for BER***************/
   //Need to resize input matrix since the equalization has eliminated no_equalizer_coeff-1 values and same at the equalized_op as well

   a_det=vector_calloc_int(current_sample_length+2*pad_length);
   vector_copy (a_det,a,current_sample_length+2*pad_length);


   mid_filter_val=ceil((double)(no_equalizer_coeff-1)/2);      //works for odd and even equalizer coefficients

   if((no_equalizer_coeff%2)!=2) //only for odd number of eq coefficients
   {
       mid_filter_val=mid_filter_val+1; //there is an anamoly ..loss of 1 bit in the viterbi..to be debugged...sometime..!!
   }

   array_resize(a_det,mid_filter_val,current_sample_length+2*pad_length-mid_filter_val);
    //Note here the array of "input" is less than 6 in comparison to output to compensate for equalization operation.
 /***************************************************************/


    if ((twoD_EQ==1) || (onesided_EQ==1))
    {//memory free only a_det adn a_merge are required. for 1D a_merge=a;
        free (a);
    }

//ber filename





for (SNR=25.0;SNR>=6.0;SNR--)
//for (SNR=13;SNR>8;) //RS 1D Model
{
    /* Noisy Input Filename Formation*/
    //snprintf(str_SNR,4,"%f",SNR);
    //dot2underscore(str_SNR,5);
    //dot2null(str_SNR,5);
     if (oneD_EQ==1)
    {

        //Noise addition
        noised_main = vector_calloc(current_sample_length+2*pad_length);
        awgn(noised_main,data_main_net,current_sample_length+2*pad_length,SNR);

        //Matched Filter output = sample length + filter length
        matched_main = vector_calloc(current_sample_length+2*pad_length+Filter_coeff);

        convolution(noised_main,h_main,current_sample_length+2*pad_length,Filter_coeff,matched_main);

        //truncate initial samples and end values from convolution process
        array_resize(matched_main,Filter_coeff-1,current_sample_length+2*pad_length);

       //Length of sample is sample_length-(Filter_coeff-1);
        current_sample_length=current_sample_length-(Filter_coeff-1);

        llr=vector_calloc(current_sample_length+2*pad_length);

        ber=oneD_equalize_oneD_SOVA(a,current_sample_length+2*pad_length,matched_main,no_equalizer_coeff,no_of_target_vals,a_det,llr_length,llr,SNR);
    }
    else if (twoD_EQ==1)
    {
        //Read Noisy File
        noised_main = vector_calloc(current_sample_length+2*pad_length);
        awgn(noised_main,data_main_net,current_sample_length+2*pad_length,SNR);

        //Matched Filter output = sample length + filter length
        matched_main = vector_calloc(current_sample_length+2*pad_length+Filter_coeff);
        convolution(noised_main,h_main,current_sample_length+2*pad_length,Filter_coeff,matched_main);

        //truncate initial samples and end values from convolution process
        array_resize(matched_main,Filter_coeff-1,current_sample_length+2*pad_length);


        /******Adj OD Track*************/
        //Read Noisy File

        noised_adj1_OD = vector_calloc(current_sample_length+2*pad_length);
        awgn(noised_adj1_OD,data_adj1_OD_net,current_sample_length+2*pad_length,SNR);

        //Matched Filter output = sample length + filter length
        matched_adj1_OD = vector_calloc(current_sample_length+2*pad_length+Filter_coeff);

        convolution(noised_adj1_OD,h_adj,current_sample_length+2*pad_length,Filter_coeff,matched_adj1_OD);

        //truncate initial samples and end values from convolution process
        array_resize(matched_adj1_OD,Filter_coeff-1,current_sample_length+2*pad_length);


        /******Adj ID Track*************/
        //Read Noisy File

        noised_adj1_ID = vector_calloc(current_sample_length+2*pad_length);
        awgn(noised_adj1_ID,data_adj1_ID_net,current_sample_length+2*pad_length,SNR);

        //Matched Filter output = sample length + filter length
        matched_adj1_ID = vector_calloc(current_sample_length+2*pad_length+Filter_coeff);

        convolution(noised_adj1_ID,h_adj,current_sample_length+2*pad_length,Filter_coeff,matched_adj1_ID);

        //truncate initial samples and end values from convolution process
        array_resize(matched_adj1_ID,Filter_coeff-1,current_sample_length+2*pad_length);


       //Length of sample is sample_length-(Filter_coeff-1);
        current_sample_length=current_sample_length-(Filter_coeff-1);

         llr=vector_calloc(current_sample_length+2*pad_length);

        ber=twoD_equalize_oneD_SOVA(a_merge,current_sample_length+2*pad_length,matched_main,matched_adj1_OD,matched_adj1_ID,no_equalizer_coeff,no_equalizer_coeff,no_equalizer_coeff,no_of_target_vals,a_det,llr_length,llr,SNR);


        free (matched_adj1_ID);
        free(noised_adj1_ID);
        free (matched_adj1_OD);
       free(noised_adj1_OD);
    }
       else if (onesided_EQ==1)
    {
        //Read Noisy File
        noised_main = vector_calloc(current_sample_length+2*pad_length);
        awgn(noised_main,data_main_net,current_sample_length+2*pad_length,SNR);

        //Matched Filter output = sample length + filter length
        matched_main = vector_calloc(current_sample_length+2*pad_length+Filter_coeff);
        convolution(noised_main,h_main,current_sample_length+2*pad_length,Filter_coeff,matched_main);

        //truncate initial samples and end values from convolution process
        array_resize(matched_main,Filter_coeff-1,current_sample_length+2*pad_length);


        /******Adj OD Track*************/
        //Read Noisy File

        noised_adj1_OD = vector_calloc(current_sample_length+2*pad_length);
        awgn(noised_adj1_OD,data_adj1_OD_net,current_sample_length+2*pad_length,SNR);

        //Matched Filter output = sample length + filter length
        matched_adj1_OD = vector_calloc(current_sample_length+2*pad_length+Filter_coeff);

        convolution(noised_adj1_OD,h_adj,current_sample_length+2*pad_length,Filter_coeff,matched_adj1_OD);

        //truncate initial samples and end values from convolution process
        array_resize(matched_adj1_OD,Filter_coeff-1,current_sample_length+2*pad_length);

        //Length of sample is sample_length-(Filter_coeff-1);
        current_sample_length=current_sample_length-(Filter_coeff-1);

        llr=vector_calloc(current_sample_length+2*pad_length);

        ber=onesided_equalize_oneD_SOVA(a_merge,current_sample_length+2*pad_length,matched_main,matched_adj1_OD,no_equalizer_coeff,no_equalizer_coeff,no_of_target_vals,a_det,llr_length,llr,SNR);

        free (matched_adj1_OD);
       free(noised_adj1_OD);
    }

    //LDPC Decoding-
    for (i=0;i<Nb_code_blocks;i++) //ignores last codeblock
	{


        //coded_start_pt=(i*ldpc.N);
        coded_start_pt=(i*ldpc.N);
        k=0;

        for (j=coded_start_pt;j<((i+1)*ldpc.N);j++)
            {
                ldpc.received_block[k]=llr[j];
                k++;
            }

         k=0;
         sample_start_pt=i*ldpc.info_length;
         for (j=sample_start_pt;j<((i+1)*ldpc.info_length);j++)
         {
             ldpc.source_data[k]=a_uncoded[j];
             k++;
         }


        res=ldpc.ldpc_decode_F();
        sum_block_error += 1 - res;
        ldpc_ber+=(long double)ldpc.info_bit_error_in_one_block;
	}
        ldpc_ber=ldpc_ber/(ldpc.info_length*Nb_code_blocks);


    free (llr);


    //Create LDPC decoded variable

//    res=ldpc.ldpc_decode_F();

//	sum_block_error += 1 - res;

    //dec_out=vector_calloc_int(llr_length));

    //free(dec_out);


    //restore the correct length back
    current_sample_length=current_sample_length+(Filter_coeff-1);

    free (matched_main);
    free(noised_main);

    strcpy(ber_filename,save_path);
    //filestore_ber(SNR,ber,ber_filename);
    //filestore_ber(ber,ber_filename);
    filestore_ber(ber,ldpc_ber,ber_filename);

    //SNR=SNR-0.1; //RS Model only

}
    /************Clear LLR length****************/
        free (llr_length);

    //Filter the noisy signal using the equalizer /coefficients


    //free (ber_filename); //Doesn't work??!!!
    if (oneD_EQ==1)
      {
          free (a); //a=a_merge
      }
    else if (twoD_EQ==1)
            {
                matrix_free (a_merge,3);

                free(data_adj1_OD_net);
                free(data_adj1_ID_net);
            }
    else if (onesided_EQ==1)
            {
                matrix_free (a_merge,2);

                free(data_adj1_OD_net);
            }

    free (h_main);
    free (h_adj);
    free(data_main_net);
    free(a_det);


    free(a_uncoded);
    free(a_op_ber);


    free(PW50_rpm_along);


    delete [] ber_filename;

    delete [] current_wd;

    return OK;

}

//Commented out for LDPC??!!Might require for corss validation i.e. debugging purposes with fixed data!!

//int fix_input()
//{
//        int *a;         //input data pointer for main track
//    int *b,*c;         //input data pointer for adj ID,OD tracks
//    int **a_merge;
//    int current_sample_length; //sample length
//    float SNR;
//    char *fname_input_main,*fname_input_adj1_OD,*fname_input_adj1_ID, *current_wd,str_SNR[5];
//    long double *h; //Main track response
//
//    int mid_filter_val;
//    long double *noised_main,*matched_main,*noised_adj1_OD,*matched_adj1_OD,*noised_adj1_ID,*matched_adj1_ID;
//
//    /****oneD_EQ==1 ******/
//
//
//    /****twoD_EQ==1*******/
//     long double *h_n1,*h_1;
//     long double PW50_cross_n1,Ty_n1;
//     long double PW50_cross_1,Ty_1;
//
//
//    /**********BER Computation*********/
//    int *a_det;
//    long double ber=0.0;
//
//    char *ber_filename;
//    /*************************************/
//
//    ber_filename= new char [200];
//
//
//
//    /**********Equalizer computation*********/
//    int no_equalizer_coeff;
//    int no_of_target_vals;
//    /*************************************/
//
//
//
//    //Get input file name
//   fname_input_main = new char [200];           //Increase depending on working directory length
//
//
//   current_wd = new char [200];
//
//    getcwd(current_wd,200);
//
//    //strcpy(fname_input_main,current_wd);
//    shiftback_wd(current_wd,200,fname_input_main,1);          //Only for GCC compiler Linux PCs
//    strcat(fname_input_main,"/C_Input_Files/input_main.txt");   //input main data file
//
//
//    current_sample_length=sample_length; //sample_length defined in init.h
//
//
//    /******Main Track*************/
//    a = vector_calloc_int(current_sample_length);
//    //Copy input information bits to variable "a"
//    fileopen(fname_input_main,a,current_sample_length);
//    /******Main Track*************/
//
//
//
//    /**********1D Equalization******************************/
//    if (oneD_EQ==1)
//    {
//         h=vector_calloc(Filter_coeff);    //_c denotes cross track _a denotes along track
//         smr_impulse_response_1D(h,PW50_along,Tx,Amplitude,0,Filter_coeff);
//
//             no_equalizer_coeff=11;
//             no_of_target_vals=3;
//    }
//    /**********1D Equalization******************************/
//
//    /********************2D Equalization********************/
//    else if (twoD_EQ==1)
//        {
//                 no_equalizer_coeff=7;
//                 no_of_target_vals=3;
//
//                fname_input_adj1_OD = new char [200];           //Increase depending on working directory length
//                fname_input_adj1_ID = new char [200];           //Increase depending on working directory length
//
//                shiftback_wd(current_wd,200,fname_input_adj1_OD,1);          //Only for GCC compiler Linux PCs
//                strcat(fname_input_adj1_OD,"/C_Input_Files/input_adj1_OD.txt");   //input main data file
//
//                /******Adjacent OD Track*************/
//                c = vector_calloc_int(current_sample_length);
//                //Copy input information bits to variable "a"
//                fileopen(fname_input_adj1_OD,c,current_sample_length);
//                /******Adjacent OD Track*************/
//
//                shiftback_wd(current_wd,200,fname_input_adj1_ID,1);          //Only for GCC compiler Linux PCs
//                strcat(fname_input_adj1_ID,"/C_Input_Files/input_adj1_ID.txt");   //input main data file
//
//                /******Adjacent ID Track*************/
//               b = vector_calloc_int(current_sample_length);
//               //Copy input information bits to variable "a"
//               fileopen(fname_input_adj1_ID,b,current_sample_length);
//               /******Adjacent ID Track*************/
//
//            /***********Cross track PW 50 and bit width deifnition**********/
//            //Define PW50_cross and Ty according to required ITI
//            //0.1984
//            PW50_cross_n1=4.2; //OD track
//            Ty_n1=5.6;
//            PW50_cross_1=4.2; //ID track
//            Ty_1=5.6;
//            //0.21
//
//            /***********Cross track PW 50 and bit width deifnition**********/
//
//            h=vector_calloc(Filter_coeff);    //_c denotes cross track _a denotes along track
//            h_n1=vector_calloc(Filter_coeff);    //_c denotes cross track _a denotes along track Adjacent ID track
//            h_1=vector_calloc(Filter_coeff);    //Adjacent OD track
//            smr_impulse_response_2D(h_n1,h,h_1,PW50_along,Tx,PW50_cross_n1,Ty_n1,PW50_cross_1,Ty_1,Amplitude,Filter_coeff);
//        }
//    /********************2D Equalization********************/
//
//    //Need to resize input matrix since the filtering at the input has eliminated Filter_Coeff-1 values and sam at the matched filter as well
//    mid_filter_val=(Filter_coeff-1)/2;  //works ok for Filter Coeff 7 and odd value
//
//    array_resize(a,mid_filter_val,current_sample_length-mid_filter_val);
//
//    if (twoD_EQ==1)
//    {
//            array_resize(b,mid_filter_val,current_sample_length-mid_filter_val);
//            array_resize(c,mid_filter_val,current_sample_length-mid_filter_val);
//    }
//
//    current_sample_length=sample_length-(Filter_coeff-1);   //To account for Filter_Coeff-1 loss at input (safe side-not required)
//
//
//    array_resize(a,mid_filter_val,current_sample_length-mid_filter_val);
//
//    if (twoD_EQ==1)
//    {
//     array_resize(b,mid_filter_val,current_sample_length-mid_filter_val);
//    array_resize(c,mid_filter_val,current_sample_length-mid_filter_val);
//
//    //Tabulate input a into a common matrix for post matched computation
//    a_merge=matrix_calloc_int(3,current_sample_length-(Filter_coeff-1));
//    twoD_array_copy_int(a_merge,c,0,current_sample_length-(Filter_coeff-1));
//    twoD_array_copy_int(a_merge,a,1,current_sample_length-(Filter_coeff-1));
//    twoD_array_copy_int(a_merge,b,2,current_sample_length-(Filter_coeff-1));
//
//        //free memory
//        free (b);
//        free (c);
//
//
//
//    }
////Note here the array of "input" is less than 6 in comparison to output to compensate for matching operation.
//
///*******************************Post detection input bits for ber *******************/
//
//
//  /***************Separate comparison data for BER***************/
//   //Need to resize input matrix since the equalization has eliminated no_equalizer_coeff-1 values and same at the equalized_op as well
//
//   a_det=vector_calloc_int(current_sample_length);
//   vector_copy (a_det,a,current_sample_length);
//
//    mid_filter_val=ceil((double)(no_equalizer_coeff-1)/2);      //works for odd and even equalizer coefficients
//
//   mid_filter_val=mid_filter_val+1; //there is an anamoly ..loss of 1 bit in the viterbi..to be debugged...sometime..!!
//   array_resize(a_det,mid_filter_val,current_sample_length-mid_filter_val);
//    //Note here the array of "input" is less than 6 in comparison to output to compensate for equalization operation.
// /***************************************************************/
//
//
//    if (twoD_EQ==1)
//    {//memory free only a_det adn a_merge are required.
//        free (a);
//    }
//
////ber filename
//
//
//
//for (SNR=22.0;SNR>=6.0;SNR--)
//{
//    /* Noisy Input Filename Formation*/
//    snprintf(str_SNR,4,"%f",SNR);
//    //dot2underscore(str_SNR,5);
//    dot2null(str_SNR,5);
//     if (oneD_EQ==1)
//    {
//        shiftback_wd(current_wd,200,fname_input_main,1);          //Only for GCC compiler Linux PCs
//        strcat(fname_input_main,"/C_Input_Files/main_");   //input main data file
//        strcat(fname_input_main,str_SNR);
//        strcat(fname_input_main,".txt");
//
//        //Read Noisy File
//
//        noised_main = vector_calloc(current_sample_length);
//        fileopen(fname_input_main,noised_main,current_sample_length);
//
//        //Matched Filter output = sample length + filter length
//        matched_main = vector_calloc(current_sample_length+Filter_coeff);
//
//        convolution(noised_main,h,current_sample_length,Filter_coeff,matched_main);
//
//        //truncate initial samples and end values from convolution process
//        array_resize(matched_main,Filter_coeff-1,current_sample_length);
//
//       //Length of sample is sample_length-(Filter_coeff-1);
//        current_sample_length=current_sample_length-(Filter_coeff-1);
//
//        ber=oneD_equalize_oneD_SOVA(a,current_sample_length,matched_main,no_equalizer_coeff,no_of_target_vals,a_det,SNR);
//    }
//    else if (twoD_EQ==1)
//    {
//        /******Main Track*************/
//        shiftback_wd(current_wd,200,fname_input_main,1);          //Only for GCC compiler Linux PCs
//        strcat(fname_input_main,"/C_Input_Files/main_");   //input main data file
//        strcat(fname_input_main,str_SNR);
//        strcat(fname_input_main,".txt");
//
//        /******Adjacent ID Track*************/
//        shiftback_wd(current_wd,200,fname_input_adj1_ID,1);          //Only for GCC compiler Linux PCs
//        strcat(fname_input_adj1_ID,"/C_Input_Files/adj1_ID_");   //input main data file
//        strcat(fname_input_adj1_ID,str_SNR);
//        strcat(fname_input_adj1_ID,".txt");
//
//        /******Adjacent OD Track*************/
//        shiftback_wd(current_wd,200,fname_input_adj1_OD,1);          //Only for GCC compiler Linux PCs
//        strcat(fname_input_adj1_OD,"/C_Input_Files/adj1_OD_");   //input main data file
//        strcat(fname_input_adj1_OD,str_SNR);
//        strcat(fname_input_adj1_OD,".txt");
//
//        /******Main Track*************/
//        //Read Noisy File
//
//        noised_main = vector_calloc(current_sample_length);
//        fileopen(fname_input_main,noised_main,current_sample_length);
//
//        //Matched Filter output = sample length + filter length
//        matched_main = vector_calloc(current_sample_length+Filter_coeff);
//
//        convolution(noised_main,h,current_sample_length,Filter_coeff,matched_main);
//
//        //truncate initial samples and end values from convolution process
//        array_resize(matched_main,Filter_coeff-1,current_sample_length);
//
//
//        /******Adj OD Track*************/
//        //Read Noisy File
//
//        noised_adj1_OD = vector_calloc(current_sample_length);
//        fileopen(fname_input_adj1_OD,noised_adj1_OD,current_sample_length);
//
//        //Matched Filter output = sample length + filter length
//        matched_adj1_OD = vector_calloc(current_sample_length+Filter_coeff);
//
//        convolution(noised_adj1_OD,h,current_sample_length,Filter_coeff,matched_adj1_OD);
//
//        //truncate initial samples and end values from convolution process
//        array_resize(matched_adj1_OD,Filter_coeff-1,current_sample_length);
//
//
//        /******Adj ID Track*************/
//        //Read Noisy File
//
//        noised_adj1_ID = vector_calloc(current_sample_length);
//        fileopen(fname_input_adj1_ID,noised_adj1_ID,current_sample_length);
//
//        //Matched Filter output = sample length + filter length
//        matched_adj1_ID = vector_calloc(current_sample_length+Filter_coeff);
//
//        convolution(noised_adj1_ID,h,current_sample_length,Filter_coeff,matched_adj1_ID);
//
//        //truncate initial samples and end values from convolution process
//        array_resize(matched_adj1_ID,Filter_coeff-1,current_sample_length);
//
//
//       //Length of sample is sample_length-(Filter_coeff-1);
//        current_sample_length=current_sample_length-(Filter_coeff-1);
//
//
//        ber=twoD_equalize_oneD_SOVA(a_merge,current_sample_length,matched_main,matched_adj1_OD,matched_adj1_ID,no_equalizer_coeff,no_equalizer_coeff,no_equalizer_coeff,no_of_target_vals,a_det,SNR);
//
//
//        free (matched_adj1_ID);
//        free(noised_adj1_ID);
//        free (matched_adj1_OD);
//       free(noised_adj1_OD);
//
//
//    }
//    free (matched_main);
//    free(noised_main);
//
//    current_sample_length=current_sample_length+(Filter_coeff-1);
//    strcpy(ber_filename,"ber_SMR_1DEQ_11tap_1DVA_ITI_20_fixed.txt");
//    //strcpy(ber_filename,"ber_SMR_2DEQ_7tap_1DVA_ITI_20__y_norm_fixed.txt");
//    //filestore_ber(SNR,ber,ber_filename);
//    filestore_ber(ber,ber_filename);
//
//}
//    //Filter the noisy signal using the equalizer /coefficients
//
//
//    //free (ber_filename); //Doesn't work??!!!
//    if (oneD_EQ==1)
//      {
//          free (h);
//          free (a);
//      }
//    else if (twoD_EQ==1)
//            {
//                matrix_free (a_merge,3);
//                free (h);
//                free (h_n1);
//                free (h_1);
//
//                delete [] fname_input_adj1_OD;
//                delete [] fname_input_adj1_ID;
//            }
//
//
//
//    free(a_det);
//
//    delete [] ber_filename;
//    delete [] fname_input_main;
//    delete [] current_wd;
//
//return OK;
//
//}
