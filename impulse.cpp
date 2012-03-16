#include "init.h"
#include "impulse.h"
#include "math_func.h"
#include <math.h>

int smr_impulse_response_1D(long double *h,double PW_50_a,double Tx_a, double Amp,float TMR_offset,int number_of_values)
{
    register int i;
    int start_val, bit_position;
    long double h_norm;
    double PW_n;


    PW_n=Tx_a/PW_50_a;          //Normalized response
    start_val=(number_of_values-1)/2; //Works for both odd and even
    start_val=-start_val;
    bit_position=start_val;

    Amp=exp(Response_constant*Response_constant*(-0.5)*pow(TMR_offset,2)/pow(P_cross,2));
    for (i=0;i<number_of_values;i++)
       {
        h[i]=Amp*exp(Response_constant*Response_constant*(-0.5)*pow(bit_position,2)*pow(PW_n,2));;
        bit_position++;
       }

    h_norm=norm(h,Filter_coeff);
    array_divide(h,h_norm,Filter_coeff);//Compute the normalized filter response
    array_multiply(h,Amp,Filter_coeff);
    return OK;
}


int smr_impulse_response_2D (long double *h_n1,long double *h,long double *h_1,double PW_50_a,double Tx_a, long double PW50_c_n1,long double Ty_c_n1,long double PW50_c_1,long double Ty_c_1, double Amp,int number_of_values)
{
    register int i;
    int start_val, bit_position;
    long double h_norm;
    double PW_n;
    long double As_c_n1, As_c_1;


    PW_n=Tx_a/PW_50_a; //Normalized main track

    start_val=(number_of_values-1)/2;
    start_val=-start_val;


    bit_position=start_val;
    for (i=0;i<number_of_values;i++)
       {
        h[i]=Amp*exp(Response_constant*Response_constant*(-0.5)*pow(bit_position,2)*pow(PW_n,2));
        bit_position++;
       }

    h_norm=norm(h,Filter_coeff);

    array_divide(h,h_norm,Filter_coeff);//Compute the normalized filter response

    array_copy(h_n1,h,Filter_coeff);      //Same as before
    array_copy(h_1,h,Filter_coeff);      //Same as before
    //Cannot compute and normalize as is since it is only As x the normalized response that is filtered by sidetracks
    As_c_n1=exp(Response_constant*Response_constant*(-0.5)*pow((Ty_c_n1/PW50_c_n1),2));
    As_c_1=exp(Response_constant*Response_constant*(-0.5)*pow((Ty_c_1/PW50_c_1),2));  //Side track amplitude for side track response
    array_multiply(h_n1,As_c_n1,Filter_coeff);
    array_multiply(h_1,As_c_1,Filter_coeff);
    return OK;
}

int main_impulse_response_parameter_init(long double *main_PW50_cross_n1,long double *main_Ty_n1,long double *main_PW50_cross_1,long double *main_Ty_1,float TMR_Main,long double main_track_pitch,long double *main_h_n1,long double *main_h,long double *main_h_1)
{


     if (ITI==20)
    {
          /***********Main Track******************************************/

        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        //0.1984
        *main_PW50_cross_n1=4.2; //OD track
        *main_Ty_n1=5.6;
        *main_PW50_cross_1=4.2; //ID track
        *main_Ty_1=5.6;
        //0.21
        /********************Mine********************************************/
  }
    else if (ITI==25)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=4.5; //OD track
        *main_Ty_n1=5.6;
        *main_PW50_cross_1=4.5; //ID track
        *main_Ty_1=5.6;
        //0.2444
        /********************Mine********************************************/
    }
    else if (ITI==27)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=4.8; //OD track
        *main_Ty_n1=5.6;
        *main_PW50_cross_1=4.8; //ID track
        *main_Ty_1=5.6;
        //0.2898
        /********************Mine********************************************/
    }
    else if (ITI==15)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=3.9; //OD track
        *main_Ty_n1=5.6;
        *main_PW50_cross_1=3.9; //ID track
        *main_Ty_1=5.6;
        //0.1532
        /********************Mine********************************************/
    }
    else if (ITI==10)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=3.5; //OD track
        *main_Ty_n1=5.6;
        *main_PW50_cross_1=3.5; //ID track
        *main_Ty_1=5.6;
        //0.0974
        /********************Mine********************************************/
    }
    else if (ITI==5)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *main_PW50_cross_n1=3.1; //OD track
        *main_Ty_n1=5.6;
        *main_PW50_cross_1=3.1; //ID track
        *main_Ty_1=5.6;
        //0.0513
        /********************Mine********************************************/
    }



    //tmr is considered towards ID side
    tmr_addition(TMR_Main,main_Ty_n1,main_Ty_1,main_track_pitch,*main_Ty_n1,*main_Ty_1);


    smr_impulse_response_2D(main_h_n1,main_h,main_h_1,PW50_along,Tx,*main_PW50_cross_n1,*main_Ty_n1,*main_PW50_cross_1,*main_Ty_1,Amplitude,Filter_coeff);
    /*******************************************************************/

    return OK;

}

int od_impulse_response_parameter_init(long double *od_PW50_cross_n1,long double *od_Ty_n1,long double *od_PW50_cross_1,long double *od_Ty_1,float TMR_OD,long double od_track_pitch,long double *od_h_n1,long double *od_h,long double *od_h_1)
{

        if (ITI==20)
    {
          /***********Main Track******************************************/

        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        //0.1984
        *od_PW50_cross_n1=4.2; //OD track
        *od_Ty_n1=5.6;
        *od_PW50_cross_1=4.2; //ID track
        *od_Ty_1=5.6;
        //0.21
        /********************Mine********************************************/
  }
    else if (ITI==25)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *od_PW50_cross_n1=4.5; //OD track
        *od_Ty_n1=5.6;
        *od_PW50_cross_1=4.5; //ID track
        *od_Ty_1=5.6;
        //0.2444
        /********************Mine********************************************/
    }
    else if (ITI==27)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *od_PW50_cross_n1=4.8; //OD track
        *od_Ty_n1=5.6;
        *od_PW50_cross_1=4.8; //ID track
        *od_Ty_1=5.6;
        //0.2898
        /********************Mine********************************************/
    }
    else if (ITI==15)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *od_PW50_cross_n1=3.9; //OD track
        *od_Ty_n1=5.6;
        *od_PW50_cross_1=3.9; //ID track
        *od_Ty_1=5.6;
        //0.1532
        /********************Mine********************************************/
    }
    else if (ITI==10)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *od_PW50_cross_n1=3.5; //OD track
        *od_Ty_n1=5.6;
        *od_PW50_cross_1=3.5; //ID track
        *od_Ty_1=5.6;
        //0.0974
        /********************Mine********************************************/
    }
    else if (ITI==5)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *od_PW50_cross_n1=3.1; //OD track
        *od_Ty_n1=5.6;
        *od_PW50_cross_1=3.1; //ID track
        *od_Ty_1=5.6;
        //0.0513
        /********************Mine********************************************/
    }



    //tmr is considered towards ID side
    tmr_addition(TMR_OD,od_Ty_n1,od_Ty_1,od_track_pitch,*od_Ty_n1,*od_Ty_1);


    /***********Cross track PW 50 and bit width definition**********/

    smr_impulse_response_2D(od_h_n1,od_h,od_h_1,PW50_along,Tx,*od_PW50_cross_n1,*od_Ty_n1,*od_PW50_cross_1,*od_Ty_1,Amplitude,Filter_coeff);
    /*******************************************************************/

return OK;

}

int id_impulse_response_parameter_init(long double *id_PW50_cross_n1,long double *id_Ty_n1,long double *id_PW50_cross_1,long double *id_Ty_1,float TMR_ID,long double id_track_pitch,long double *id_h_n1,long double *id_h,long double *id_h_1)
{

         if (ITI==20)
    {
          /***********Main Track******************************************/

        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        //0.1984
        *id_PW50_cross_n1=4.2; //OD track
        *id_Ty_n1=5.6;
        *id_PW50_cross_1=4.2; //ID track
        *id_Ty_1=5.6;
        //0.21
        /********************Mine********************************************/
  }
    else if (ITI==25)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *id_PW50_cross_n1=4.5; //OD track
        *id_Ty_n1=5.6;
        *id_PW50_cross_1=4.5; //ID track
        *id_Ty_1=5.6;
        //0.2444
        /********************Mine********************************************/
    }
    else if (ITI==27)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *id_PW50_cross_n1=4.8; //OD track
        *id_Ty_n1=5.6;
        *id_PW50_cross_1=4.8; //ID track
        *id_Ty_1=5.6;
        //0.2898
        /********************Mine********************************************/
    }
    else if (ITI==15)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *id_PW50_cross_n1=3.9; //OD track
        *id_Ty_n1=5.6;
        *id_PW50_cross_1=3.9; //ID track
        *id_Ty_1=5.6;
        //0.1532
        /********************Mine********************************************/
    }
    else if (ITI==10)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *id_PW50_cross_n1=3.5; //OD track
        *id_Ty_n1=5.6;
        *id_PW50_cross_1=3.5; //ID track
        *id_Ty_1=5.6;
        //0.0974
        /********************Mine********************************************/
    }
    else if (ITI==5)
    {
        /***********Cross track PW 50 and bit width deifnition**********/
        //Define PW50_cross and Ty according to required ITI
        *id_PW50_cross_n1=3.1; //OD track
        *id_Ty_n1=5.6;
        *id_PW50_cross_1=3.1; //ID track
        *id_Ty_1=5.6;
        //0.0513
        /********************Mine********************************************/
    }





    //tmr is considered towards ID side
    tmr_addition(TMR_ID,id_Ty_n1,id_Ty_1,id_track_pitch,*id_Ty_n1,*id_Ty_1);




    smr_impulse_response_2D(id_h_n1,id_h,id_h_1,PW50_along,Tx,*id_PW50_cross_n1,*id_Ty_n1,*id_PW50_cross_1,*id_Ty_1,Amplitude,Filter_coeff);
    /*******************************************************************/

    return OK;

}



// Offtrack is towards ID side tracks i.e. read less from OD and more from ID
    /********************TMR Addition for Mine specs************************/
            //Offtrack moves towards OD track
            //5% TMR
            //Ty_n1=5.3;
            //Ty_1=5.9;

            //10% TMR
            //Ty_n1=5.0;
            //Ty_1=6.2;

            //15% TMR
            //Ty_n1=4.8;
            //Ty_1=6.4;

            //20% TMR
            //*Ty_n1=4.48;
            //*Ty_1=6.72;

            //For 1-sided look inside if statement
    /********************TMR Addition for Mine specs************************/


        /********************TMR Addition for One Sided specs************************/
            //5% TMR
            //Ty_n1=5.25;
            //Ty_1=5.95;

            //10% TMR
            //Ty_n1=5.5;
            //Ty_1=5.7;

            //15% TMR
            //Ty_n1=5.75;
            //Ty_1=5.45;

            //20% TMR
            //*Ty_n1=6.0;
            //*Ty_1=5.2;

    /********************TMR Addition for Mine specs************************/


int tmr_addition(float percentage_tmr,long double *Ty_n1,long double *Ty_1,long double track_pitch_main,long double track_pitch_OD,long double track_pitch_ID)
{
    long double val;

    val=percentage_tmr/100;

    *Ty_n1=track_pitch_OD + (val*track_pitch_main);
    *Ty_1=track_pitch_ID - (val*track_pitch_main);

    return OK;
}



