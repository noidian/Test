#ifndef IMPULSE
#define IMPULSE

#define Filter_coeff 7      //Defines matched filter for impulse response coefficients

/******************SMR***************************************/
//Along track values : _a denotes along track
#define PW50_along 5.2
#define P_cross 4.8 //ITI 20 =4.2 ITI 27= 4.8
#define Tx 5.2
#define Amplitude 1
#define Response_constant 1.34898   //Hwang paper
//Cross track values can vary depending on ITI
/******************SMR***************************************/




int impulse_response_init (long double *h,double PW_50_a,double Tx_a, double PW_50_c, double Ty_c, double Amp,int number_of_values);


int smr_impulse_response_1D(long double *h,double PW_50_a,double Tx_a, double Amp,float TMR_offset,int number_of_values);


int smr_impulse_response_2D (long double *h_n1,long double *h,long double *h_1,double PW_50_a,double Tx_a, long double PW50_c_n1,long double Ty_n1,long double PW50_c_1,long double Ty_1, double Amp,int number_of_values);

int main_impulse_response_parameter_init(long double *main_PW50_cross_n1,long double *main_Ty_n1,long double *main_PW50_cross_1,long double *main_Ty_1,float TMR_Main,long double main_track_pitch,long double *main_h_n1,long double *main_h,long double *main_h_1);

int od_impulse_response_parameter_init(long double *od_PW50_cross_n1,long double *od_Ty_n1,long double *od_PW50_cross_1,long double *od_Ty_1,float TMR_OD,long double od_track_pitch,long double *od_h_n1,long double *od_h,long double *od_h_1);

int id_impulse_response_parameter_init(long double *id_PW50_cross_n1,long double *id_Ty_n1,long double *id_PW50_cross_1,long double *id_Ty_1,float TMR_ID,long double id_track_pitch,long double *id_h_n1,long double *id_h,long double *id_h_1);

int tmr_addition(float percentage_tmr,long double *Ty_n1,long double *Ty_1,long double track_pitch_main,long double track_pitch_OD,long double track_pitch_ID);


#endif
